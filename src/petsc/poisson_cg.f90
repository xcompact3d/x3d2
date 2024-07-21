!!! PETSc-based implementation of the iterative Poisson solver

module m_cg_types
  !! Types module providing the context type required by the PETSc matrix-free
  !! operator.

  use m_base_backend, only: base_backend_t
  use m_allocator, only: field_t
  use m_poisson_cg, only: laplace_operator_t

  implicit none

  private

  type, public :: mat_ctx_t
    class(base_backend_t), pointer :: backend
    type(laplace_operator_t) :: lapl
    class(field_t), pointer :: xfield
    class(field_t), pointer :: ffield
  end type mat_ctx_t

  interface mat_ctx_t
    module procedure init_ctx
  end interface mat_ctx_t

contains

  function init_ctx(backend, lapl, dir) result(ctx)

    class(base_backend_t), target, intent(in) :: backend
    class(laplace_operator_t), intent(in) :: lapl
    integer, intent(in) :: dir
    type(mat_ctx_t) :: ctx

    ctx%xfield => backend%allocator%get_block(dir)
    ctx%ffield => backend%allocator%get_block(dir)
    ctx%backend => backend
    ctx%lapl = lapl

  end function init_ctx

end module m_cg_types

submodule(m_poisson_cg) m_petsc_poisson_cg
  !! Module implementing a Poisson solver based on the (preconditioned)
  !! Conjugate Gradient method using PETSc.

  use petsc

  use m_cg_types

  use m_common, only: dp, DIR_X, DIR_C, CELL
  use m_base_backend, only: base_backend_t

  implicit none

  type, extends(field_t) :: petsc_field_t
    !! Field extension to wrap PETSc vector data in a field interface
  end type petsc_field_t

  type, extends(poisson_solver_t) :: petsc_poisson_cg_t
    !! Conjugate Gradient based Poisson solver using PETSc as a backend.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.

    type(mat_ctx_t) :: ctx
    type(tKSP) :: ksp  ! The solver
    type(tMat) :: Amat ! The operator matrix
    type(tMat) :: Pmat ! The preconditioner matrix
    type(tVec) :: fvec ! The RHS vector
    type(tVec) :: pvec ! The solution vector

  contains
    procedure :: solve => solve_petsc
    procedure :: create_operator
    procedure :: create_preconditioner
    procedure :: create_vectors
    procedure :: create_solver
  end type petsc_poisson_cg_t

  interface MatCreateShell
    !! Defines the interface to the external (PETSc) function to create a
    !! matrix-free operator.
    subroutine MatCreateShell(comm, nrow_l, ncol_l, nrow_g, ncol_g, ctx, M, &
                              ierr)
      use petsc
      use m_cg_types
      integer :: comm
      integer :: nrow_l ! Local number of rows
      integer :: ncol_l ! Local number of columns
      integer :: nrow_g ! Global number of rows
      integer :: ncol_g ! Global number of columns
      type(mat_ctx_t) :: ctx ! The shell matrix context
      type(tMat) :: M   ! The matrix object
      integer :: ierr
    end subroutine MatCreateShell
  end interface MatCreateShell

  interface MatShellSetContext
    !! Defines the interface to the external (PETSc) function to store
    !! application-dependent information required by the matrix-free operator.
    subroutine MatShellSetContext(M, ctx, ierr)
      use petsc
      use m_cg_types
      type(tMat) :: M      ! The matrix object
      type(mat_ctx_t) :: ctx ! The shell matrix context
      integer :: ierr
    end subroutine MatShellSetContext
  end interface MatShellSetContext

  interface MatShellGetContext
    !! Defines the interface to the external (PETSc) function to retrieve
    !! application-dependent information from the matrix-free operator.
    subroutine MatShellGetContext(M, ctx, ierr)
      use petsc
      use m_cg_types
      type(tMat) :: M      ! The matrix object
      type(mat_ctx_t) :: ctx ! The shell matrix context
      integer :: ierr
    end subroutine MatShellGetContext
  end interface MatShellGetContext

  interface MatShellSetOperation
    !! Defines the interface to the external (PETSc) function to set the
    !! matrix-free operator procedure that evaluates `f = Mx`.
    subroutine MatShellSetOperation(M, OP, fn, ierr)
      use petsc
      type(tMat) :: M
      integer :: OP
      interface
        subroutine fn(M, x, f, ierr)
          use petsc
          type(tMat) :: M ! The operator
          type(tVec) :: x ! The input vector
          type(tVec) :: f ! The output vector
          integer :: ierr ! The error code
        end subroutine fn
      end interface
      integer :: ierr
    end subroutine MatShellSetOperation
  end interface MatShellSetOperation

  type(mat_ctx_t), save :: ctx_global ! XXX: This sucks!
  
contains

  module subroutine solve_petsc(self, p, f, backend)
    class(petsc_poisson_cg_t) :: self
    class(field_t), intent(inout) :: p ! Pressure solution
    class(field_t), intent(in) :: f    ! Poisson RHS
    class(base_backend_t), intent(in) :: backend

    integer :: ierr

    call copy_field_to_vec(self%fvec, f, backend)
    call KSPSolve(self%ksp, self%fvec, self%pvec, ierr)
    call copy_vec_to_field(p, self%pvec, backend)
    
  end subroutine solve_petsc

  module subroutine init_solver(solver, backend) 
    !! Public constructor for the poisson_cg_t type.
    class(poisson_solver_t), allocatable, intent(out) :: solver
    class(base_backend_t), target, intent(in) :: backend

    allocate (petsc_poisson_cg_t :: solver)
    
    select type (solver)
    type is (petsc_poisson_cg_t)
      call init_petsc_cg(solver, backend)
    class default
      ! This should be impossible
      error stop "Failure in allocating PETSc Poisson solver -- this indicates a serious problem"
    end select
  end subroutine init_solver

  subroutine init_petsc_cg(self, backend)
    !! Private constructor for the poisson_cg_t type.
    type(petsc_poisson_cg_t), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend

    integer :: n ! Local problem size

    integer :: ierr

    call PetscInitialize(ierr)

    ! Determine local problem size
    n = product(backend%mesh%get_dims(CELL))

    ctx_global = mat_ctx_t(backend, self%lapl, DIR_X)

    ! Initialise preconditioner and operator matrices
    ! XXX: Add option to use preconditioner as operator (would imply low-order
    !      solution)?
    call self%create_operator(n)
    call self%create_preconditioner(n)

    ! Initialise RHS and solution vectors
    call self%create_vectors(n)

    ! Create the linear system
    call self%create_solver()

  end subroutine init_petsc_cg

  subroutine create_operator(self, n)
    class(petsc_poisson_cg_t) :: self
    integer, intent(in) :: n

    integer :: ierr

    call MatCreateShell(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, &
      PETSC_DETERMINE, self%ctx, self%Amat, ierr)
    call MatShellSetContext(self%Amat, self%ctx, ierr) ! Is this necessary?
    call MatShellSetOperation(self%Amat, MATOP_MULT, poissmult_petsc, ierr)
    call MatSetUp(self%Amat, ierr)

    call MatAssemblyBegin(self%Amat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Amat, MAT_FINAL_ASSEMBLY, ierr)
    
  end subroutine create_operator

  subroutine create_preconditioner(self, n)
    class(petsc_poisson_cg_t) :: self
    integer :: n

    integer :: ierr

    type(tVec) :: v

    call MatCreate(PETSC_COMM_WORLD, self%Pmat, ierr)
    call MatSetSizes(self%Pmat, n, n, PETSC_DECIDE, PETSC_DECIDE, ierr)
    call MatSetFromOptions(self%Pmat, ierr)
    call MatSetUp(self%Pmat, ierr)

    !! Create an identity matrix
    call create_vec(v, n)
    call VecSet(v, 1.0_dp, ierr)
    call MatDiagonalSet(self%Pmat, v, INSERT_VALUES, ierr)
    
    call MatAssemblyBegin(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine create_preconditioner

  subroutine create_vectors(self, n)
    class(petsc_poisson_cg_t) :: self
    integer, intent(in) :: n

    call create_vec(self%fvec, n)
    call create_vec(self%pvec, n)
    
  end subroutine create_vectors

  subroutine create_vec(v, n)
    type(tVec), intent(out) :: v
    integer, intent(in) :: n

    integer :: ierr

    call VecCreate(PETSC_COMM_WORLD, v, ierr)
    call VecSetSizes(v, n, PETSC_DETERMINE, ierr)
    call VecSetFromOptions(v, ierr)

    call VecAssemblyBegin(v, ierr)
    call VecAssemblyEnd(v, ierr)
    
  end subroutine create_vec

  subroutine create_solver(self)
    class(petsc_poisson_cg_t) :: self

    integer :: ierr

    call KSPCreate(PETSC_COMM_WORLD, self%ksp, ierr)
    call KSPSetOperators(self%ksp, self%Amat, self%Pmat, ierr)
    call KSPSetFromOptions(self%ksp, ierr)

  end subroutine create_solver

  subroutine poissmult_petsc(M, x, f, ierr)
    !! Computes the action of the Poisson operator, i.e. `f = Mx` where `M` is
    !! the discrete Laplacian.
    type(tMat) :: M ! The operator
    type(tVec) :: x ! The input vector
    type(tVec) :: f ! The output vector
    integer :: ierr ! The error code

    type(mat_ctx_t) :: ctx

    ! XXX: Fixme
    ! call MatShellGetContext(M, ctx, ierr)
    ! print *, ctx%foo
    ctx = ctx_global

    call copy_vec_to_field(ctx%xfield, x, ctx%backend)
    call ctx%lapl%apply(ctx%ffield, ctx%xfield, ctx%backend)
    call copy_field_to_vec(f, ctx%ffield, ctx%backend)

  end subroutine poissmult_petsc

  subroutine copy_vec_to_field(f, v, backend)
    !! Copies the contents of a PETSc vector into an x3d2 field object
    ! XXX: This can be avoided if a field can wrap the vector memory
    class(field_t), intent(inout) :: f
    type(tVec) :: v
    class(base_backend_t), intent(in) :: backend

    real(dp), dimension(:), pointer :: vdata
    real(dp), dimension(:, :, :), pointer :: vdata3d
    integer :: ierr

    integer, dimension(3) :: dims
    integer :: nx, ny, nz

    dims = backend%mesh%get_dims(CELL)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    ! Local copy
    call VecGetArrayReadF90(v, vdata, ierr)
    if (nx*ny*nz /= size(vdata)) then
      print *, "Vector and field sizes are incompatible (padding?)"
      stop 1
    end if
    vdata3d(1:nx, 1:ny, 1:nz) => vdata(:) ! Get a 3D representation of the vector
    call backend%set_field_data(f, vdata3d, DIR_C)
    call VecRestoreArrayReadF90(v, vdata, ierr)
    nullify (vdata3d)

    ! Halo exchange

  end subroutine copy_vec_to_field

  subroutine copy_field_to_vec(v, f, backend)
    !! Copies the contents of an x3d2 field object into a PETSc vector
    ! XXX: This can be avoided if a field can wrap the vector memory
    type(tVec) :: v
    class(field_t), intent(in) :: f
    class(base_backend_t), intent(in) :: backend

    real(dp), dimension(:), pointer :: vdata
    real(dp), dimension(:, :, :), pointer :: vdata3d
    integer :: ierr

    integer, dimension(3) :: dims
    integer :: nx, ny, nz
    
    dims = backend%mesh%get_dims(CELL)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    ! Local copy
    call VecGetArrayF90(v, vdata, ierr)
    if (nx*ny*nz /= size(vdata)) then
      print *, "Vector and field sizes are incompatible (padding?)"
      stop 1
    end if
    vdata3d(1:nx, 1:ny, 1:nz) => vdata(:) ! Get a 3D representation of the vector
    call backend%get_field_data(vdata3d, f, DIR_C)
    call VecRestoreArrayF90(v, vdata, ierr)
    nullify (vdata3d)

    ! Halo exchange
    call VecAssemblyBegin(v, ierr)
    call VecAssemblyEnd(v, ierr)

  end subroutine copy_field_to_vec

end submodule m_petsc_poisson_cg
