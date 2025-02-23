!!! PETSc-based implementation of the iterative Poisson solver

module m_cg_types
  !! Types module providing the context type required by the PETSc matrix-free
  !! operator.

  use m_common, only: CELL
  use m_base_backend, only: base_backend_t
  use m_allocator, only: field_t
  use m_base_poisson_cg, only: laplace_operator_t

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

    ctx%xfield => backend%allocator%get_block(dir, CELL)
    ctx%ffield => backend%allocator%get_block(dir, CELL)
    ctx%backend => backend
    ctx%lapl = lapl

  end function init_ctx

end module m_cg_types

module m_poisson_cg_backend
  !! Module implementing a Poisson solver based on the (preconditioned)
  !! Conjugate Gradient method using PETSc.

  use m_base_poisson_cg, only: poisson_solver_t, poisson_precon_t, laplace_operator_t

  use petsc

  use m_cg_types

  use m_common, only: dp, DIR_Z, DIR_C, CELL
  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: init_solver
  public :: init_precon

  type, extends(poisson_solver_t), public :: poisson_solver_impl
    !! Conjugate Gradient based Poisson solver using PETSc as a backend.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.
    type(mat_ctx_t) :: ctx
    type(tKSP) :: ksp  ! The solver
    type(tMat) :: Amat ! The operator matrix
    type(tVec) :: fvec ! The RHS vector
    type(tVec) :: pvec ! The solution vector
    type(laplace_operator_t) :: lapl
    class(poisson_precon_t), allocatable :: precon
  contains
    procedure, public :: solve => solve_petsc
    procedure :: create_operator
    procedure :: create_vectors
    procedure :: create_solver
  end type poisson_solver_impl

  type, extends(poisson_precon_t), public :: poisson_precon_impl
    !! The PETSc implementation of the Poisson preconditioner, implements a 2nd
    !! order finite difference approximation of the Laplacian.
    type(tMat) :: Pmat ! The preconditioner matrix
    type(tDM) :: da
  contains
    procedure :: init => init_precon_petsc
    procedure :: apply => petsc_apply_precon
  end type poisson_precon_impl

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

  function init_precon(backend) result(precon)
    ! Constructs the PETSc preconditioner implementation
    class(base_backend_t), intent(in) :: backend
    class(poisson_precon_t), allocatable :: precon

    allocate (poisson_precon_impl :: precon)
    select type (precon)
    type is (poisson_precon_impl)
      call precon%init(backend)
    class default
      error stop "IMPOSSIBLE"
    end select

  end function init_precon

  subroutine init_precon_petsc(self, backend)
    ! Initialise the PETSc implementation of the preconditioner object
#include "petsc/finclude/petscmat.h"

    class(poisson_precon_impl), intent(out) :: self
    class(base_backend_t), intent(in) :: backend

    type(tMatNullSpace) :: nsp
    integer :: ierr

    integer, parameter :: nnb = 26 ! Number of neighbours (27-point stencil has 26 neighbours)

    integer, dimension(3) :: dims
    integer :: i, j, k
    real(dp) :: dx, dy, dz

    MatStencil :: row(4, 1)
    MatStencil :: col(4, 27)
    real(dp), dimension(nnb + 1) :: coeffs

    logical :: initialised

    integer, dimension(3), parameter :: stencil1d = [1, -2, 1]
    integer, dimension(3, 3, 3) :: stencil3d
    integer, dimension(3, 3, 3) :: stencil3d_x, stencil3d_y, stencil3d_z

    integer, dimension(:), allocatable :: procx, procy, procz
    integer, dimension(:), allocatable :: nxglobal, nyglobal, nzglobal
    integer, dimension(:), allocatable :: lx, ly, lz
    integer, dimension(:, :, :), allocatable :: procgrid
    integer, parameter :: dof = 1 ! Variables per point in the linear system (P)
    integer, parameter :: stencil_width = 1
    integer :: ctr
    integer :: ii, jj, kk
    integer :: ifirst, jfirst, kfirst
    integer :: ilast, jlast, klast

    ! Ensure PETSc is initialised
    call PetscInitialized(initialised, ierr)
    if (.not. initialised) then
      if (backend%mesh%par%nrank == 0) then
        print *, "Initialising PETSc"
      end if
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    end if
    if (backend%mesh%par%nrank == 0) then
      print *, "PETSc Initialised"
    end if

    ! Create an explicit preconditioner matrix
    associate (mesh => backend%mesh)

      allocate(procgrid(mesh%par%nproc_dir(1), mesh%par%nproc_dir(2), mesh%par%nproc_dir(3)))
      procgrid = mesh%par%compute_global_rank_layout()
      
      procx = procgrid(:, mesh%par%nrank_dir(2) + 1, mesh%par%nrank_dir(3) + 1)
      procy = procgrid(mesh%par%nrank_dir(1) + 1, :, mesh%par%nrank_dir(3) + 1)
      procz = procgrid(mesh%par%nrank_dir(1) + 1, mesh%par%nrank_dir(2) + 1, :)
      deallocate(procgrid)

      dims = mesh%get_dims(CELL)
      allocate(nxglobal(mesh%par%nproc))
      allocate(nyglobal(mesh%par%nproc))
      allocate(nzglobal(mesh%par%nproc))
      nxglobal = 0; nyglobal = 0; nzglobal = 0
      nxglobal(mesh%par%nrank + 1) = dims(1)
      nyglobal(mesh%par%nrank + 1) = dims(2)
      nzglobal(mesh%par%nrank + 1) = dims(3)

      call MPI_Allreduce(MPI_IN_PLACE, nxglobal, mesh%par%nproc, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, nyglobal, mesh%par%nproc, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, nzglobal, mesh%par%nproc, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      
      lx = nxglobal(procx + 1)
      ly = nyglobal(procy + 1)
      lz = nzglobal(procz + 1)

      deallocate(nxglobal)
      deallocate(nyglobal)
      deallocate(nzglobal)
      
      dims = mesh%get_global_dims(CELL)
      call DMDACreate3d(PETSC_COMM_WORLD, &
                        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX, &
                        dims(1), dims(2), dims(3), &
                        mesh%par%nproc_dir(1), mesh%par%nproc_dir(2), mesh%par%nproc_dir(3), &
                        dof, &
                        stencil_width, &
                        lx, ly, lz, &
                        self%da, &
                        ierr)
      call DMSetFromOptions(self%da, ierr)
      call DMSetUp(self%da, ierr)
    end associate

    call DMCreateMatrix(self%da, self%Pmat, ierr)
    call MatSetFromOptions(self%Pmat, ierr)
    call MatSetUp(self%Pmat, ierr)

    ! Set up stencils
    do k = 1, 3
      do j = 1, 3
        stencil3d(:, j, k) = stencil1d
      end do
    end do
    stencil3d(:, 2, 1) = 2 * stencil3d(:, 2, 1)
    stencil3d(:, 1, 2) = 2 * stencil3d(:, 1, 2)
    stencil3d(:, 2, 2) = 4 * stencil3d(:, 2, 2)
    stencil3d(:, 3, 2) = 2 * stencil3d(:, 3, 2)
    stencil3d(:, 2, 3) = 2 * stencil3d(:, 2, 3)

    stencil3d_x = stencil3d
    stencil3d_y = reshape(stencil3d, shape=[3, 3, 3], order=[2, 1, 3])
    stencil3d_z = reshape(stencil3d, shape=[3, 3, 3], order=[3, 2, 1])

    ! Set the Poisson coefficients
    associate (mesh => backend%mesh)
      dims = mesh%get_dims(CELL)
      dx = mesh%geo%d(1); dy = mesh%geo%d(2); dz = mesh%geo%d(3)

      call DMDAGetCorners(self%da, ifirst, jfirst, kfirst, ilast, jlast, klast, ierr)
      ilast = ifirst + (ilast - 1)
      jlast = jfirst + (jlast - 1)
      klast = kfirst + (klast - 1)
      do k = kfirst, klast
        do j = jfirst, jlast
          do i = ifirst, ilast

            row(MatStencil_i, 1) = i
            row(MatStencil_j, 1) = j
            row(MatStencil_k, 1) = k
            ctr = 1
            do kk = -1, 1
              do jj = -1, 1
                do ii = -1, 1
                  col(MatStencil_i, ctr) = i + ii
                  col(MatStencil_j, ctr) = j + jj
                  col(MatStencil_k, ctr) = k + kk
                  ctr = ctr + 1
                end do
              end do
            end do
            coeffs = reshape(stencil3d_x / dx**2 + stencil3d_y / dy**2 + stencil3d_z / dz**2, shape=[27]) / 16.0_dp

            ! Push to matrix
            call MatSetValuesStencil(self%Pmat, 1, row, nnb + 1, col, &
                                     coeffs, INSERT_VALUES, ierr)
          end do
        end do
      end do
    end associate

    call MatAssemblyBegin(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)

    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_VEC, nsp, ierr)
    call MatSetnullSpace(self%Pmat, nsp, ierr)
    call MatNullSpaceDestroy(nsp, ierr)

  end subroutine init_precon_petsc

  module subroutine petsc_apply_precon(self, p, b, backend)
    class(poisson_precon_impl) :: self
    class(field_t), intent(in) :: p
    class(field_t), intent(inout) :: b
    class(base_backend_t), intent(in) :: backend

    type(tVec) :: pVec, bVec
    integer :: ierr

    integer :: n

    n = product(backend%mesh%get_dims(CELL))
    call create_vec(pVec, n)
    call create_vec(bVec, n)

    call copy_field_to_vec(pVec, p, backend)
    call MatMult(self%PMat, pVec, bVec, ierr)
    call copy_vec_to_field(b, bVec, backend)

    call VecDestroy(pVec, ierr)
    call VecDestroy(bVec, ierr)

  end subroutine petsc_apply_precon

  module subroutine solve_petsc(self, p, f, backend)
    class(poisson_solver_impl) :: self
    class(field_t), intent(inout) :: p ! Pressure solution
    class(field_t), intent(in) :: f    ! Poisson RHS
    class(base_backend_t), intent(in) :: backend

    integer :: ierr

    ctx_global = mat_ctx_t(backend, self%lapl, DIR_Z)
    call copy_field_to_vec(self%fvec, f, backend)
    call KSPSolve(self%ksp, self%fvec, self%pvec, ierr)
    call copy_vec_to_field(p, self%pvec, backend)

  end subroutine solve_petsc

  module subroutine init_solver(solver, backend, mesh)
    !! Public constructor for the poisson_cg_t type.
    class(poisson_solver_t), allocatable, intent(out) :: solver
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), intent(in) :: mesh

    allocate (poisson_solver_impl :: solver)

    select type (solver)
    type is (poisson_solver_impl)
      solver%precon = init_precon(backend)
      solver%lapl = laplace_operator_t(backend, mesh)
      call init_petsc_cg(solver, backend)
    class default
      ! This should be impossible
      error stop "Failure in allocating PETSc Poisson solver -- this indicates a serious problem"
    end select
  end subroutine init_solver

  subroutine init_petsc_cg(self, backend)
    !! Private constructor for the poisson_cg_t type.
    type(poisson_solver_impl), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend

    integer :: n ! Local problem size

    integer :: ierr

    logical :: initialised

    call PetscInitialized(initialised, ierr)
    if (.not. initialised) then
      if (backend%mesh%par%nrank == 0) then
        print *, "Initialising PETSc"
      end if
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    end if
    if (backend%mesh%par%nrank == 0) then
      print *, "PETSc Initialised"
    end if

    ! Determine local problem size
    n = product(backend%mesh%get_dims(CELL))

    ! Initialise preconditioner and operator matrices
    ! XXX: Add option to use preconditioner as operator (would imply low-order
    !      solution)?
    call self%create_operator(n)

    ! Initialise RHS and solution vectors
    call self%create_vectors(n)

    ! Create the linear system
    call self%create_solver()

  end subroutine init_petsc_cg

  subroutine create_operator(self, n)
    ! Set the PETSc MATVEC to use the x3d2 high-order Laplacian operator
    class(poisson_solver_impl) :: self
    integer, intent(in) :: n ! The local problem size

    type(tMatNullSpace) :: nsp
    integer :: ierr

    call MatCreateShell(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, &
                        PETSC_DETERMINE, self%ctx, self%Amat, ierr)
    call MatShellSetContext(self%Amat, self%ctx, ierr) ! Is this necessary?
    call MatShellSetOperation(self%Amat, MATOP_MULT, poissmult_petsc, ierr)
    call MatSetUp(self%Amat, ierr)

    call MatAssemblyBegin(self%Amat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Amat, MAT_FINAL_ASSEMBLY, ierr)

    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_VEC, nsp, ierr)
    call MatSetnullSpace(self%Amat, nsp, ierr)
    call MatNullSpaceDestroy(nsp, ierr)

  end subroutine create_operator

  subroutine create_vectors(self, n)
    ! Allocates the pressure and forcing vectors.

    class(poisson_solver_impl) :: self ! The Poisson solver
    integer, intent(in) :: n          ! The local vector size

    call create_vec(self%fvec, n)
    call create_vec(self%pvec, n)

  end subroutine create_vectors

  subroutine create_vec(v, n)
    ! Utility subroutine to allocate a PETSc vector.

    type(tVec), intent(out) :: v ! The vector
    integer, intent(in) :: n     ! The local vector size

    integer :: ierr

    call VecCreate(PETSC_COMM_WORLD, v, ierr)
    call VecSetSizes(v, n, PETSC_DETERMINE, ierr)
    call VecSetFromOptions(v, ierr)

    call VecAssemblyBegin(v, ierr)
    call VecAssemblyEnd(v, ierr)

  end subroutine create_vec

  subroutine create_solver(self)
    ! Sets up the PETSc linear solver.

    class(poisson_solver_impl) :: self

    integer :: ierr

    associate (precon => self%precon)
      select type (precon)
      type is (poisson_precon_impl)
        call KSPCreate(PETSC_COMM_WORLD, self%ksp, ierr)
        call KSPSetOperators(self%ksp, self%Amat, precon%Pmat, ierr)
        call KSPSetFromOptions(self%ksp, ierr)
        call KSPSetInitialGuessNonzero(self%ksp, PETSC_TRUE, ierr)
      class default
        error stop "Poisson preconditioner type is wrong"
      end select
    end associate

  end subroutine create_solver

  subroutine poissmult_petsc(M, x, f, ierr)
    !! Computes the action of the Poisson operator, i.e. `f = Mx` where `M` is
    !! the discrete Laplacian.
    type(tMat) :: M ! The operator
    type(tVec) :: x ! The input vector
    type(tVec) :: f ! The output vector
    integer :: ierr ! The error code

    type(mat_ctx_t) :: ctx

    ierr = 0; 
    ! XXX: Fixme
    ! call MatShellGetContext(M, ctx, ierr)
    associate (matrix => M); end associate ! Silence unused argument
    ! print *, ctx%foo
    ctx = ctx_global

    call copy_vec_to_field(ctx%xfield, x, ctx%backend)
    call ctx%lapl%apply(ctx%ffield, ctx%xfield)
    call copy_field_to_vec(f, ctx%ffield, ctx%backend)

  end subroutine poissmult_petsc

  subroutine copy_vec_to_field(f, v, backend)
    !! Copies the contents of a PETSc vector into an x3d2 field object
    ! XXX: This can be avoided if a field can wrap the vector memory
    class(field_t), intent(inout) :: f ! The destination field
    type(tVec) :: v                    ! The source vector
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
    type(tVec) :: v                 ! The destination vector.
    class(field_t), intent(in) :: f ! The source field.
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

end module m_poisson_cg_backend
