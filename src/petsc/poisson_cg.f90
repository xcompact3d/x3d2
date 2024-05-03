!!! PETSc-based implementation of the iterative Poisson solver

module m_cg_types
  !! Types module providing the context type required by the PETSc matrix-free
  !! operator.

  use m_base_backend, only: base_backend_t
  use m_allocator, only: field_t
  use m_poisson_cg, only: laplace_operator_t
  use m_tdsops, only: dirps_t

  implicit none

  private
  public :: mat_ctx_t

  type, public :: mat_ctx_t
    class(base_backend_t), pointer :: backend
    class(laplace_operator_t), pointer :: lapl
    class(field_t), pointer :: xfield
    class(field_t), pointer :: ffield
    class(dirps_t), pointer :: xdirps, ydirps, zdirps
  end type mat_ctx_t

  interface mat_ctx_t
    module procedure init_ctx
  end interface mat_ctx_t

contains

  function init_ctx(backend, lapl, dir) result(ctx)

    class(base_backend_t), pointer, intent(in) :: backend
    class(laplace_operator_t), pointer, intent(in) :: lapl
    integer, intent(in) :: dir
    type(mat_ctx_t) :: ctx

    ctx%xfield => backend%allocator%get_block(dir)
    ctx%ffield => backend%allocator%get_block(dir)
    ctx%backend => backend
    ctx%lapl => lapl

  end function init_ctx

end module m_cg_types

submodule(m_poisson_cg) m_petsc_poisson_cg
  !! Module implementing a Poisson solver based on the (preconditioned)
  !! Conjugate Gradient method using PETSc.

  use petsc

  use m_cg_types

  use m_common, only: dp, DIR_X, DIR_C
  use m_base_backend, only: base_backend_t

  implicit none

  type, extends(field_t) :: petsc_field_t
    !! Field extension to wrap PETSc vector data in a field interface
  end type petsc_field_t

  type, extends(poisson_cg_t) :: petsc_poisson_cg_t
    !! Conjugate Gradient based Poisson solver using PETSc as a backend.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.

    type(mat_ctx_t) :: ctx
    type(tMat) :: A ! The operator matrix
    type(tMat) :: P ! The preconditioner matrix

  contains
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

contains

  module function init_cg(xdirps, ydirps, zdirps, backend) &
                  result(poisson_cg)
    !! Public constructor for the poisson_cg_t type.
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps ! X/Y/Z discretisation operators
    class(poisson_cg_t), allocatable :: poisson_cg
    class(base_backend_t), pointer, intent(in) :: backend

    allocate (petsc_poisson_cg_t :: poisson_cg)

    call init_lapl(poisson_cg%lapl)
    
    select type (poisson_cg)
    type is (petsc_poisson_cg_t)
      call init_petsc_cg(poisson_cg, xdirps, ydirps, zdirps, backend)
    class default
      ! This should be impossible
      print *, "Failure in allocating PETSc Poisson solver -- this indicates a serious problem"
      stop 1
    end select
  end function init_cg

  subroutine init_petsc_cg(self, xdirps, ydirps, zdirps, backend)
    !! Private constructor for the poisson_cg_t type.
    type(petsc_poisson_cg_t), intent(inout) :: self
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps ! X/Y/Z discretisation operators
    class(base_backend_t), pointer, intent(in) :: backend

    integer :: nx, ny, nz, n ! Local problem size

    ! Determine local problem size
    nx = xdirps%n
    ny = ydirps%n
    nz = zdirps%n
    n = nx*ny*nz

    self%ctx = mat_ctx_t(backend, self%lapl, DIR_X)

    ! Initialise preconditioner and operator matrices
    ! XXX: Add option to use preconditioner as operator (would imply low-order
    !      solution)?
    call create_matrix(n, "assemled", backend, self%lapl, self%P)
    call create_matrix(n, "matfree", backend, self%lapl, self%A)
  end subroutine init_petsc_cg

  subroutine create_matrix(nlocal, mat_type, backend, lapl, M)
    !! Creates either a matrix object given the local problem size.
    !! The matrix can be either "assembled" - suitable for preconditioners, or
    !! "matfree" - for use as a high-order operator.
    integer, intent(in) :: nlocal            ! The local problem size
    character(len=*), intent(in) :: mat_type ! The desired type of matrix - valid values
                                             ! are "assembled" or "matfree"
    class(base_backend_t), pointer, intent(in) :: backend  ! The compute backend
    class(laplace_operator_t), pointer, intent(in) :: lapl ! The Laplacian operator
    type(tMat), intent(out) :: M             ! The matrix object

    type(mat_ctx_t) :: ctx

    integer :: ierr

    if (mat_type == "assembled") then
      call MatCreate(PETSC_COMM_WORLD, M, ierr)
      call MatSetSizes(M, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE, ierr)
      call MatSetFromOptions(M, ierr)
    else
      call MatCreateShell(PETSC_COMM_WORLD, nlocal, nlocal, PETSC_DETERMINE, &
                          PETSC_DETERMINE, ctx, M, ierr)
      call MatShellSetContext(M, ctx, ierr) ! Is this necessary?
      call MatShellSetOperation(M, MATOP_MULT, poissmult_petsc, ierr)
    end if
    call MatSetUp(M, ierr)

  end subroutine create_matrix

  subroutine poissmult_petsc(M, x, f, ierr)
    !! Computes the action of the Poisson operator, i.e. `f = Mx` where `M` is
    !! the discrete Laplacian.
    type(tMat) :: M ! The operator
    type(tVec) :: x ! The input vector
    type(tVec) :: f ! The output vector
    integer :: ierr ! The error code

    type(mat_ctx_t) :: ctx

    call MatShellGetContext(M, ctx, ierr)

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

    integer :: nx, ny, nz

    nx = backend%nx_loc
    ny = backend%ny_loc
    nz = backend%nz_loc

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

    integer :: nx, ny, nz
    
    nx = backend%nx_loc
    ny = backend%ny_loc
    nz = backend%nz_loc

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
