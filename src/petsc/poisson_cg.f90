!!! PETSc-based implementation of the iterative Poisson solver

module m_cg_types
  !! Types module providing the context type required by the PETSc matrix-free
  !! operator.

  implicit none

  private

  type, public :: mat_ctx
  end type mat_ctx

end module m_cg_types

submodule(m_poisson_cg) m_petsc_poisson_cg
  !! Module implementing a Poisson solver based on the (preconditioned)
  !! Conjugate Gradient method using PETSc.

  use petsc

  use m_cg_types

  implicit none

  type, extends(field_t) :: petsc_field_t
    !! Field extension to wrap PETSc vector data in a field interface
  end type petsc_field_t

  type, extends(poisson_cg_t) :: petsc_poisson_cg_t
    !! Conjugate Gradient based Poisson solver using PETSc as a backend.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.

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
      type(mat_ctx) :: ctx ! The shell matrix context
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
      type(mat_ctx) :: ctx ! The shell matrix context
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
      type(mat_ctx) :: ctx ! The shell matrix context
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

  module function init_cg(xdirps, ydirps, zdirps) result(poisson_cg)
    !! Public constructor for the poisson_cg_t type.
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps ! X/Y/Z discretisation operators
    class(poisson_cg_t), allocatable :: poisson_cg

    allocate (petsc_poisson_cg_t :: poisson_cg)
    select type (poisson_cg)
    type is (petsc_poisson_cg_t)
      call init_petsc_cg(poisson_cg, xdirps, ydirps, zdirps)
    class default
      ! This should be impossible
      print *, "Failure in allocating PETSc Poisson solver -- this indicates a serious problem"
      stop 1
    end select
  end function init_cg

  subroutine init_petsc_cg(self, xdirps, ydirps, zdirps)
    !! Private constructor for the poisson_cg_t type.
    type(petsc_poisson_cg_t), intent(inout) :: self
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps ! X/Y/Z discretisation operators

    integer :: nx, ny, nz, n ! Local problem size

    ! Determine local problem size
    nx = xdirps%n
    ny = ydirps%n
    nz = zdirps%n
    n = nx*ny*nz

    ! Initialise preconditioner and operator matrices
    ! XXX: Add option to use preconditioner as operator (would imply low-order
    !      solution)?
    call create_matrix(n, "assemled", self%P)
    call create_matrix(n, "matfree", self%A)
  end subroutine init_petsc_cg

  subroutine create_matrix(nlocal, mat_type, M)
    !! Creates either a matrix object given the local problem size.
    !! The matrix can be either "assembled" - suitable for preconditioners, or
    !! "matfree" - for use as a high-order operator.
    integer, intent(in) :: nlocal            ! The local problem size
    character(len=*), intent(in) :: mat_type ! The desired type of matrix - valid values
                                             ! are "assembled" or "matfree"
    type(tMat), intent(out) :: M             ! The matrix object

    type(mat_ctx) :: ctx

    integer :: ierr

    if (mat_type == "assembled") then
      call MatCreate(PETSC_COMM_WORLD, M, ierr)
      call MatSetSizes(M, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE, ierr)
      call MatSetFromOptions(M, ierr)
    else
      ! TODO: How do we get the ctx?

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

    type(mat_ctx) :: ctx

    type(petsc_field_t), allocatable :: xfield ! Field wrapper for the input vector
    type(petsc_field_t), allocatable :: ffield ! Field wrapper for the output vector

    call MatShellGetContext(M, ctx, ierr)

    ! Wrap x and f in a field
    call poissmult(xfield, ffield)
    ! Extract field into f

  end subroutine poissmult_petsc

end submodule m_petsc_poisson_cg
