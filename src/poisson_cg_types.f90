module m_base_poisson_cg

  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, &
                      RDR_X2Y, RDR_X2Z, RDR_X2C, &
                      RDR_Y2X, RDR_Z2X, RDR_C2X, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      CELL, &
                      BC_PERIODIC
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t, dirps_t

  implicit none

  private

  type, public :: laplace_operator_t
    !! Operator that computes the Laplacian of a field.
    private
    type(dirps_t) :: xdirps, ydirps, zdirps
  contains
    procedure :: apply => poissmult
    procedure :: poissmult_dirx
  end type laplace_operator_t

  interface laplace_operator_t
    !! Public constructor for the laplace_operator_t type.
    procedure init_lapl
  end interface laplace_operator_t

  type, abstract, public :: poisson_precon_t
    !! Wrapper definition of the Poisson preconditioner.
    !! User code should use this class which will instantiate the
    !! backend-specific code (determined at compile time).
    private
  contains
    private
    ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
    procedure(apply_precon), public, deferred :: apply
  end type poisson_precon_t

  abstract interface
    subroutine apply_precon(self, p, b, backend)
      ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
      import poisson_precon_t
      import field_t
      import base_backend_t
      implicit none
      class(poisson_precon_t) :: self
      class(field_t), intent(in) :: p    ! Pressure solution
      class(field_t), intent(inout) :: b ! The evaluated matrix-vector product
      class(base_backend_t), intent(in) :: backend
    end subroutine apply_precon
  end interface

  type, abstract, public :: poisson_solver_t
    !! Base definition of the backend-specific implementation of the iterative
    !! Poisson solver.
  contains
    private
    ! Solves the Poisson problem
    procedure(solve), public, deferred :: solve
  end type poisson_solver_t

  abstract interface
    subroutine solve(self, p, f, backend)
      import poisson_solver_t
      import field_t
      import base_backend_t
      implicit none
      class(poisson_solver_t) :: self
      class(field_t), intent(inout) :: p
      class(field_t), intent(in) :: f
      class(base_backend_t), intent(in) :: backend
    end subroutine solve
  end interface

contains

  function init_lapl(backend, mesh) result(lapl)
    !! Public constructor for the laplace_operator_t type.
    type(laplace_operator_t) :: lapl
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh

    integer :: nx, ny, nz
    real(dp) :: dx, dy, dz
    integer :: bcx1, bcxn
    integer :: bcy1, bcyn
    integer :: bcz1, bczn

    lapl%xdirps%dir = DIR_X; lapl%ydirps%dir = DIR_Y; lapl%zdirps%dir = DIR_Z

    nx = mesh%get_n(DIR_X, CELL)
    ny = mesh%get_n(DIR_Y, CELL)
    nz = mesh%get_n(DIR_Z, CELL)
    dx = mesh%geo%d(DIR_X)
    dy = mesh%geo%d(DIR_Y)
    dz = mesh%geo%d(DIR_Z)

    ! TODO: Add more sensible BCs
    bcx1 = BC_PERIODIC; bcxn = BC_PERIODIC
    bcy1 = BC_PERIODIC; bcyn = BC_PERIODIC
    bcz1 = BC_PERIODIC; bczn = BC_PERIODIC

    call backend%alloc_tdsops(lapl%xdirps%der2nd, nx, dx, &
                              "second-deriv", "compact6", bcx1, bcxn)

    call backend%alloc_tdsops(lapl%ydirps%der2nd, ny, dy, &
                              "second-deriv", "compact6", bcy1, bcyn)

    call backend%alloc_tdsops(lapl%zdirps%der2nd, nz, dz, &
                              "second-deriv", "compact6", bcz1, bczn)
  end function init_lapl

  subroutine poissmult(self, f, p, backend)
    !! Computes the action of the Laplace operator, i.e. `f = Ax` where `A` is
    !! the discrete Laplacian.
    class(laplace_operator_t) :: self
    class(field_t), intent(inout) :: f ! The output field
    class(field_t), intent(in) :: p    ! The input field
    class(base_backend_t), intent(in) :: backend

    class(field_t), pointer :: f_x, p_x
    integer :: reorder_op, reorder_op2x

    if (p%dir /= f%dir) then
      error stop "Currently orientations of P and F must match"
    end if

    if (f%dir == DIR_X) then
      call self%poissmult_dirx(f, p, backend)
    else
      f_x => backend%allocator%get_block(DIR_X, CELL)
      p_x => backend%allocator%get_block(DIR_X, CELL)

      if (f%dir == DIR_Y) then
        reorder_op2x = RDR_Y2X
        reorder_op = RDR_X2Y
      else if (f%dir == DIR_Z) then
        reorder_op2x = RDR_Z2X
        reorder_op = RDR_X2Z
      else if (f%dir == DIR_C) then
        reorder_op2x = RDR_C2X
        reorder_op = RDR_X2C
      else
        error stop "Unsupported Poisson orientation"
      end if
      call backend%reorder(p_x, p, reorder_op2x)

      call self%poissmult_dirx(f_x, p_x, backend)

      call backend%reorder(f, f_x, reorder_op)

      call backend%allocator%release_block(f_x)
      call backend%allocator%release_block(p_x)
    end if

  end subroutine poissmult

  subroutine poissmult_dirx(self, f, p, backend)
    !! Computes the action of the Laplace operator, i.e. `f = Ax` where `A` is
    !! the discrete Laplacian.
    !
    !! XXX: This requires fields in the DIR_X orientation due to use of
    !! sumintox.
    class(laplace_operator_t) :: self
    class(field_t), intent(inout) :: f ! The output field
    class(field_t), intent(in) :: p    ! The input field
    class(base_backend_t), intent(in) :: backend

    ! Compute d2pdx2
    call compute_der2nd(f, p, backend, self%xdirps%der2nd)

    ! Compute d2pdy2, d2pdz2 and accumulate
    call compute_and_acc_der2nd(f, p, backend, self%ydirps%der2nd, RDR_X2Y)
    call compute_and_acc_der2nd(f, p, backend, self%zdirps%der2nd, RDR_X2Z)

  end subroutine poissmult_dirx

  subroutine compute_and_acc_der2nd(f, p, backend, tdsops, reorder_op)
    !! Accumulates 2nd derivatives into the Laplacian

    class(field_t), intent(inout) :: f ! The Laplacian
    class(field_t), intent(in) :: p    ! The pressure field
    class(base_backend_t), intent(in) :: backend
    class(tdsops_t), intent(in) :: tdsops        ! The tridiagonal operator
    integer, intent(in) :: reorder_op  ! The reordering operation

    class(field_t), pointer :: p_i ! P in operation order
    class(field_t), pointer :: f_i ! F in operation order

    integer :: DIR

    if (reorder_op == RDR_X2Y) then
      DIR = DIR_Y
    else if (reorder_op == RDR_X2Z) then
      DIR = DIR_Z
    else
      error stop "Unsupported reordering operation"
    end if

    p_i => backend%allocator%get_block(DIR)
    f_i => backend%allocator%get_block(DIR)
    call backend%reorder(p_i, p, reorder_op)

    call compute_der2nd(f_i, p_i, backend, tdsops)
    if (reorder_op == RDR_X2Y) then
      call backend%sum_yintox(f, f_i)
    else if (reorder_op == RDR_X2Z) then
      call backend%sum_zintox(f, f_i)
    else
      error stop "Unsupported reordering operation"
    end if

    call backend%allocator%release_block(p_i)
    call backend%allocator%release_block(f_i)

  end subroutine compute_and_acc_der2nd

  subroutine compute_der2nd(d2fdx2, f, backend, tdsops)
    !! Computes the 2nd derivative of a field
    class(field_t), intent(inout) :: d2fdx2      ! The 2nd derivative
    class(field_t), intent(in) :: f              ! The field for derivative
    class(base_backend_t), intent(in) :: backend ! The backend implementation of operations
    class(tdsops_t), intent(in) :: tdsops        ! The tridiagonal operator

    call backend%tds_solve(d2fdx2, f, tdsops)

  end subroutine compute_der2nd

end module m_base_poisson_cg
