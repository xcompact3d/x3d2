module m_tridiagsolv
  implicit none

  type :: tridiagsolv
     private
     real, allocatable :: fwd(:), bwd(:), updiag(:)
   contains
     procedure, public :: solve
  end type tridiagsolv

  type, extends(tridiagsolv) :: periodic_tridiagsolv
     real :: gamma = - 1.
     real, allocatable :: q(:)
   contains
     procedure, public :: solve => solve_periodic
  end type periodic_tridiagsolv

  interface tridiagsolv
     module procedure construct
  end interface tridiagsolv

  interface periodic_tridiagsolv
     module procedure construct_periodic
  end interface periodic_tridiagsolv

contains

  function construct(low, up) result(solver)
    type(tridiagsolv) :: solver
    real, intent(in) :: low(:), up(:)
    real, allocatable :: fwd(:), i_bwd(:)
    integer :: i, n

    n = size(low) + 1
    ! Allocate solver's coefficient arrays.
    allocate(solver%fwd(n), solver%bwd(n), solver%updiag(n))
    solver%updiag = up
    ! Allocate extra temp arrays for readability.  These only live
    ! within this the function scop.
    allocate(fwd(n), i_bwd(n))

    i_bwd = 1. ! Initialise bwd array to diagonal.

    !i_bwd is the /inverse/ of the coefficients for the bwd step.
    do i = 2, n
       fwd(i) = low(i - 1) / i_bwd(i - 1)
       i_bwd(i) = i_bwd(i) - up(i - 1) * fwd(i)
    end do
    solver%fwd = fwd
    solver%bwd = 1. / i_bwd
  end function construct

  function construct_periodic(low, up) result(solver)
    real, intent(in) :: low(:), up(:)
    type(periodic_tridiagsolv) :: solver

    real, allocatable :: u(:), u_vec(:, :), q_vec(:,:)
    real, allocatable :: uu(:)
    integer :: i, n
    n = size(low) + 1
    ! Runs non-periodic (base type) constructor which sets up forward
    ! and backward coefficient arrays.
    solver%tridiagsolv = tridiagsolv(low, up)
    ! q member array is used in the Sermann-Morrison formula

    ! solve() method expects a rank 2 array of shape (SZ X n) as
    ! input, where SZ is the number of pencils in the pencil group.
    ! Therefore artificially work a temp array u_vec and q_vec of rank
    ! 2.
    allocate(u(n), u_vec(1, n), q_vec(1, n))
    u = [solver%gamma, (0., i=2, n-1), up(1)]
    u_vec = transpose(reshape(u, [n, 1], pad = u))
    uu = u_vec(1, :)
    call solver%tridiagsolv%solve(u_vec, q_vec)
    solver%q = q_vec(1, :)
  end function construct_periodic

  pure subroutine solve(self, f, df)
    class(tridiagsolv), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    integer :: i, j, n

    df = f
    n = size(self%fwd)
    if (size(f, 2) /= n) error stop
    do j = 2, n
       !$omp simd
       do i = 1, size(f, 1)
          df(i, j) = df(i, j) - df(i, j-1) * self%fwd(j)
       end do
       !$omp end simd
    end do
    !$omp simd
    do i = 1, size(f, 1)
       df(i, n) = df(i, n) * self%bwd(n)
    end do
    !$omp end simd

    do j = n-1, 1, -1
       !$omp simd
       do i = 1, size(f, 1)
          df(i, j) = (df(i, j) - df(i, j+1) * self%updiag(j)) &
               & * self%bwd(j)
       end do
       !$omp end simd
    end do
  end subroutine solve

  pure subroutine solve_periodic(self, f, df)
    class(periodic_tridiagsolv), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    real :: y(size(f, 1), size(f, 2))
    integer :: i, n
    real :: alpha

    call self%tridiagsolv%solve(f, y)

    n = size(f, 2)
    alpha = self%updiag(1)
    select type (self)
    type is (periodic_tridiagsolv)
       !$omp simd
       do i = 1, size(f, 1)
          df(i, :) = y(i, :) - ((y(i, 1) - alpha * y(i, n)) &
               & / (1. + self%q(1) - alpha * self%q(n))) * self%q
       end do
       !$omp end simd
    class default
       error stop
    end select
  end subroutine solve_periodic
end module m_tridiagsolv
