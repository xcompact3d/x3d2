module m_tridiagsolv
  implicit none

  type :: tridiagsolv
     private
     real, allocatable :: low(:), up(:), diag(:)
   contains
     procedure, public :: solve
  end type tridiagsolv

  type, extends(tridiagsolv) :: periodic_tridiagsolv
     real, parameter :: gamma = 1.
     real, allocatable :: q(:)
   contains
     procedure, public :: solve => solve_periodic
  end type periodic_tridiagsolv

  ! Storage array for intermediate result
  real, allocatable :: m_d

contains

  subroutine ensure_work_arrays(n)
    integer, intent(in) :: n
    if(.not. allocated(m_a)) then
       allocate(m_a(n), m_b(n), m_c(n))
       return
    end if
    if(n > size(m_d)) m_d = reshape(m_d, [n])

  end subroutine ensure_work_arrays

  function construct(low, up, n) result(solver)
    type(tridiagsolv) :: solver
    real :: low(:), up(:)
    integer, intent(in) :: n
    integer :: i
    solver%n = n
    solver%up = up
    solver%low = low
    solver%diag = [ &
         & (1. - low(i) * up(i - 1), i = 2, n) &
         & ]
    call ensure_work_arrays(n)
  end function construct

  function construct_periodic(low, up, n) result(solver)
    real, intent(in) :: low, up
    integer, intent(in) :: n
    real :: u(n)
    integer :: i
    solver%n = n
    solver%up = [(up, i = 1, n)]
    solver%low = [(low, i = 1, n)]
    solver%diag = [ &
         & 1. - solver%gamma, &
         & (1., i = 2, n - 1), &
         & 1. - (up * low) / solver%gamma &
         & ]
    u = [solver%gamma, (0., i=2, n-1), up]
    solver%q = solver%tridiagsolv%solve(u)
  end function construct_periodic

  pure function solve(u)
    real, intent(in) :: u(:)
    real :: w
    integer :: i, n

    if (self%n /= size(u)) error stop
    do i = 2, self%n
       m_d(i) = d(i) - self%low(i) * d(i-1)
    end do
    solve(n) = m_d(n) / self%diag(n)
    do i = self%n-1, 1, -1
       solve(i) = (m_d(i) - self%up(i) * solve(i+1)) / self%diag(i)
    end do
  end function solve

  pure function solve_periodic(self, u)
    class(periodic_tridiagsolv), intent(in) :: self
    real, intent(in) :: u(:)
    real :: y(size(u))

    y = self%tridiagsolve%solve(u)
    associate ( low => self%low, n => self%n )
      solve = y - ((y(1) - low * y(n)) &
           & / (1. + q(1) - low * q(n))) * q
    end associate
  end function solve_periodic
end module m_tridiagsolv
