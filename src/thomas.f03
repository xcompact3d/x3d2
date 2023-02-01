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
     real, allocatable :: q(:, :)
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

  function construct_periodic(low, up, n, SZ) result(solver)
    real, intent(in) :: low, up
    integer, intent(in) :: n, SZ
    real :: u(n), u_vec(SZ, n)
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
    ! SZ = 4, n =8
    ! u_vec(1, :) = gamma 0 0 0 0 0 0 up
    ! u_vec(2, :) = gamma 0 0 0 0 0 0 up
    ! u_vec(3, :) = gamma 0 0 0 0 0 0 up
    ! u_vec(4, :) = gamma 0 0 0 0 0 0 up
    u_vec = transpose(reshape(u, [n, SZ], pad = u))
    solver%q = solver%tridiagsolv%solve(u_vec)
  end function construct_periodic

  pure subroutine solve(self, f, df)
    class(tridiagsolv), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    integer :: i, n

    if (self%n /= size(f, 2)) error stop
    do j = 2, self%n
       !$omp simd
       do i = 1, size(u, 1)
          m_f(j) = f(i, j) - self%low(j) * f(i, j-1)
       end do
       !$omp end simd
    end do
    !$omp simd
    do i = 1, size(f, 1)
       df(i, n) = m_f(n) / self%diag(n)
    end do
    !$omp end simd

    do j = self%n-1, 1, -1
       !$omp simd
       do i = 1, size(f, 1)
          df(i, j) = (m_f(j) - self%up(j) * df(i, j+1)) / self%diag(i)
       end do
       !$omp end simd
    end do
  end subroutine solve

  pure subroutine solve_periodic(self, f, df)
    class(periodic_tridiagsolv), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    real :: y(size(f, 1), size(f, 2))

    y = self%tridiagsolve%solve(u)
    associate ( low => self%low, n => self%n )
      !$omp simd
      do i = 1, size(f, 1)
         df(i, :) = y(i, :) - ((y(i, 1) - low * y(i, n)) &
              & / (1. + q(i, 1) - low * q(i, n))) * q
      end do
      !$omp end simd
    end associate
  end subroutine solve_periodic
end module m_tridiagsolv
