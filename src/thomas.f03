module m_tridiagsolv
  implicit none

  type :: tridiagsolv
     private
     real, allocatable :: fwd(:), bwd(:), updiag(:)
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
  interface tridiagsolv
     module procedure construct
  end interface tridiagsolv

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
