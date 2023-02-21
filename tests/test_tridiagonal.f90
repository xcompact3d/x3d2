program test_tridiagonal
  use iso_fortran_env, only: stderr => error_unit

  use m_tridiagsolv, only: tridiagsolv, periodic_tridiagsolv

  implicit none

  class(tridiagsolv), allocatable :: solver

  integer :: i
  real :: low(5)
  real :: up(5)
  real :: f(2, 6), df(2, 6), expected_vec(2, 6)
  real :: expected(6)

  real, parameter :: tol = 0.0001
  logical :: allpass

  allpass = .true.

  ! DIRICHLET BOUNDARY CONDITIONS - Solve for x in

  ! 1    2.   0     0     0     0         x1    0
  ! 0.25 1    0.25  0     0     0         x2    1
  ! 0    0.33 1     0.33  0     0     X   x3  = 2
  ! 0    0    0.33  1     0.33  0         x4    3
  ! 0    0    0     0.25  1     0.25      x5    2
  ! 0    0    0     0     2     1         x6    1


  ! Solution computed with numpy.linalg.solve()
  expected = [ &
       & -1.71428, &
       & 0.85714, &
       & 2.28571, &
       & 1.28571, &
       & 2.85714, &
       & -4.71428 &
       & ]
  expected_vec(1, :) = expected
  expected_vec(2, :) = expected

  low = [1. / 4., 1. / 3., 1. / 3., 1. / 4., 2.]
  up = [2., 1. / 4., 1. / 3., 1. / 3., 1. / 4.]
  solver = tridiagsolv(low, up)

  f(1, :) = [0., 1., 3., 3., 2., 1.]
  f(2, :) = [0., 1., 3., 3., 2., 1.]
  call solver%solve(f, df)

  if (.not. all(abs(df - expected_vec) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Tridiagonal system is solved correctly (DIRICHLET)... failed'
  else
     write(stderr, '(a)') 'Tridiagonal system is solved correctly (DIRICHLET)... passed'
  end if
  deallocate(solver)

  ! PERIODIC BOUNDARY CONDITIONS - Solve for x in

  ! 1    0.33 0     0     0     0.33      x1    0
  ! 0.33 1    0.33  0     0     0         x2    1
  ! 0    0.33 1     0.33  0     0     X   x3  = 2
  ! 0    0    0.33  1     0.33  0         x4    3
  ! 0    0    0     0.33  1     0.33      x5    2
  ! 0.33 0    0     0     0.33  1         x6    1

  expected = [-0.375, 0.375, 2.25, 1.875, 1.125, 0.75]
  expected_vec(1, :) = expected
  expected_vec(2, :) = expected

  low = [(1. / 3., i = 1, 5)]
  up = [(1. / 3., i = 1, 5)]
  allocate(periodic_tridiagsolv::solver)
  solver = periodic_tridiagsolv(low, up)

  f(1, :) = [0., 1., 3., 3., 2., 1.]
  f(2, :) = [0., 1., 3., 3., 2., 1.]
  call solver%solve(f, df)

  if (.not. all(abs(df - expected_vec) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Tridiagonal system is solved correctly (PERIODIC)... failed'
  else
     write(stderr, '(a)') 'Tridiagonal system is solved correctly (PERIODIC)... passed'
  end if
  
end program test_tridiagonal
