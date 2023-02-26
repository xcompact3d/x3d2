program test_diffengine
  use iso_fortran_env, only: stderr => error_unit

  use m_diffengine, only: diffengine_t

  implicit none

  type(diffengine_t) :: diffeng
  integer, parameter :: n = 20
  integer :: i
  real, allocatable :: f(:), expected(:)
  real, allocatable :: f_vec(:, :), expected_vec(:, :), df(:, :)
  real :: dx
  logical :: allpass
  integer, parameter :: SZ = 4
  real, parameter :: tol = 0.01

  dx = 2. * acos(-1.) / (n - 1)
  f = [(sin((i-1)*dx), i=1,n)]
  f_vec = transpose(reshape(f, [n, SZ], pad = f))
  allpass = .true.

  ! First derivative with periodic boundary conditions
  expected = [(cos((i-1)*dx), i=1,n)]
  expected_vec = transpose(reshape(expected, [n, SZ], pad = expected))
  diffeng = diffengine_t("compact6", length=n, order=1)
  allocate(df, source=f_vec)
  call diffeng%diff(f_vec, df, dx)
  if (.not. all(abs(df - expected_vec) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'First derivatives (periodic) are computed correctly... failed'
  else
     write(stderr, '(a)') 'First derivatives (periodic) are computed correctly... passed'
  end if

  if (allpass) then
     write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
     error stop 'SOME TESTS FAILED.'
  end if
end program test_diffengine
