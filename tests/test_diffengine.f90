program test_diffengine
  use iso_fortran_env, only: stderr => error_unit

  use m_diffengine, only: diffengine_t

  implicit none

  type(diffengine_t) :: diffeng
  integer, parameter :: n = 20
  integer :: i, j
  real, allocatable :: f_vec(:, :), expected_vec(:, :), df(:, :)
  real :: dx
  logical :: allpass = .true.
  integer, parameter :: SZ = 4
  real, parameter :: tol = 0.01

  dx = 2. * acos(-1.) / (n - 1)
  allocate(f_vec(SZ,n))
  do j = 1, n
     do i = 1, SZ
        f_vec(i, j) = sin((j-1)*dx)
     end do
  end do

  ! First derivative with periodic boundary conditions
  allocate(expected_vec(SZ,n))
  do j = 1, n
     do i = 1, SZ
        expected_vec(i, j) = cos((j-1)*dx)
     end do
  end do

  diffeng = diffengine_t("compact6", length=n, order=1, dx=dx)
  allocate(df, source=f_vec)
  call diffeng%diff(f_vec, df)
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
