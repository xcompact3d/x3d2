program test_poisson_bc_code
  !! Checks poisson_bc_code on every periodicity combination, and
  !! spectral_dirs on the layout of each case the FFT Poisson solver runs.

  use m_poisson_fft, only: poisson_bc_code, spectral_dirs

  implicit none

  logical :: test_pass

  test_pass = .true.

  call test_bc_code([.true., .true., .true.], '000')
  call test_bc_code([.true., .false., .true.], '010')
  call test_bc_code([.false., .true., .true.], '100')
  call test_bc_code([.false., .false., .true.], '110')
  call test_bc_code([.true., .true., .false.], '001')
  call test_bc_code([.true., .false., .false.], '011')
  call test_bc_code([.false., .true., .false.], '101')
  call test_bc_code([.false., .false., .false.], '111')

  call test_spectral_dirs('000', [1, 2, 3])
  call test_spectral_dirs('010', [1, 2, 3])
  call test_spectral_dirs('100', [2, 1, 3])
  call test_spectral_dirs('110', [3, 1, 2])

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine test_bc_code(periodic, expect)
    logical, intent(in) :: periodic(3)
    character(len=3), intent(in) :: expect

    character(len=3) :: got

    got = poisson_bc_code(periodic)
    print *, "poisson_bc_code", periodic
    print *, "- Expect: ", expect
    print *, "- Got: ", got
    if (got /= expect) test_pass = .false.

  end subroutine test_bc_code

  subroutine test_spectral_dirs(bc_code, expect)
    character(len=3), intent(in) :: bc_code
    integer, intent(in) :: expect(3)

    integer :: got(3)

    got = spectral_dirs(bc_code)
    print *, "spectral_dirs ", bc_code
    print *, "- Expect: ", expect
    print *, "- Got: ", got
    if (any(got /= expect)) test_pass = .false.

  end subroutine test_spectral_dirs

end program test_poisson_bc_code
