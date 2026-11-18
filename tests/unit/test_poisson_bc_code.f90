program test_poisson_bc_code
  !! Checks poisson_bc_code on every periodicity combination, and
  !! spectral_dirs on the layout of each case the FFT Poisson solver runs.
  !! Also checks that the FFT_BC_* class the CUDA solver derives from
  !! poisson_bc_code matches the one x3d2-memcheck gets from
  !! classify_boundaries.

  use m_mesh, only: classify_boundaries, FFT_BC_000, FFT_BC_010, &
                    FFT_BC_100, FFT_BC_110, FFT_BC_UNSUPPORTED
  use m_poisson_fft, only: poisson_bc_code, spectral_dirs, fft_bc_class

  implicit none

  logical :: test_pass

  test_pass = .true.

  call test_bc_code([.true., .true., .true.], '000', FFT_BC_000)
  call test_bc_code([.true., .false., .true.], '010', FFT_BC_010)
  call test_bc_code([.false., .true., .true.], '100', FFT_BC_100)
  call test_bc_code([.false., .false., .true.], '110', FFT_BC_110)
  call test_bc_code([.true., .true., .false.], '001', FFT_BC_UNSUPPORTED)
  call test_bc_code([.true., .false., .false.], '011', FFT_BC_UNSUPPORTED)
  call test_bc_code([.false., .true., .false.], '101', FFT_BC_UNSUPPORTED)
  call test_bc_code([.false., .false., .false.], '111', FFT_BC_UNSUPPORTED)

  call test_spectral_dirs('000', [1, 2, 3])
  call test_spectral_dirs('010', [1, 2, 3])
  call test_spectral_dirs('100', [2, 1, 3])
  call test_spectral_dirs('110', [3, 1, 2])

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine test_bc_code(periodic, expect, expect_class)
    logical, intent(in) :: periodic(3)
    character(len=3), intent(in) :: expect
    integer, intent(in) :: expect_class

    character(len=3) :: got
    character(len=8) :: bcs(3, 2)
    logical :: periodic_mesh(3)
    integer :: dir, got_class, mesh_class

    got = poisson_bc_code(periodic)
    print *, "poisson_bc_code", periodic
    print *, "- Expect: ", expect
    print *, "- Got: ", got
    if (got /= expect) test_pass = .false.

    ! The class the CUDA solver uses, and the one m_mesh gives from the BC
    ! strings of the same periodicity
    got_class = fft_bc_class(got)
    do dir = 1, 3
      if (periodic(dir)) then
        bcs(dir, :) = 'periodic'
      else
        bcs(dir, :) = 'neumann'
      end if
    end do
    call classify_boundaries(bcs(1, :), bcs(2, :), bcs(3, :), &
                             periodic_mesh, mesh_class)
    print *, "- Expect class: ", expect_class
    print *, "- Got class (fft_bc_class): ", got_class
    print *, "- Got class (classify_boundaries): ", mesh_class
    if (got_class /= expect_class) test_pass = .false.
    if (mesh_class /= expect_class) test_pass = .false.
    if (any(periodic_mesh .neqv. periodic)) test_pass = .false.

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
