program test_memcheck_stretched_y
  !! Check stretched_y_matrix_bytes against the 010-case, non-uniform
  !! y-stretching Poisson coefficient matrices device allocation.

  use m_common, only: i8, nbytes
  use m_memory_estimate, only: stretched_y_matrix_bytes

  implicit none

  logical :: test_pass

  test_pass = .true.

  call test_uniform_is_zero()
  call test_non_010_is_zero()
  call test_top_bottom_not_lowmem()
  call test_top_bottom_lowmem()
  call test_centred_equals_top_bottom()
  call test_odd_ny_spec()
  call test_ng2_halves_ny_spec()

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine test_uniform_is_zero()
    integer(i8) :: bytes

    bytes = stretched_y_matrix_bytes(.true., 'uniform', .false., &
                                      [512, 256, 256], 1)
    call check("uniform, 010, is zero", bytes, 0_i8)
  end subroutine test_uniform_is_zero

  subroutine test_non_010_is_zero()
    integer(i8) :: bytes

    bytes = stretched_y_matrix_bytes(.false., 'top-bottom', .false., &
                                      [512, 256, 256], 1)
    call check("top-bottom, not 010, is zero", bytes, 0_i8)
  end subroutine test_non_010_is_zero

  subroutine test_top_bottom_not_lowmem()
    integer(i8) :: bytes, expect

    bytes = stretched_y_matrix_bytes(.true., 'top-bottom', .false., &
                                      [512, 256, 256], 1)
    expect = 20_i8*257_i8*256_i8*256_i8*int(nbytes, i8)
    call check("top-bottom, lowmem_fft=F", bytes, expect)
  end subroutine test_top_bottom_not_lowmem

  subroutine test_top_bottom_lowmem()
    integer(i8) :: bytes, expect

    bytes = stretched_y_matrix_bytes(.true., 'top-bottom', .true., &
                                      [512, 256, 256], 1)
    expect = 10_i8*257_i8*256_i8*256_i8*int(nbytes, i8)
    call check("top-bottom, lowmem_fft=T", bytes, expect)
  end subroutine test_top_bottom_lowmem

  subroutine test_centred_equals_top_bottom()
    integer(i8) :: bytes_centred, bytes_top_bottom

    bytes_centred = stretched_y_matrix_bytes(.true., 'centred', .true., &
                                              [512, 256, 256], 1)
    bytes_top_bottom = stretched_y_matrix_bytes(.true., 'top-bottom', .true., &
                                                 [512, 256, 256], 1)
    call check("centred equals top-bottom", bytes_centred, bytes_top_bottom)
  end subroutine test_centred_equals_top_bottom

  subroutine test_odd_ny_spec()
    !! cdims (16,9,8): nx_spec=9, ny_spec=9, nz_spec=8 - pins the symmetric
    !! integer division ny_spec/2=4.
    integer(i8) :: bytes, expect

    bytes = stretched_y_matrix_bytes(.true., 'bottom', .true., [16, 9, 8], 1)
    expect = 2_i8*9_i8*9_i8*8_i8*5_i8*int(nbytes, i8)
    call check("odd ny_spec, bottom", bytes, expect)

    bytes = stretched_y_matrix_bytes(.true., 'top-bottom', .true., &
                                      [16, 9, 8], 1)
    expect = 4_i8*9_i8*4_i8*8_i8*5_i8*int(nbytes, i8)
    call check("odd ny_spec, top-bottom", bytes, expect)
  end subroutine test_odd_ny_spec

  subroutine test_ng2_halves_ny_spec()
    integer(i8) :: bytes, expect

    bytes = stretched_y_matrix_bytes(.true., 'top-bottom', .true., &
                                      [512, 256, 256], 2)
    expect = 4_i8*257_i8*64_i8*256_i8*5_i8*int(nbytes, i8)
    call check("ng=2 halves ny_spec", bytes, expect)
  end subroutine test_ng2_halves_ny_spec

  subroutine check(label, actual, expect)
    character(*), intent(in) :: label
    integer(i8), intent(in) :: actual, expect

    print *, label
    print *, "- Expect: ", expect
    print *, "- Got: ", actual
    if (actual /= expect) then
      test_pass = .false.
    end if
  end subroutine check

end program test_memcheck_stretched_y
