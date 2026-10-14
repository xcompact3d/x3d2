program test_mem_estimate_fit
  !! Check the padding maths behind the GPU memory estimator
  !! (verification/test_memory_estimate.f90) against hand-worked values,
  !! independently of any solver build or hardware query.

  use m_common, only: dp, i8
  use m_test_utils, only: padded_dim, padded_cells

  implicit none

  logical :: test_pass

  test_pass = .true.

  call test_padded_dim()
  call test_padded_cells()

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine expect_int(test, expect, got)
    character(len=*), intent(in) :: test
    integer, intent(in) :: expect, got

    print *, test
    print *, "- Expect: ", expect
    print *, "- Got: ", got
    if (got /= expect) test_pass = .false.
  end subroutine expect_int

  subroutine test_padded_dim()
    ! n already a multiple of sz: left unchanged.
    call expect_int("padded_dim(16, 16)", 16, padded_dim(16, 16))
    ! n one above a multiple of sz: rounds up to the next multiple.
    call expect_int("padded_dim(17, 16)", 32, padded_dim(17, 16))
    ! n not aligned to sz: rounds up to the next multiple.
    call expect_int("padded_dim(100, 16)", 112, padded_dim(100, 16))
  end subroutine test_padded_dim

  subroutine test_padded_cells()
    integer(i8) :: got

    ! Distinct x, y, z, sz so a swapped argument would be caught: x pads
    ! 100->112, y pads 50->64, z is untouched at 30.
    got = padded_cells([100, 50, 30], 16)
    print *, "padded_cells([100, 50, 30], 16)"
    print *, "- Expect: ", 215040_i8
    print *, "- Got: ", got
    if (got /= 215040_i8) test_pass = .false.
  end subroutine test_padded_cells

end program test_mem_estimate_fit
