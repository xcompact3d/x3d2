program test_reorder_map
  !! Check that reorder mappings return the expected values

  use m_common
  
  implicit none

  logical :: test_pass

  test_pass = .true.
  
  call test_get_dirs_from_rdr()
  call test_get_rdr_from_dirs()

  if (.not. test_pass) then
    error stop "FAIL"
  end if
  
contains

  subroutine test_get_dirs_from_rdr()
    
    call test_one_get_dirs_from_rdr("RDR_X2Y", RDR_X2Y, DIR_X, DIR_Y)
    call test_one_get_dirs_from_rdr("RDR_X2Z", RDR_X2Z, DIR_X, DIR_Z)
    call test_one_get_dirs_from_rdr("RDR_X2C", RDR_X2C, DIR_X, DIR_C)
    
    call test_one_get_dirs_from_rdr("RDR_Y2X", RDR_Y2X, DIR_Y, DIR_X)
    call test_one_get_dirs_from_rdr("RDR_Y2Z", RDR_Y2Z, DIR_Y, DIR_Z)
    call test_one_get_dirs_from_rdr("RDR_Y2C", RDR_Y2C, DIR_Y, DIR_C)
    
    call test_one_get_dirs_from_rdr("RDR_Z2X", RDR_Z2X, DIR_Z, DIR_X)
    call test_one_get_dirs_from_rdr("RDR_Z2Y", RDR_Z2Y, DIR_Z, DIR_Y)
    call test_one_get_dirs_from_rdr("RDR_Z2C", RDR_Z2C, DIR_Z, DIR_C)
    
    call test_one_get_dirs_from_rdr("RDR_C2X", RDR_C2X, DIR_C, DIR_X)
    call test_one_get_dirs_from_rdr("RDR_C2Y", RDR_C2Y, DIR_C, DIR_Y)
    call test_one_get_dirs_from_rdr("RDR_C2Z", RDR_C2Z, DIR_C, DIR_Z)
    
  end subroutine test_get_dirs_from_rdr

  subroutine test_one_get_dirs_from_rdr(test, rdr, expect_from, expect_to)
    character(len=*), intent(in) :: test
    integer, intent(in) :: rdr, expect_from, expect_to

    integer :: dir_from, dir_to

    call get_dirs_from_rdr(dir_from, dir_to, rdr)
    print *, test
    print *, "- Expect: ", expect_from, expect_to
    print *, "- Got: ", dir_from, dir_to
    if ((dir_from /= expect_from) .or. (dir_to /= expect_to)) then
      test_pass = .false.
    end if

  end subroutine test_one_get_dirs_from_rdr

  subroutine test_get_rdr_from_dirs()

    call test_one_get_rdr_from_dirs("X->Y", DIR_X, DIR_Y, RDR_X2Y)
    call test_one_get_rdr_from_dirs("X->Z", DIR_X, DIR_Z, RDR_X2Z)
    call test_one_get_rdr_from_dirs("X->C", DIR_X, DIR_C, RDR_X2C)

    call test_one_get_rdr_from_dirs("Y->X", DIR_Y, DIR_X, RDR_Y2X)
    call test_one_get_rdr_from_dirs("Y->Z", DIR_Y, DIR_Z, RDR_Y2Z)
    call test_one_get_rdr_from_dirs("Y->C", DIR_Y, DIR_C, RDR_Y2C)

    call test_one_get_rdr_from_dirs("Z->X", DIR_Z, DIR_X, RDR_Z2X)
    call test_one_get_rdr_from_dirs("Z->Y", DIR_Z, DIR_Y, RDR_Z2Y)
    call test_one_get_rdr_from_dirs("Z->C", DIR_Z, DIR_C, RDR_Z2C)

    call test_one_get_rdr_from_dirs("C->X", DIR_C, DIR_X, RDR_C2X)
    call test_one_get_rdr_from_dirs("C->Y", DIR_C, DIR_Y, RDR_C2Y)
    call test_one_get_rdr_from_dirs("C->Z", DIR_C, DIR_Z, RDR_C2Z)
    
  end subroutine test_get_rdr_from_dirs

  subroutine test_one_get_rdr_from_dirs(test, from, to, expect_rdr)

    character(len=*), intent(in) :: test
    integer, intent(in) :: from, to, expect_rdr
    integer :: rdr

    rdr = get_rdr_from_dirs(from, to)
    print *, test
    print *, "- Expect: ", expect_rdr
    print *, "- Got: ", rdr
    if (rdr /= expect_rdr) then
      test_pass = .false.
    end if

  end subroutine test_one_get_rdr_from_dirs
  
end program test_reorder_map
