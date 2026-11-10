program test_memcheck_mirror_010
  !! Check mirror_buffer_bytes_010 against the dim2 mirror buffers that
  !! init_dim2_mirror allocates for the 010 case.

  use m_common, only: i8, nbytes
  use m_memory_estimate, only: mirror_buffer_bytes_010

  implicit none

  logical :: test_pass

  test_pass = .true.

  call test_single_rank_is_zero()
  call test_two_ranks_unstretched()
  call test_two_ranks_stretched_is_zero()

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine test_single_rank_is_zero()
    call check("1 rank, not stretched", &
               mirror_buffer_bytes_010([33, 32, 64], 1, .false.), 0_i8)
    call check("1 rank, stretched", &
               mirror_buffer_bytes_010([33, 32, 64], 1, .true.), 0_i8)
  end subroutine test_single_rank_is_zero

  subroutine test_two_ranks_unstretched()
    !! c_mirror_dev and c_slab_send_dev hold (nx, ny, nz) complex values each,
    !! c_plane_send_dev and c_plane_recv_dev hold (nx, nz) complex values each.
    integer(i8) :: expect

    expect = (2_i8*33_i8*16_i8*64_i8 + 2_i8*33_i8*64_i8)*2_i8*int(nbytes, i8)
    call check("2 ranks, not stretched", &
               mirror_buffer_bytes_010([33, 16, 64], 2, .false.), expect)
  end subroutine test_two_ranks_unstretched

  subroutine test_two_ranks_stretched_is_zero()
    call check("2 ranks, stretched", &
               mirror_buffer_bytes_010([33, 16, 64], 2, .true.), 0_i8)
  end subroutine test_two_ranks_stretched_is_zero

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

end program test_memcheck_mirror_010
