program test_snapshot_append
  !! Unit test for restarting a snapshot file: count_steps_up_to finds the
  !! steps to keep, and an append session drops the later ones before
  !! writing new steps after them. read_string_attribute returns the stored
  !! vtk.xml unchanged, so a restart can tell whether its layout matches.
  use m_common, only: dp, i8
  use m_io_backend, only: count_steps_up_to, read_string_attribute
  use m_io_session, only: writer_session_t
  use m_test_utils, only: initialise_mpi, finalise_test, global_all
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier
  implicit none

  character(len=*), parameter :: test_file = "test_snapshot_append.bp"
  character(len=*), parameter :: vtk_xml = '<?xml version="1.0"?>'// &
                                 new_line('a')//'<VTKFile/>'

  integer :: ierr, irank, nproc
  logical :: allpass = .true.

  call initialise_mpi(irank, nproc)

  ! The run being restarted wrote five steps
  call write_steps([0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp])

  call check_count(-1.0_dp, 0, "limit before the first step")
  call check_count(0.0_dp, 1, "limit equal to the first time")
  call check_count(0.25_dp, 3, "limit between two steps")
  call check_count(1.0_dp, 5, "limit after the last step")
  call check_missing_variable()
  call check_attribute()

  ! The restart keeps the first three steps and writes two more
  call write_steps([0.25_dp, 0.3_dp], append_after_steps=3)

  ! Expect 0, 0.1, 0.2, 0.25, 0.3. Had the old 0.3 and 0.4 been kept, the
  ! count would stop at 0.3 and the total would be seven.
  call check_count(0.27_dp, 4, "after append, limit inside new steps")
  call check_count(1.0_dp, 5, "after append, total steps")

  ! A restart that keeps nothing replaces the file
  call write_steps([0.5_dp], append_after_steps=0)
  call check_count(1.0_dp, 1, "after append with no kept steps")

  call global_all(allpass)

  if (irank == 0) call execute_command_line("rm -rf "//test_file)

  call finalise_test(allpass, irank)

contains

  subroutine write_steps(times, append_after_steps)
    !! Write one step per time, as the snapshot manager does
    real(dp), intent(in) :: times(:)
    integer, intent(in), optional :: append_after_steps

    type(writer_session_t) :: writer
    real(dp) :: field(4, 3, 1)
    integer(i8) :: global_dims(3), local_start(3), local_count(3)
    integer :: n

    global_dims = [4_i8, 3_i8, int(nproc, i8)]
    local_start = [0_i8, 0_i8, int(irank, i8)]
    local_count = [4_i8, 3_i8, 1_i8]

    if (present(append_after_steps)) then
      call writer%open(test_file, MPI_COMM_WORLD, append_after_steps)
    else
      call writer%open(test_file, MPI_COMM_WORLD)
    end if
    call writer%write_attribute("vtk.xml", vtk_xml)

    do n = 1, size(times)
      if (n > 1) call writer%begin_step()
      field = times(n)
      call writer%write_data("time", times(n))
      call writer%write_data("u", field, global_dims, local_start, &
                             local_count)
      call writer%end_step()
    end do

    call writer%close()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine write_steps

  subroutine check_count(limit, expected, label)
    real(dp), intent(in) :: limit
    integer, intent(in) :: expected
    character(len=*), intent(in) :: label

    integer :: n

    n = count_steps_up_to(test_file, "time", limit, MPI_COMM_WORLD)
    if (n /= expected) then
      allpass = .false.
      if (irank == 0) print *, label, ": got ", n, ", expected ", expected
    end if
  end subroutine check_count

  subroutine check_missing_variable()
    integer :: n

    n = count_steps_up_to(test_file, "not_a_variable", 1.0_dp, &
                          MPI_COMM_WORLD)
    if (n /= 0) then
      allpass = .false.
      if (irank == 0) print *, "missing variable: got ", n, ", expected 0"
    end if
  end subroutine check_missing_variable

  subroutine check_attribute()
    character(len=:), allocatable :: value

    value = read_string_attribute(test_file, "vtk.xml", MPI_COMM_WORLD)
    if (value /= vtk_xml) then
      allpass = .false.
      if (irank == 0) print *, "vtk.xml read back as: [", value, "]"
    end if

    value = read_string_attribute(test_file, "not_an_attribute", &
                                  MPI_COMM_WORLD)
    if (value /= '') then
      allpass = .false.
      if (irank == 0) print *, "missing attribute read as: [", value, "]"
    end if
  end subroutine check_attribute

end program test_snapshot_append
