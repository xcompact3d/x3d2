program test_io_session
  !! Unit test for parallel I/O session functionality using single shared file
  use mpi
  use m_common, only: dp, i8
  use m_io_session, only: writer_session_t, reader_session_t
  use iso_fortran_env, only: stderr => error_unit
  implicit none

  integer, parameter, dimension(3) :: local_dims = [8, 6, 1]
  integer(i8), dimension(3) :: global_dims, local_start, local_count
  character(len=*), parameter :: test_file = "test_parallel_session.bp"

  integer :: ierr, irank, nproc
  real(dp), dimension(local_dims(1), local_dims(2), local_dims(3)) :: &
    write_data, read_data
  integer :: timestep_write = 42, timestep_read = 0
  real(dp) :: time_write = 3.14_dp, time_read = 0.0_dp
  
  type(writer_session_t) :: writer
  type(reader_session_t) :: reader
  
  logical :: allpass = .true.
  integer :: i, j, k

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! setup global domain dimensions
  global_dims = [int(local_dims(1), i8), int(local_dims(2), i8), int(nproc, i8)]
  local_start = [0_i8, 0_i8, int(irank, i8)]
  local_count = [int(local_dims(1), i8), int(local_dims(2), i8), 1_i8]

  ! create unique data per rank
  do k = 1, local_dims(3)
    do j = 1, local_dims(2)
      do i = 1, local_dims(1)
        write_data(i, j, k) = real(i + j*10 + (irank+1)*100, dp)
      end do
    end do
  end do

  ! parallel write to single file
  call writer%open(test_file, MPI_COMM_WORLD)
  if (irank == 0) then
    call writer%write_data("timestep", timestep_write)
    call writer%write_data("time", time_write)
  end if
  call writer%write_data("velocity_u", write_data, global_dims, &
    local_start, local_count)
  call writer%close()
  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! parallel read from single file  
  call reader%open(test_file, MPI_COMM_WORLD)
  if (irank == 0) then
    call reader%read_data("timestep", timestep_read)
    call reader%read_data("time", time_read)
  end if
  call reader%read_data("velocity_u", read_data, local_start, &
    local_count, global_dims)
  call reader%close()
  
  ! validation (only rank 0 checks scalar data)
  if (irank == 0) then
    if (timestep_read /= timestep_write) then
      allpass = .false.
      print *, "Timestep mismatch:", timestep_read, "!=", timestep_write
    end if
    
    if (abs(time_read - time_write) > 1e-12_dp) then
      allpass = .false.
      print *, "Time mismatch:", time_read, "!=", time_write
    end if
  end if
  
  ! all ranks validate their local data
  do k = 1, local_dims(3)
    do j = 1, local_dims(2)
      do i = 1, local_dims(1)
        if (abs(read_data(i,j,k) - write_data(i,j,k)) > 1e-12_dp) then
          allpass = .false.
          print *, "Rank", irank, "Data mismatch at:", i, j, k
          print *, "  Expected:", write_data(i,j,k)
          print *, "  Got:     ", read_data(i,j,k)
          exit
        end if
      end do
      if (.not. allpass) exit
    end do
    if (.not. allpass) exit
  end do

  call MPI_Allreduce(MPI_IN_PLACE, allpass, 1, MPI_LOGICAL, &
    MPI_LAND, MPI_COMM_WORLD, ierr)
  
  ! cleanup
  if (irank == 0) call execute_command_line("rm -rf " // test_file)
  
  call MPI_Finalize(ierr)

  if (allpass) then
    if (irank == 0) write (stderr, &
      '(a)') 'PARALLEL I/O SESSION TEST PASSED SUCCESSFULLY.'
  else
    error stop 'PARALLEL I/O SESSION TEST FAILED.'
  end if

end program test_io_session
