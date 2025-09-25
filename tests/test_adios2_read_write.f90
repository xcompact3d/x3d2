program test_adios2
  use mpi
  use m_adios2_io, only: adios2_writer_t, adios2_reader_t, &
                         adios2_mode_write, adios2_mode_read, &
                         adios2_file_t
  use m_common, only: dp
  use iso_fortran_env, only: int64, real64, stderr => error_unit
  implicit none

  ! ADIOS2 handlers
  type(adios2_writer_t) :: adios2_writer
  type(adios2_reader_t) :: adios2_reader
  type(adios2_file_t) :: file

  ! MPI variables
  integer :: ierr, irank, isize
  integer(kind=int64), dimension(2) :: shape_dims, start_dims, count_dims
  integer(kind=int64), dimension(2) :: sel_start, sel_count
  real(kind=dp), dimension(:, :), allocatable :: data_write, data_read

  ! application variables
  integer :: i, j, rank_id, inx = 3, iny = 4
  logical :: allpass = .true.
  real(kind=dp) :: expected

  ! launch MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

  ! data initialization
  allocate (data_write(inx, iny))

  ! initialize data (unique per rank)
  do j = 1, iny
    do i = 1, inx
      data_write(i, j) = irank*inx*iny + (j - 1)*inx + (i - 1)
    end do
  end do

  ! global shape and local offsets for parallel I/O
  shape_dims = [isize*inx, iny]
  start_dims = [irank*inx, 0]
  count_dims = [inx, iny]

  ! write data
  call adios2_writer%init(MPI_COMM_WORLD, "test_io_write")
  file = adios2_writer%open("test_output.bp", adios2_mode_write)
  call adios2_writer%begin_step(file)
  call adios2_writer%write_data("data2D", data_write, file, &
                                shape_dims, start_dims, count_dims)
  call adios2_writer%end_step(file)
  call adios2_writer%close(file)

  if (allocated(data_write)) deallocate (data_write)

  ! read data (rank 0 only)
  if (irank == 0) then
    call adios2_reader%init(MPI_COMM_SELF, "test_io_read")
    file = adios2_reader%open("test_output.bp", adios2_mode_read)
    call adios2_reader%begin_step(file)

    sel_start = [0, 0]
    sel_count = [shape_dims(1), shape_dims(2)]

    ! read entire dataset
    allocate (data_read(sel_count(1), sel_count(2)))
    call adios2_reader%read_data("data2D", data_read, file, &
                                 sel_start, sel_count)

    ! verify data
    do j = 1, sel_count(2)
      do i = 1, sel_count(1)
        ! calculate expected value based on original data pattern
        rank_id = (i - 1)/inx ! determine which rank wrote this data
        expected = rank_id*inx*iny + (j - 1)*inx + mod(i - 1, inx)
        if (abs(data_read(i, j) - expected) > 1.e-8_dp) then
          allpass = .false.
          print *, "Data mismatch at (", i, ",", j, "): ", data_read(i, j)
        end if
      end do
    end do

    call adios2_reader%end_step(file)
    call adios2_reader%close(file)

    if (allocated(data_read)) deallocate (data_read)
  end if

  ! Cleanup and finalize
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)

  ! Test result
  if (irank == 0) then
    if (allpass) then
      write (stderr, '(a)') 'ADIOS2 TEST PASSED SUCCESSFULLY.'
    else
      error stop 'ADIOS2 TEST FAILED.'
    end if
  end if

end program test_adios2
