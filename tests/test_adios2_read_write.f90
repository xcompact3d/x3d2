program test_adios2
  use mpi
  use m_io_backend, only: allocate_io_writer, allocate_io_reader
  use m_io_base, only: io_writer_t, io_reader_t, io_file_t, io_mode_write, &
                       io_mode_read
  use m_common, only: dp, i8, is_sp
  use iso_fortran_env, only: stderr => error_unit
  implicit none

  ! ADIOS2 handlers
  class(io_writer_t), allocatable :: adios2_writer
  class(io_reader_t), allocatable :: adios2_reader
  class(io_file_t), allocatable :: file

  ! MPI variables
  integer :: ierr, irank, isize
  integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
  integer(i8), dimension(3) :: sel_start, sel_count
  real(dp), dimension(:, :, :), allocatable :: data_write, data_read

  ! application variables
  integer :: i, j, k, rank_id, inx = 3, iny = 4, inz = 2
  logical :: allpass = .true.
  real(dp) :: expected, tolerance

  ! launch MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

  ! data initialisation
  allocate (data_write(inx, iny, inz))

  ! initialize data (unique per rank)
  do k = 1, inz
    do j = 1, iny
      do i = 1, inx
        data_write(i, j, k) = real(irank*inx*iny*inz + (k - 1)*inx*iny &
                           + (j - 1)*inx + (i - 1), dp)
      end do
    end do
  end do

  ! global shape and local offsets for parallel I/O
  shape_dims = [int(isize*inx, i8), int(iny, i8), int(inz, i8)]
  start_dims = [int(irank*inx, i8), 0_i8, 0_i8]
  count_dims = [int(inx, i8), int(iny, i8), int(inz, i8)]

  ! write data
  call allocate_io_writer(adios2_writer)
  call adios2_writer%init(MPI_COMM_WORLD, "test_io_write")
  file = adios2_writer%open("test_output.bp", io_mode_write, MPI_COMM_WORLD)
  call file%begin_step()
  call adios2_writer%write_data("data3D", data_write, file, &
                                shape_dims, start_dims, count_dims)
  call file%end_step()
  call file%close()

  if (allocated(data_write)) deallocate (data_write)
  call adios2_writer%finalise()

  ! read data (rank 0 only)
  if (irank == 0) then
    call allocate_io_reader(adios2_reader)
    call adios2_reader%init(MPI_COMM_SELF, "test_io_read")
    ! Note: file is automatically deallocated when going out of scope
    file = adios2_reader%open("test_output.bp", io_mode_read, MPI_COMM_SELF)
    call file%begin_step()

    sel_start = [0, 0, 0]
    sel_count = [shape_dims(1), shape_dims(2), shape_dims(3)]

    ! read entire dataset
    allocate (data_read(shape_dims(1), shape_dims(2), shape_dims(3)))
    call adios2_reader%read_data("data3D", data_read, file, &
                                 start_dims=sel_start, count_dims=sel_count)

    if (is_sp) then
      tolerance = 1.e-6_dp
    else
      tolerance = 1.e-8_dp
    end if

    ! verify data
    do k = 1, sel_count(3)
      do j = 1, sel_count(2)
        do i = 1, sel_count(1)
          ! calculate expected value based on original data pattern
          rank_id = (i - 1)/inx ! determine which rank wrote this data
          expected = rank_id*inx*iny*inz + (k - 1)*inx*iny + &
                     (j - 1)*inx + mod(i - 1, inx)
          if (abs(data_read(i, j, k) - expected) > tolerance) then
            allpass = .false.
            print *, "Data mismatch at (", i, ",", j, ",", k, "): ", &
              data_read(i, j, k), " expected: ", expected
          end if
        end do
      end do
    end do

    call file%end_step()
    call file%close()

    if (allocated(data_read)) deallocate (data_read)
    call adios2_reader%finalise()
    if (allocated(adios2_reader)) deallocate(adios2_reader)
  end if

  ! cleanup writer
  if (allocated(adios2_writer)) deallocate(adios2_writer)

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
