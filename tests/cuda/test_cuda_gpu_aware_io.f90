program test_cuda_gpu_aware_io
  !! Test for GPU-aware ADIOS2 I/O with CUDA backend
  !! Tests:
  !! 1. Write directly from device memory using GPU-aware ADIOS2 path
  !! 2. Verify data integrity by reading back and comparing
  use mpi
  use m_common, only: dp, i8, DIR_C
  use m_allocator, only: field_t
  use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
  use m_io_backend, only: allocate_io_writer, allocate_io_reader
  use m_io_base, only: io_writer_t, io_reader_t, io_file_t, &
                       io_mode_write, io_mode_read
  use iso_fortran_env, only: stderr => error_unit
  implicit none

  ! Test variables
  type(cuda_allocator_t) :: allocator
  class(field_t), pointer :: cuda_field
  class(io_writer_t), allocatable :: writer
  class(io_reader_t), allocatable :: reader
  class(io_file_t), allocatable :: file

  ! MPI variables
  integer :: ierr, irank, isize, dummy_backend
  integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
  integer(i8), dimension(3) :: sel_start, sel_count
  real(dp), dimension(:, :, :), allocatable :: data_write, data_read

  ! Test parameters
  integer :: i, j, k, rank_id
  integer, dimension(3) :: field_shape
  logical :: allpass = .true.
  real(dp) :: expected, tolerance

  ! Initialise MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

  ! Initialise allocator with mesh dimensions
  allocator = cuda_allocator_t([16, 16, 16], 2)

  ! Get a CUDA field and populate it via device assignment
  cuda_field => allocator%get_block(DIR_C)
  field_shape = cuda_field%get_shape()

  if (irank == 0) then
    write(stderr, '(a,3i4)') 'Field dimensions (DIR_C): ', field_shape
  end if

  allocate(data_write(field_shape(1), field_shape(2), field_shape(3)))

  ! Initialise data with unique pattern per rank
  do k = 1, field_shape(3)
    do j = 1, field_shape(2)
      do i = 1, field_shape(1)
        data_write(i, j, k) = real(irank*1000 + (k-1)*field_shape(1)*field_shape(2) &
                                   + (j-1)*field_shape(1) + (i-1), dp)
      end do
    end do
  end do

  ! Copy host -> device directly via cudafor assignment
  select type (cuda_field)
  type is (cuda_field_t)
    cuda_field%data_d = data_write
  end select

  ! GPU-aware ADIOS2 I/O - write directly from device memory
  shape_dims = [int(isize*field_shape(1), i8), int(field_shape(2), i8), int(field_shape(3), i8)]
  start_dims = [int(irank*field_shape(1), i8), 0_i8, 0_i8]
  count_dims = [int(field_shape(1), i8), int(field_shape(2), i8), int(field_shape(3), i8)]

  if (irank == 0) then
    write(stderr, '(a)') 'Testing GPU-aware write from device memory...'
  end if

  ! Write from device using GPU-aware path
  ! dummy_backend is passed as class(*) — the GPU path dispatches on field type
  dummy_backend = 0
  call allocate_io_writer(writer)
  call writer%init(MPI_COMM_WORLD, "test_cuda_gpu_aware")
  file = writer%open("test_cuda_gpu_output.bp", io_mode_write, MPI_COMM_WORLD)
  call file%begin_step()

  call writer%write_field_from_solver("velocity_x", cuda_field, file, dummy_backend, &
                                      shape_dims, start_dims, count_dims)

  call file%end_step()
  call file%close()
  call writer%finalise()

  if (irank == 0) then
    write(stderr, '(a)') 'GPU-aware write completed'
  end if

  deallocate(data_write)

  ! Read and verify (rank 0 only)
  if (irank == 0) then
    write(stderr, '(a)') 'Verifying GPU-aware write...'
    call allocate_io_reader(reader)
    call reader%init(MPI_COMM_SELF, "test_cuda_gpu_read")
    file = reader%open("test_cuda_gpu_output.bp", io_mode_read, MPI_COMM_SELF)
    call file%begin_step()

    sel_start = [0_i8, 0_i8, 0_i8]
    sel_count = [shape_dims(1), shape_dims(2), shape_dims(3)]

    allocate(data_read(sel_count(1), sel_count(2), sel_count(3)))
    call reader%read_data("velocity_x", data_read, file, &
                         start_dims=sel_start, count_dims=sel_count)

    call file%end_step()
    call file%close()
    call reader%finalise()

    ! Verify data
    tolerance = 1.0e-12_dp
    do k = 1, sel_count(3)
      do j = 1, sel_count(2)
        do i = 1, sel_count(1)
          rank_id = (i-1)/field_shape(1)
          expected = real(rank_id*1000 + (k-1)*field_shape(1)*field_shape(2) + &
                         (j-1)*field_shape(1) + mod(i-1, field_shape(1)), dp)
          if (abs(data_read(i, j, k) - expected) > tolerance) then
            allpass = .false.
            write(stderr, '(a,3i4,a,f12.4,a,f12.4)') &
              'ERROR: GPU-aware I/O mismatch at (', i, j, k, '): ', &
              data_read(i, j, k), ' expected: ', expected
          end if
        end do
      end do
    end do

    deallocate(data_read)
  end if

  ! Cleanup
  call allocator%release_block(cuda_field)
  call allocator%destroy()

  ! Finalise
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)

  ! Report final results
  if (irank == 0) then
    if (allpass) then
      write(stderr, '(a)') 'CUDA GPU-AWARE I/O TEST PASSED!'
    else
      error stop 'CUDA GPU-AWARE I/O TEST FAILED.'
    end if
  end if

end program test_cuda_gpu_aware_io
