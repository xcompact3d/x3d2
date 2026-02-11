program test_cuda_gpu_aware_io
  !! Simplified test for CUDA backend copy functions and ADIOS2 I/O
  !! Tests:
  !! 1. Backend copy_data_to_f (host -> device)
  !! 2. Backend copy_f_to_data (device -> host)
  !! 3. ADIOS2 write/read with host arrays (like test_adios2_read_write)
  use mpi
  use m_common, only: dp, i8, DIR_C
  use m_allocator, only: field_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_io_backend, only: allocate_io_writer, allocate_io_reader
  use m_io_base, only: io_writer_t, io_reader_t, io_file_t, &
                       io_mode_write, io_mode_read
  use iso_fortran_env, only: stderr => error_unit
  implicit none

  ! Test variables
  type(cuda_allocator_t) :: allocator
  type(cuda_backend_t) :: backend
  class(field_t), pointer :: cuda_field
  class(io_writer_t), allocatable :: writer
  class(io_reader_t), allocatable :: reader
  class(io_file_t), allocatable :: file
  
  ! MPI variables
  integer :: ierr, irank, isize
  integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
  integer(i8), dimension(3) :: sel_start, sel_count
  real(dp), dimension(:, :, :), allocatable :: data_write, data_from_device, data_read
  
  ! Test parameters
  integer :: i, j, k, rank_id
  integer, dimension(3) :: field_shape
  logical :: allpass = .true.
  real(dp) :: expected, tolerance
  
  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)
  
  ! Initialize allocator and backend with mesh dimensions
  allocator = cuda_allocator_t([16, 16, 16], 2)
  
  ! Test 1: Backend copy to/from device
  ! Get a CUDA field and use its shape for the data array
  cuda_field => allocator%get_block(DIR_C)
  field_shape = cuda_field%get_shape()
  
  if (irank == 0) then
    write(stderr, '(a,3i4)') 'Field dimensions (DIR_C): ', field_shape
  end if
  
  allocate(data_write(field_shape(1), field_shape(2), field_shape(3)))
  
  ! Initialize data with unique pattern per rank
  do k = 1, field_shape(3)
    do j = 1, field_shape(2)
      do i = 1, field_shape(1)
        data_write(i, j, k) = real(irank*1000 + (k-1)*field_shape(1)*field_shape(2) &
                                   + (j-1)*field_shape(1) + (i-1), dp)
      end do
    end do
  end do
  
  ! Copy host -> device
  call backend%copy_data_to_f(cuda_field, data_write)
  
  ! Copy device -> host to verify round-trip
  allocate(data_from_device(field_shape(1), field_shape(2), field_shape(3)))
  call backend%copy_f_to_data(data_from_device, cuda_field)
  
  ! Verify backend copy operations
  tolerance = 1.0e-12_dp
  do k = 1, field_shape(3)
    do j = 1, field_shape(2)
      do i = 1, field_shape(1)
        if (abs(data_from_device(i,j,k) - data_write(i,j,k)) > tolerance) then
          allpass = .false.
          write(stderr, '(a,i2,a,3i4,a,f12.4,a,f12.4)') &
            'Rank ', irank, ' ERROR: Backend copy mismatch at (', i, j, k, '): ', &
            data_from_device(i, j, k), ' expected: ', data_write(i, j, k)
        end if
      end do
    end do
  end do
  
  if (irank == 0 .and. allpass) then
    write(stderr, '(a)') 'Backend copy test: PASSED'
  end if
  
  ! Test 2: ADIOS2 I/O (using simple write_data like test_adios2_read_write)
  shape_dims = [int(isize*field_shape(1), i8), int(field_shape(2), i8), int(field_shape(3), i8)]
  start_dims = [int(irank*field_shape(1), i8), 0_i8, 0_i8]
  count_dims = [int(field_shape(1), i8), int(field_shape(2), i8), int(field_shape(3), i8)]
  
  ! Write data
  call allocate_io_writer(writer)
  call writer%init(MPI_COMM_WORLD, "test_cuda_io")
  file = writer%open("test_cuda_output.bp", io_mode_write, MPI_COMM_WORLD)
  call file%begin_step()
  call writer%write_data("data3D", data_write, file, &
                         shape_dims, start_dims, count_dims)
  call file%end_step()
  call file%close()
  call writer%finalise()
  
  deallocate(data_write, data_from_device)
  call allocator%release_block(cuda_field)
  call allocator%destroy()
  
  ! Read and verify (rank 0 only)
  if (irank == 0) then
    call allocate_io_reader(reader)
    call reader%init(MPI_COMM_SELF, "test_cuda_io_read")
    file = reader%open("test_cuda_output.bp", io_mode_read, MPI_COMM_SELF)
    call file%begin_step()
    
    sel_start = [0_i8, 0_i8, 0_i8]
    sel_count = [shape_dims(1), shape_dims(2), shape_dims(3)]
    
    allocate(data_read(sel_count(1), sel_count(2), sel_count(3)))
    call reader%read_data("data3D", data_read, file, &
                         start_dims=sel_start, count_dims=sel_count)
    
    call file%end_step()
    call file%close()
    call reader%finalise()
    
    ! Verify data - calculate expected values based on actual field_shape
    do k = 1, sel_count(3)
      do j = 1, sel_count(2)
        do i = 1, sel_count(1)
          rank_id = (i-1)/field_shape(1)
          expected = real(rank_id*1000 + (k-1)*field_shape(1)*field_shape(2) + &
                         (j-1)*field_shape(1) + mod(i-1, field_shape(1)), dp)
          if (abs(data_read(i, j, k) - expected) > tolerance) then
            allpass = .false.
            write(stderr, '(a,3i4,a,f12.4,a,f12.4)') &
              'ERROR: ADIOS2 I/O mismatch at (', i, j, k, '): ', &
              data_read(i, j, k), ' expected: ', expected
          end if
        end do
      end do
    end do
    
    if (allpass) then
      write(stderr, '(a)') 'ADIOS2 I/O test: PASSED'
    end if
    
    deallocate(data_read)
  end if
  
  ! Finalize
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)
  
  ! Report final results
  if (irank == 0) then
    if (allpass) then
      write(stderr, '(a)') 'CUDA I/O TEST PASSED SUCCESSFULLY.'
    else
      error stop 'CUDA I/O TEST FAILED.'
    end if
  end if

end program test_cuda_gpu_aware_io
