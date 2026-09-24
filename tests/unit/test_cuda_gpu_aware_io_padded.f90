program test_cuda_gpu_aware_io_padded
  !! Regression test for GPU-aware ADIOS2 I/O on a DIR_X solver field whose
  !! padded storage does not match its true (nx, ny, nz) extent (global
  !! dims are not multiples of SZ in x and y). This is exactly the shape
  !! snapshot_manager/checkpoint_manager hand to write_field_from_solver:
  !! a raw DIR_X field, not the small all-multiples-of-SZ DIR_C field used
  !! by test_cuda_gpu_aware_io.f90.
  !!
  !! Writes through the same write_field_from_solver entry the snapshot
  !! manager uses, reads the result back, and compares against
  !! backend%get_field_data + slicing to (1:nx, 1:ny, 1:nz). Before the
  !! reorder/pack fix this fails (the raw padded field was put as if it
  !! were a contiguous Cartesian array); after the fix it passes exactly.
  !!
  !! It also covers staging buffer reuse: DIR_X and DIR_C fields in both
  !! precisions, several fields per step, and two steps. Each write gets
  !! new values in a field an earlier write used, so a write that read the
  !! shared staging buffer too late would see the wrong values.
  !! An optional first argument overrides the output file name.
  use mpi
  use m_common, only: dp, i8, DIR_C, DIR_X, VERT, is_sp
  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t
  use m_mesh, only: mesh_t
  use m_io_backend, only: allocate_io_writer, allocate_io_reader
  use m_io_base, only: io_writer_t, io_reader_t, io_file_t, &
                       io_mode_write, io_mode_read
  use m_test_utils, only: initialise_mpi, finalise_test
  use iso_fortran_env, only: stderr => error_unit
  implicit none

  integer, parameter :: n_steps = 2, n_writes = 10

  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend
  type(mesh_t), target :: mesh
  class(field_t), pointer :: u_x, u_c, field
  class(io_writer_t), allocatable :: writer
  class(io_reader_t), allocatable :: reader
  class(io_file_t), allocatable :: file

  integer :: ierr, irank, isize
  integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
  real(dp), dimension(:, :, :), allocatable :: data_cart, data_read
  real(dp), dimension(:, :, :, :), allocatable :: ref_full

  integer, dimension(3) :: dims_global, nproc_dir, dims_padded_c, dims_local
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  character(len=256) :: filename
  character(len=16) :: write_names(n_writes)
  logical :: write_sp(n_writes), write_dir_x(n_writes)
  integer :: nx, ny, nz, i, j, k, istep, iwrite, arg_len
  logical :: allpass = .true.
  real(dp) :: tolerance, expected

  call initialise_mpi(irank, isize)

  filename = "test_cuda_gpu_padded_output.bp"
  call get_command_argument(1, filename, length=arg_len)
  if (arg_len == 0) filename = "test_cuda_gpu_padded_output.bp"

  ! Global dims are deliberately not multiples of SZ (32) in x and y, to
  ! exercise the padding that the bug ignored.
  dims_global = [45, 38, 20]
  nproc_dir = [1, 1, isize]
  L_global = [1.0_dp, 1.0_dp, 1.0_dp]
  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)
  call runtime%init(mesh)
  backend => runtime%backend

  dims_local = mesh%get_dims(VERT)
  nx = dims_local(1); ny = dims_local(2); nz = dims_local(3)

  u_x => runtime%allocator%get_block(DIR_X)
  u_c => runtime%allocator%get_block(DIR_C)

  dims_padded_c = runtime%allocator%get_padded_dims(DIR_C)
  allocate (data_cart(dims_padded_c(1), dims_padded_c(2), dims_padded_c(3)))
  allocate (ref_full(dims_padded_c(1), dims_padded_c(2), dims_padded_c(3), &
                     n_steps))

  ! Alternate DIR_X/DIR_C sources and precisions over more writes than the
  ! initial staging pool (8 buffers) holds.
  do iwrite = 1, n_writes
    write (write_names(iwrite), '(A,I0)') 'field_', iwrite
    write_dir_x(iwrite) = mod(iwrite, 2) == 1
    write_sp(iwrite) = mod((iwrite - 1)/2, 2) == 1
  end do

  shape_dims = [int(nx, i8), int(ny, i8), int(nz, i8)*int(isize, i8)]
  start_dims = [0_i8, 0_i8, int(irank, i8)*int(nz, i8)]
  count_dims = [int(nx, i8), int(ny, i8), int(nz, i8)]

  call allocate_io_writer(writer)

  if (.not. writer%supports_device_field_write(u_x, backend)) then
    if (irank == 0) then
      write (stderr, '(a)') 'GPU-aware ADIOS2 not available - skipping test'
    end if
    call runtime%allocator%release_block(u_x)
    call runtime%allocator%release_block(u_c)
    call MPI_Finalize(ierr)
    stop
  end if

  call writer%init(MPI_COMM_WORLD, "test_cuda_gpu_aware_padded")
  file = writer%open(trim(filename), io_mode_write, MPI_COMM_WORLD)

  do istep = 1, n_steps
    ! New data every step, so a staging buffer reused too early would show
    call random_number(data_cart)
    call backend%set_field_data(u_x, data_cart, DIR_C)

    ! Reference: the same data read back through the host path
    call backend%get_field_data(ref_full(:, :, :, istep), u_x, DIR_C)

    call file%begin_step()
    do iwrite = 1, n_writes
      if (write_dir_x(iwrite)) then
        field => u_x
      else
        field => u_c
      end if
      ! Write iwrite holds ref + (iwrite - 1)
      call backend%set_field_data(field, data_cart + real(iwrite - 1, dp), &
                                  DIR_C)
      call writer%write_field_from_solver( &
        trim(write_names(iwrite)), field, file, backend, &
        shape_dims, start_dims, count_dims, use_sp=write_sp(iwrite))
    end do
    call file%end_step()
  end do

  call file%close()
  call writer%finalise()

  ! Verify (rank 0 only, reading its own local block back)
  if (irank == 0) then
    call allocate_io_reader(reader)
    call reader%init(MPI_COMM_SELF, "test_cuda_gpu_read_padded")
    file = reader%open(trim(filename), io_mode_read, MPI_COMM_SELF)
    allocate (data_read(nx, ny, nz))

    do istep = 1, n_steps
      call file%begin_step()
      do iwrite = 1, n_writes
        call reader%read_data(trim(write_names(iwrite)), data_read, file, &
                              start_dims=start_dims, count_dims=count_dims)

        if (is_sp .or. write_sp(iwrite)) then
          tolerance = 1.0e-5_dp
        else
          tolerance = 1.0e-12_dp
        end if

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              expected = ref_full(i, j, k, istep) + real(iwrite - 1, dp)
              if (abs(data_read(i, j, k) - expected) > tolerance) then
                allpass = .false.
                write (stderr, '(a,i0,a,3i4,a,f14.6,a,f14.6)') &
                  'ERROR: step ', istep, ' '//trim(write_names(iwrite))// &
                  ' mismatch at (', i, j, k, '): ', data_read(i, j, k), &
                  ' expected: ', expected
              end if
            end do
          end do
        end do
      end do
      call file%end_step()
    end do

    call file%close()
    call reader%finalise()
    deallocate (data_read)
  end if

  call runtime%allocator%release_block(u_x)
  call runtime%allocator%release_block(u_c)

  call finalise_test(allpass, irank)

end program test_cuda_gpu_aware_io_padded
