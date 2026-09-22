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

  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend
  type(mesh_t), target :: mesh
  class(field_t), pointer :: u_x
  class(io_writer_t), allocatable :: writer
  class(io_reader_t), allocatable :: reader
  class(io_file_t), allocatable :: file

  integer :: ierr, irank, isize
  integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
  real(dp), dimension(:, :, :), allocatable :: data_cart, ref_full, data_read

  integer, dimension(3) :: dims_global, nproc_dir, dims_padded_c, dims_local
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  integer :: nx, ny, nz, i, j, k
  logical :: allpass = .true.
  real(dp) :: tolerance

  call initialise_mpi(irank, isize)

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

  dims_padded_c = runtime%allocator%get_padded_dims(DIR_C)
  allocate (data_cart(dims_padded_c(1), dims_padded_c(2), dims_padded_c(3)))
  allocate (ref_full(dims_padded_c(1), dims_padded_c(2), dims_padded_c(3)))

  call random_number(data_cart)
  call backend%set_field_data(u_x, data_cart, DIR_C)

  ! Reference: the same reorder the fix uses (DIR_X -> DIR_C), sliced to
  ! the true field extent.
  call backend%get_field_data(ref_full, u_x, DIR_C)

  shape_dims = [int(nx, i8), int(ny, i8), int(nz, i8)*int(isize, i8)]
  start_dims = [0_i8, 0_i8, int(irank, i8)*int(nz, i8)]
  count_dims = [int(nx, i8), int(ny, i8), int(nz, i8)]

  call allocate_io_writer(writer)

  if (.not. writer%supports_device_field_write(u_x, backend)) then
    if (irank == 0) then
      write (stderr, '(a)') 'GPU-aware ADIOS2 not available - skipping test'
    end if
    call runtime%allocator%release_block(u_x)
    call MPI_Finalize(ierr)
    stop
  end if

  call writer%init(MPI_COMM_WORLD, "test_cuda_gpu_aware_padded")
  file = writer%open("test_cuda_gpu_padded_output.bp", io_mode_write, &
                     MPI_COMM_WORLD)
  call file%begin_step()

  call writer%write_field_from_solver("field_x", u_x, file, backend, &
                                      shape_dims, start_dims, count_dims)

  call file%end_step()
  call file%close()
  call writer%finalise()

  ! Verify (rank 0 only, reading its own local block back)
  if (irank == 0) then
    call allocate_io_reader(reader)
    call reader%init(MPI_COMM_SELF, "test_cuda_gpu_read_padded")
    file = reader%open("test_cuda_gpu_padded_output.bp", io_mode_read, &
                       MPI_COMM_SELF)
    call file%begin_step()

    allocate (data_read(nx, ny, nz))
    call reader%read_data("field_x", data_read, file, &
                          start_dims=start_dims, count_dims=count_dims)

    call file%end_step()
    call file%close()
    call reader%finalise()

    if (is_sp) then
      tolerance = 1.0e-5_dp
    else
      tolerance = 1.0e-12_dp
    end if

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (abs(data_read(i, j, k) - ref_full(i, j, k)) > tolerance) then
            allpass = .false.
            write (stderr, '(a,3i4,a,f14.6,a,f14.6)') &
              'ERROR: GPU-aware padded I/O mismatch at (', i, j, k, '): ', &
              data_read(i, j, k), ' expected: ', ref_full(i, j, k)
          end if
        end do
      end do
    end do

    deallocate (data_read)
  end if

  call runtime%allocator%release_block(u_x)

  call finalise_test(allpass, irank)

end program test_cuda_gpu_aware_io_padded
