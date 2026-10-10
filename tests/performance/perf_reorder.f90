program perf_reorder
  !! Benchmark of the data reorderings between the x, y, z pencil and the
  !! Cartesian layouts (backend%reorder).
  use m_mpi, only: MPI_Wtime

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_device, &
                               backend_label
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      RDR_X2C, RDR_C2X
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test, &
                          write_perf_metric, write_perf_summary, &
                          write_device_bw_metric

  implicit none

  integer, parameter :: n_iters = 100, n_warmup = 10
  real(dp), parameter :: consumed_bw = 2.0_dp
  integer :: n_glob, ndof, nrank, nproc
  integer :: mem_clock_rt, mem_bus_width
  logical :: has_device_bw_info

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(field_t), pointer :: u_x, u_y, u_z, u_c, u_out
  real(dp), allocatable :: data(:, :, :)
  integer :: dims(3)

  call initialise_mpi(nrank, nproc)

  if (backend_is_device) then
    n_glob = 512
  else
    n_glob = 256
  end if
  ndof = n_glob**3

  mesh = mesh_t([n_glob, n_glob, n_glob], [1, 1, 1], &
                [1._dp, 1._dp, 1._dp], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])
  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator
  call backend%get_device_bw_info(mem_clock_rt, mem_bus_width, &
                                  has_device_bw_info)

  u_x => allocator%get_block(DIR_X)
  u_y => allocator%get_block(DIR_Y)
  u_z => allocator%get_block(DIR_Z)
  u_c => allocator%get_block(DIR_C)

  dims = u_x%get_shape()
  allocate (data(dims(1), dims(2), dims(3)))
  call random_number(data)
  call backend%set_field_data(u_x, data, DIR_X)
  call backend%reorder(u_y, u_x, RDR_X2Y)
  call backend%reorder(u_z, u_x, RDR_X2Z)
  call backend%reorder(u_c, u_x, RDR_X2C)

  call run_case('x2y', DIR_Y, u_x, RDR_X2Y)
  call run_case('x2z', DIR_Z, u_x, RDR_X2Z)
  call run_case('y2x', DIR_X, u_y, RDR_Y2X)
  call run_case('y2z', DIR_Z, u_y, RDR_Y2Z)
  call run_case('z2x', DIR_X, u_z, RDR_Z2X)
  call run_case('z2y', DIR_Y, u_z, RDR_Z2Y)
  call run_case('x2c', DIR_C, u_x, RDR_X2C)
  call run_case('c2x', DIR_X, u_c, RDR_C2X)

  if (has_device_bw_info) then
    call write_device_bw_metric(mem_clock_rt, mem_bus_width)
  end if

  call allocator%release_block(u_x)
  call allocator%release_block(u_y)
  call allocator%release_block(u_z)
  call allocator%release_block(u_c)
  call finalise_test(.true., nrank)

contains

  subroutine run_case(case_name, dir_out, u_in, rdr)
    character(len=*), intent(in) :: case_name
    integer, intent(in) :: dir_out, rdr
    class(field_t), intent(in) :: u_in

    integer :: iter
    real(dp) :: tstart, tend
    character(len=:), allocatable :: label

    label = trim(backend_label)//'_reorder_'//case_name
    if (nrank == 0) print *, 'Performance test:', label

    u_out => allocator%get_block(dir_out)

    do iter = 1, n_warmup
      call backend%reorder(u_out, u_in, rdr)
    end do
    call backend%sync()

    tstart = MPI_Wtime()
    do iter = 1, n_iters
      call backend%reorder(u_out, u_in, rdr)
    end do
    call backend%sync()
    tend = MPI_Wtime()

    call allocator%release_block(u_out)

    call write_perf_metric(label, tend - tstart, n_iters, ndof, consumed_bw)
    if (has_device_bw_info) then
      call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw, &
                              mem_clock_rt, mem_bus_width)
    else
      call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw)
    end if
  end subroutine run_case

end program perf_reorder
