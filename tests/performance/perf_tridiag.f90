program perf_tridiag
  !! Benchmark of the distributed tridiagonal solver (backend%tds_solve),
  !! including the halo exchange, on a domain decomposed in x.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier, MPI_Wtime

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_device, &
                               backend_label, backend_sz
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, BC_PERIODIC, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test, report_perf_minmax

  implicit none

  integer, parameter :: n_glob = 1024, n_warmup = 10
  integer :: n, n_groups, n_iters, ndof
  integer :: nrank, nproc, ierr
  integer :: mem_clock_rt, mem_bus_width
  logical :: has_device_bw_info
  real(dp) :: periodic_bw

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(field_t), pointer :: u_field, du_field
  real(dp), allocatable :: u_data(:, :, :)
  ! Program scope, so it is never deallocated: NVHPC (25.3) crashes in
  ! pgf90_dealloc_poly03 when a class(tdsops_t) holding a cuda_tdsops_t is.
  class(tdsops_t), allocatable :: tdsops

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Performance benchmark with', nproc, 'ranks'

  if (backend_is_device) then
    n_groups = 512*512/backend_sz
    n_iters = 100
  else
    n_groups = 64*64/backend_sz
    n_iters = 1000
  end if
  periodic_bw = 5.0_dp

  mesh = mesh_t([n_glob, backend_sz, n_groups], [nproc, 1, 1], &
                [2*pi, 1._dp, 1._dp], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])
  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator

  n = mesh%get_n(DIR_X, VERT)
  ndof = n*n_groups*backend_sz

  u_field => allocator%get_block(DIR_X)
  du_field => allocator%get_block(DIR_X)

  call backend%get_device_bw_info(mem_clock_rt, mem_bus_width, &
                                  has_device_bw_info)

  call run_case('periodic', 2*pi/n_glob, periodic_bw)

  call allocator%release_block(u_field)
  call allocator%release_block(du_field)
  call finalise_test(.true., nrank)

contains

  subroutine run_case(case_name, delta, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: delta, consumed_bw

    integer :: iter, j, dims(3)
    real(dp) :: tstart, tend
    character(len=:), allocatable :: label

    dims = u_field%get_shape()
    allocate (u_data(dims(1), dims(2), dims(3)))
    u_data = 0._dp
    do j = 1, n
      u_data(:, j, :) = sin((j - 1 + mesh%par%nrank_dir(DIR_X)*n)*delta)
    end do
    call backend%set_field_data(u_field, u_data, DIR_X)

    call backend%alloc_tdsops(tdsops, n, delta, operation='second-deriv', &
                              scheme='compact6', &
                              bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    do iter = 1, n_warmup
      call backend%tds_solve(du_field, u_field, tdsops)
    end do
    call backend%sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    tstart = MPI_Wtime()
    do iter = 1, n_iters
      call backend%tds_solve(du_field, u_field, tdsops)
    end do
    call backend%sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    tend = MPI_Wtime()

    label = trim(backend_label)//'_tridiag_'//case_name
    if (has_device_bw_info) then
      call report_perf_minmax(label, tend - tstart, n_iters, ndof, &
                              consumed_bw, mem_clock_rt, mem_bus_width)
    else
      call report_perf_minmax(label, tend - tstart, n_iters, ndof, &
                              consumed_bw)
    end if

    deallocate (u_data)
  end subroutine run_case

end program perf_tridiag
