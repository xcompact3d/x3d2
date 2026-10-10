program perf_penta
  !! Benchmark of the compact10_penta (Lele 10th-order pentadiagonal)
  !! first-derivative solve: alpha=0.5, beta=0.05, non-periodic, each rank
  !! solving its own lines with zero halos.
  use m_backend_runtime, only: backend_runtime_t, backend_is_device, &
                               backend_label, backend_sz, penta_benchmark
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, BC_PERIODIC
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test, report_perf_minmax

  implicit none

  ! Lele 10 penta 1st-deriv: reads 7 values, writes 1 -> bandwidth factor 8
  real(dp), parameter :: penta_bw = 8.0_dp
  integer, parameter :: n_glob = 1024, n_halo = 4, n_warmup = 10
  integer :: n, n_groups, n_iters, ndof, nrank, nproc, j
  integer :: mem_clock_rt, mem_bus_width
  logical :: has_device_bw_info
  real(dp) :: dx, elapsed
  character(len=:), allocatable :: label

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(tdsops_t), allocatable :: tdsops
  real(dp), allocatable :: u(:, :, :), u_recv_s(:, :, :), u_recv_e(:, :, :)

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) then
    print *, 'Performance benchmark for compact10_penta (Lele 10 penta 1st-deriv)'
    print *, 'Ranks:', nproc
  end if

  if (backend_is_device) then
    n_groups = 512*512/backend_sz
    n_iters = 100
  else
    n_groups = 64*64/backend_sz
    n_iters = 100
  end if

  ! The backend is only needed to build backend-specific tdsops
  mesh = mesh_t([backend_sz, backend_sz, backend_sz], [1, 1, 1], &
                [1._dp, 1._dp, 1._dp], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])
  call runtime%init(mesh)
  backend => runtime%backend
  call backend%get_device_bw_info(mem_clock_rt, mem_bus_width, &
                                  has_device_bw_info)

  n = n_glob/nproc
  ndof = n*n_groups*backend_sz
  dx = 2*pi/n_glob

  allocate (u(backend_sz, n, n_groups))
  allocate (u_recv_s(backend_sz, n_halo, n_groups))
  allocate (u_recv_e(backend_sz, n_halo, n_groups))
  do j = 1, n
    u(:, j, :) = sin((j - 1 + nrank*n)*dx)
  end do
  u_recv_s = 0._dp
  u_recv_e = 0._dp

  ! Periodic BCs only set alpha/beta/coeffs; the non-periodic kernel is used
  call backend%alloc_tdsops(tdsops, n, dx, operation='first-deriv', &
                            scheme='compact10_penta', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  elapsed = penta_benchmark(u, u_recv_s, u_recv_e, tdsops, .false., &
                            n_warmup, n_iters)

  label = trim(backend_label)//'_penta_penta_lele10_nonper'
  if (has_device_bw_info) then
    call report_perf_minmax(label, elapsed, n_iters, ndof, penta_bw, &
                            mem_clock_rt, mem_bus_width)
  else
    call report_perf_minmax(label, elapsed, n_iters, ndof, penta_bw)
  end if

  call finalise_test(.true., nrank)

end program perf_penta
