program perf_cuda_thom

  use cudafor
  use m_common, only: dp, nbytes, pi, BC_PERIODIC, BC_DIRICHLET
  use m_cuda_common, only: SZ
  use m_cuda_exec_thom, only: exec_thom_tds_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

  implicit none

  integer :: i, j, k
  integer :: n_glob, n, n_block
  integer :: n_iters, n_warmup
  integer :: ndof, ierr

  real(dp) :: tstart, tend
  real(dp), allocatable, dimension(:, :, :) :: u
  real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
  real(dp) :: dx_per, dx
  real(dp) :: achievedBW
  integer :: memClockRt, memBusWidth

  type(cuda_tdsops_t) :: tdsops
  type(dim3) :: blocks, threads

  ! Performance configuration
  n_glob = 512
  n_block = 512*512/SZ
  n_iters = 1000
  n_warmup = 10
  n = n_glob
  ndof = n_glob*n_block*SZ

  allocate (u(SZ, n, n_block))
  allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  blocks = dim3(n_block, 1, 1)
  threads = dim3(SZ, 1, 1)

  ! Get device BW info
  ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
  ierr = cudaDeviceGetAttribute(memBusWidth, &
                                cudaDevAttrGlobalMemoryBusWidth, 0)

  ! Periodic case benchmark
  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx_per)
      end do
    end do
  end do
  u_dev = u

  tdsops = cuda_tdsops_init(n, dx_per, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  ! Warmup
  do i = 1, n_warmup
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  ierr = cudaDeviceSynchronize()

  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  ierr = cudaDeviceSynchronize()
  call cpu_time(tend)

  achievedBW = 6.0_dp*n_iters*ndof*nbytes/(tend - tstart)
  print '(a, f10.6, a, f10.3, a)', &
    'PERF_METRIC: cuda_thom_periodic time=', tend - tstart, &
    's bw=', achievedBW/real(2**30, dp), ' GiB/s'

  ! Dirichlet case benchmark
  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx)
      end do
    end do
  end do
  u_dev = u

  tdsops = cuda_tdsops_init(n, dx, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

  ! Warmup
  do i = 1, n_warmup
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  ierr = cudaDeviceSynchronize()

  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  ierr = cudaDeviceSynchronize()
  call cpu_time(tend)

  achievedBW = 4.0_dp*n_iters*ndof*nbytes/(tend - tstart)
  print '(a, f10.6, a, f10.3, a)', &
    'PERF_METRIC: cuda_thom_dirichlet time=', tend - tstart, &
    's bw=', achievedBW/real(2**30, dp), ' GiB/s'

  ! Device BW reference
  print '(a, f10.3, a)', 'PERF_METRIC: device_bw ref=', &
    2.0_dp*memBusWidth/8._dp*memClockRt*1000/real(2**30, dp), ' GiB/s'

end program perf_cuda_thom
