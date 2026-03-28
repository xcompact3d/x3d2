program perf_thom

#ifdef CUDA
  use cudafor
#else
  use omp_lib
#endif
  use m_common, only: dp, pi, BC_PERIODIC, BC_DIRICHLET
#ifdef CUDA
  use m_cuda_common, only: SZ
  use m_cuda_exec_thom, only: exec_thom_tds_compact
  use m_cuda_tdsops, only: tdsops_t => cuda_tdsops_t, tdsops_init => cuda_tdsops_init
#else
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_exec_thom, only: exec_thom_tds_compact
#endif
  use m_test_utils, only: write_perf_metric, write_perf_summary, write_device_bw_metric

  implicit none

#ifdef CUDA
  character(len=*), parameter :: backend_label = 'cuda'
  real(dp), parameter :: periodic_bw = 6.0_dp
  real(dp), parameter :: dirichlet_bw = 4.0_dp
#else
  character(len=*), parameter :: backend_label = 'omp'
  real(dp), parameter :: periodic_bw = 3.0_dp
  real(dp), parameter :: dirichlet_bw = 3.0_dp
#endif

  integer :: i, j, k
  integer :: n_glob, n, n_groups
  integer :: n_iters, n_warmup
  integer :: ndof
  real(dp) :: dx_per, dx

#ifdef CUDA
  integer :: ierr
  integer :: memClockRt, memBusWidth
  real(dp), allocatable :: u(:, :, :)
  real(dp), device, allocatable :: u_dev(:, :, :), du_dev(:, :, :)
  type(dim3) :: blocks, threads
#else
  real(dp), allocatable :: u(:, :, :), du(:, :, :)
#endif

  type(tdsops_t) :: tdsops

  call configure_benchmark()
  call allocate_fields()
  call setup_backend()

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  call run_case('periodic', dx_per, BC_PERIODIC, BC_PERIODIC, periodic_bw)
  call run_case('dirichlet', dx, BC_DIRICHLET, BC_DIRICHLET, dirichlet_bw)

#ifdef CUDA
  call write_device_bw_metric(memClockRt, memBusWidth)
#endif

contains

  subroutine configure_benchmark()
    n_glob = 512
#ifdef CUDA
    n_groups = 512*512/SZ
    n_iters = 1000
#else
    n_groups = 128*128/SZ
    n_iters = 500
#endif
    n_warmup = 10
    n = n_glob
    ndof = n_glob*n_groups*SZ
  end subroutine configure_benchmark

  subroutine allocate_fields()
#ifdef CUDA
    allocate (u(SZ, n, n_groups))
    allocate (u_dev(SZ, n, n_groups), du_dev(SZ, n, n_groups))
#else
    allocate (u(SZ, n, n_groups), du(SZ, n, n_groups))
#endif
  end subroutine allocate_fields

  subroutine setup_backend()
#ifdef CUDA
    blocks = dim3(n_groups, 1, 1)
    threads = dim3(SZ, 1, 1)

    ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
    ierr = cudaDeviceGetAttribute(memBusWidth, &
                                  cudaDevAttrGlobalMemoryBusWidth, 0)
#endif
  end subroutine setup_backend

  subroutine run_case(case_name, delta, bc_start, bc_end, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: delta, consumed_bw
    integer, intent(in) :: bc_start, bc_end

    integer :: iter
    real(dp) :: tstart, tend

    call initialise_input(delta)

#ifdef CUDA
    u_dev = u
#endif

    tdsops = tdsops_init(n, delta, &
                         operation='second-deriv', scheme='compact6', &
                         bc_start=bc_start, bc_end=bc_end)

    do iter = 1, n_warmup
      call run_kernel()
    end do
    call sync_backend()

    call start_timer(tstart)
    do iter = 1, n_iters
      call run_kernel()
    end do
    call sync_backend()
    call stop_timer(tend)

    call write_perf_metric(trim(backend_label)//'_thom_'//trim(case_name), &
                           tend - tstart, n_iters, ndof, consumed_bw)
#ifdef CUDA
    call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw, &
                            memClockRt, memBusWidth)
#else
    call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw)
#endif
  end subroutine run_case

  subroutine initialise_input(delta)
    real(dp), intent(in) :: delta

    do k = 1, n_groups
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = sin((j - 1)*delta)
        end do
      end do
    end do
  end subroutine initialise_input

  subroutine run_kernel()
#ifdef CUDA
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
#else
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
#endif
  end subroutine run_kernel

  subroutine sync_backend()
#ifdef CUDA
    integer :: ierr

    ierr = cudaDeviceSynchronize()
#endif
  end subroutine sync_backend

  subroutine start_timer(t)
    real(dp), intent(out) :: t

#ifdef CUDA
    call cpu_time(t)
#else
    t = omp_get_wtime()
#endif
  end subroutine start_timer

  subroutine stop_timer(t)
    real(dp), intent(out) :: t

#ifdef CUDA
    call cpu_time(t)
#else
    t = omp_get_wtime()
#endif
  end subroutine stop_timer

end program perf_thom
