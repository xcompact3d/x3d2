program perf_cuda_penta
  !! Performance benchmark for the compact10_penta (Lele 10, pentadiagonal)
  !! first-derivative scheme.
  !!
  !! Uses the non-periodic single-GPU pentadiagonal Thomas solver implemented
  !! in der_penta_full. Scheme coefficients (Lele 1992, Table 1, 10th-order):
  !!   alpha=0.5, beta=0.05  (LHS pentadiag)
  !!   a=17/12, b=101/150, c=1/100  (RHS stencil)
  use cudafor
  use mpi

  use m_common, only: dp, nbytes, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_penta_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init
  use m_test_utils, only: write_perf_minmax_metrics, write_perf_minmax_summary, &
                          write_device_bw_metric

  implicit none

  character(len=*), parameter :: backend_label = 'cuda'
  ! Lele 10 penta 1st-deriv: reads 7 values, writes 1 → 8 words = bandwidth factor 8
  real(dp), parameter :: penta_bw = 8.0_dp

  real(dp), allocatable :: u(:, :, :)
  real(dp), device, allocatable :: u_dev(:, :, :), du_dev(:, :, :)
  real(dp), device, allocatable :: u_recv_s_dev(:, :, :), u_recv_e_dev(:, :, :)

  type(cuda_tdsops_t) :: tdsops

  integer :: n, n_block, n_halo, n_iters, n_warmup, n_glob, ndof, nrank, nproc
  integer :: pprev, pnext
  integer :: ierr
  integer :: ndevs, devnum
  integer :: memClockRt, memBusWidth
  type(dim3) :: blocks, threads
  real(dp) :: dx

  call initialise_mpi()
  call select_device()
  call configure_benchmark()
  call allocate_fields()
  call setup_backend()
  call run_case('penta_lele10_nonper', dx, penta_bw)
  call finalise()

contains

  subroutine initialise_mpi()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (nrank == 0) then
      print *, 'Performance benchmark for compact10_penta (Lele 10 penta 1st-deriv)'
      print *, 'Scheme: alpha=0.5 beta=0.05, non-periodic, single-GPU Thomas'
      print *, 'Ranks:', nproc
    end if

    pnext = modulo(nrank - nproc + 1, nproc)
    pprev = modulo(nrank - 1, nproc)
  end subroutine initialise_mpi

  subroutine configure_benchmark()
    n_glob = 1024
    n = n_glob/nproc
    n_block = 512*512/SZ
    n_iters = 100
    n_warmup = 10
    n_halo = 4
    ndof = n*n_block*SZ
    dx = 2*pi/n_glob
  end subroutine configure_benchmark

  subroutine allocate_fields()
    allocate (u(SZ, n, n_block))
    allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
    ! Halo arrays: for single-GPU non-periodic, halos are zero (Dirichlet ghost cells)
    allocate (u_recv_s_dev(SZ, n_halo, n_block))
    allocate (u_recv_e_dev(SZ, n_halo, n_block))
    u_recv_s_dev = 0._dp
    u_recv_e_dev = 0._dp
  end subroutine allocate_fields

  subroutine select_device()
    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(devnum)
  end subroutine select_device

  subroutine setup_backend()
    integer :: i, j, k

    do k = 1, n_block
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = sin((j - 1 + nrank*n)*dx)
        end do
      end do
    end do
    u_dev = u

    ! compact10_penta: Lele 10 pentadiagonal first derivative, periodic BCs
    ! used here only to set alpha/beta/coeffs; the kernel ignores periodicity
    ! at boundary rows (halo = 0 for non-periodic benchmark).
    tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                              scheme='compact10_penta', &
                              bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    blocks = dim3(n_block, 1, 1)
    threads = dim3(SZ, 1, 1)

    ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, devnum)
    ierr = cudaDeviceGetAttribute(memBusWidth, cudaDevAttrGlobalMemoryBusWidth, devnum)
  end subroutine setup_backend

  subroutine run_case(case_name, delta, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: delta
    real(dp), intent(in) :: consumed_bw

    integer :: iter
    real(dp) :: tstart, tend
    real(dp) :: achieved_bw, achieved_bw_min, achieved_bw_max

    call initialise_input(delta)

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

    call collect_perf_stats(tend - tstart, consumed_bw, achieved_bw, &
                            achieved_bw_min, achieved_bw_max)

    if (nrank == 0) then
      call write_perf_minmax_metrics(trim(backend_label)//'_penta_'//trim(case_name), &
                                     tend - tstart, &
                                     trim(backend_label)//'_penta_'//trim(case_name)//'_bw', &
                                     achieved_bw_min, achieved_bw_max)
      call write_perf_minmax_summary(achieved_bw_min, achieved_bw_max, &
                                     memClockRt, memBusWidth)
      call write_device_bw_metric(memClockRt, memBusWidth)
    end if
  end subroutine run_case

  subroutine initialise_input(delta)
    real(dp), intent(in) :: delta

    integer :: i, j, k

    do k = 1, n_block
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = sin((j - 1 + nrank*n)*delta)
        end do
      end do
    end do
    u_dev = u
  end subroutine initialise_input

  subroutine run_kernel()
    ! No halo exchange needed for single-GPU non-periodic solve
    call exec_dist_penta_compact(du_dev, u_dev, u_recv_s_dev, u_recv_e_dev, &
                                 tdsops, blocks, threads)
  end subroutine run_kernel

  subroutine collect_perf_stats(time, consumed_bw, achieved_bw, achieved_bw_min, &
                                achieved_bw_max)
    real(dp), intent(in) :: time
    real(dp), intent(in) :: consumed_bw
    real(dp), intent(out) :: achieved_bw
    real(dp), intent(out) :: achieved_bw_min
    real(dp), intent(out) :: achieved_bw_max

    achieved_bw = consumed_bw*n_iters*ndof*nbytes/time

    call MPI_Allreduce(achieved_bw, achieved_bw_max, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(achieved_bw, achieved_bw_min, 1, MPI_X3D2_DP, &
                       MPI_MIN, MPI_COMM_WORLD, ierr)
  end subroutine collect_perf_stats

  subroutine sync_backend()
    ierr = cudaDeviceSynchronize()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine sync_backend

  subroutine start_timer(t)
    real(dp), intent(out) :: t

    call cpu_time(t)
  end subroutine start_timer

  subroutine stop_timer(t)
    real(dp), intent(out) :: t

    call cpu_time(t)
  end subroutine stop_timer

  subroutine finalise()
    call MPI_Finalize(ierr)
  end subroutine finalise

end program perf_cuda_penta
