program perf_cuda_tridiag
  use cudafor
  use mpi

  use m_common, only: dp, nbytes, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_tds_compact
  use m_cuda_sendrecv, only: sendrecv_fields
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init
  use m_test_utils, only: write_perf_minmax_metrics, write_perf_minmax_summary, &
                          write_device_bw_metric

  implicit none

  character(len=*), parameter :: backend_label = 'cuda'
  real(dp), parameter :: periodic_bw = 6.0_dp

  real(dp), allocatable :: u(:, :, :)
  real(dp), device, allocatable :: u_dev(:, :, :), du_dev(:, :, :)
  real(dp), device, allocatable :: u_recv_s_dev(:, :, :), u_recv_e_dev(:, :, :)
  real(dp), device, allocatable :: u_send_s_dev(:, :, :), u_send_e_dev(:, :, :)
  real(dp), device, allocatable :: du_send_s_dev(:, :, :), du_send_e_dev(:, :, :)
  real(dp), device, allocatable :: du_recv_s_dev(:, :, :), du_recv_e_dev(:, :, :)

  type(cuda_tdsops_t) :: tdsops

  integer :: n, n_block, n_halo, n_iters, n_warmup, n_glob, ndof, nrank, nproc
  integer :: pprev, pnext
  integer :: ierr
  integer :: ndevs, devnum
  integer :: memClockRt, memBusWidth
  type(dim3) :: blocks, threads
  real(dp) :: dx_per

  call initialise_mpi()
  call configure_benchmark()
  call allocate_fields()
  call setup_backend()
  call run_case('periodic', dx_per, periodic_bw)
  call finalise_benchmark()

contains

  subroutine initialise_mpi()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (nrank == 0) print *, 'Performance benchmark with', nproc, 'ranks'

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
    dx_per = 2*pi/n_glob
  end subroutine configure_benchmark

  subroutine allocate_fields()
    allocate (u(SZ, n, n_block))
    allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))

    allocate (u_send_s_dev(SZ, n_halo, n_block), u_send_e_dev(SZ, n_halo, n_block))
    allocate (u_recv_s_dev(SZ, n_halo, n_block), u_recv_e_dev(SZ, n_halo, n_block))
    allocate (du_send_s_dev(SZ, 1, n_block), du_send_e_dev(SZ, 1, n_block))
    allocate (du_recv_s_dev(SZ, 1, n_block), du_recv_e_dev(SZ, 1, n_block))
  end subroutine allocate_fields

  subroutine setup_backend()
    integer :: i, j, k

    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(devnum)

    do k = 1, n_block
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
        end do
      end do
    end do
    u_dev = u

    tdsops = cuda_tdsops_init(n, dx_per, operation='second-deriv', &
                              scheme='compact6', &
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
      call write_perf_minmax_metrics(trim(backend_label)//'_tridiag_'//trim(case_name), &
                                     tend - tstart, &
                                     trim(backend_label)//'_tridiag_'//trim(case_name)//'_bw', &
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
    u_send_s_dev(:, :, :) = u_dev(:, 1:n_halo, :)
    u_send_e_dev(:, :, :) = u_dev(:, n - n_halo + 1:n, :)

    call sendrecv_fields(u_recv_s_dev, u_recv_e_dev, &
                         u_send_s_dev, u_send_e_dev, &
                         SZ*n_halo*n_block, nproc, pprev, pnext)

    call exec_dist_tds_compact(du_dev, u_dev, u_recv_s_dev, u_recv_e_dev, &
                               du_send_s_dev, du_send_e_dev, &
                               du_recv_s_dev, du_recv_e_dev, &
                               tdsops, nproc, pprev, pnext, blocks, threads)
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

  subroutine finalise_benchmark()
    call MPI_Finalize(ierr)
  end subroutine finalise_benchmark

end program perf_cuda_tridiag
