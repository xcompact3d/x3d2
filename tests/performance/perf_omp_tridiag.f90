program perf_omp_tridiag
  use mpi
  use omp_lib

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_omp_common, only: SZ
  use m_omp_sendrecv, only: sendrecv_fields
  use m_omp_exec_dist, only: exec_dist_tds_compact
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_test_utils, only: report_perf

  implicit none

  real(dp), allocatable, dimension(:, :, :) :: u, du
  real(dp), allocatable, dimension(:, :, :) :: u_recv_s, u_recv_e, &
                                               u_send_s, u_send_e
  real(dp), allocatable, dimension(:, :, :) :: send_s, send_e, &
                                               recv_s, recv_e

  type(tdsops_t) :: tdsops

  integer :: n, n_groups, j, n_halo, n_iters, n_warmup
  integer :: n_glob
  integer :: nrank, nproc, pprev, pnext
  integer :: ierr

  real(dp) :: dx_per, tstart, tend
  real(dp) :: achievedBW, achievedBWmin, achievedBWmax

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Performance benchmark with', nproc, 'ranks'

  pnext = modulo(nrank - nproc + 1, nproc)
  pprev = modulo(nrank - 1, nproc)

  ! Performance configuration
  n_glob = 1024
  n = n_glob/nproc
  n_groups = 64*64/SZ
  n_iters = 1000
  n_warmup = 10

  allocate (u(SZ, n, n_groups), du(SZ, n, n_groups))

  dx_per = 2*pi/n_glob

  do j = 1, n
    u(:, j, :) = sin(((j - 1) + nrank*n)*dx_per)
  end do

  n_halo = 4
  allocate (u_send_s(SZ, n_halo, n_groups))
  allocate (u_send_e(SZ, n_halo, n_groups))
  allocate (u_recv_s(SZ, n_halo, n_groups))
  allocate (u_recv_e(SZ, n_halo, n_groups))
  allocate (send_s(SZ, 1, n_groups), send_e(SZ, 1, n_groups))
  allocate (recv_s(SZ, 1, n_groups), recv_e(SZ, 1, n_groups))

  tdsops = tdsops_init(n, dx_per, operation='second-deriv', &
                       scheme='compact6', &
                       bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  ! Warmup
  call run_iterations(n_warmup)

  ! Timed run
  tstart = omp_get_wtime()
  call run_iterations(n_iters)
  tend = omp_get_wtime()

  ! MPI-reduced BW metrics
  achievedBW = 5._dp*n_iters*n*n_groups*SZ*dp/(tend - tstart)
  call MPI_Allreduce(achievedBW, achievedBWmax, 1, MPI_X3D2_DP, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(achievedBW, achievedBWmin, 1, MPI_X3D2_DP, &
                     MPI_MIN, MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    print '(a, f10.6, a)', 'PERF_METRIC: omp_tridiag_periodic time=', &
      tend - tstart, 's'
    print '(a, f10.3, a)', 'PERF_METRIC: omp_tridiag_bw_min bw=', &
      achievedBWmin/real(2**30, dp), ' GiB/s'
    print '(a, f10.3, a)', 'PERF_METRIC: omp_tridiag_bw_max bw=', &
      achievedBWmax/real(2**30, dp), ' GiB/s'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine run_iterations(niters)
    integer, intent(in) :: niters
    integer :: iters, i, j, k

    do iters = 1, niters
      !$omp parallel do
      do k = 1, n_groups
        do j = 1, 4
          !$omp simd
          do i = 1, SZ
            u_send_s(i, j, k) = u(i, j, k)
            u_send_e(i, j, k) = u(i, n - n_halo + j, k)
          end do
          !$omp end simd
        end do
      end do
      !$omp end parallel do

      call sendrecv_fields(u_recv_s, u_recv_e, u_send_s, u_send_e, &
                           SZ*n_halo*n_groups, nproc, pprev, pnext)

      call exec_dist_tds_compact(du, u, u_recv_s, u_recv_e, &
                                 send_s, send_e, recv_s, recv_e, &
                                 tdsops, nproc, pprev, pnext, n_groups)
    end do
  end subroutine run_iterations

end program perf_omp_tridiag
