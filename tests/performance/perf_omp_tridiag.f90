program perf_omp_tridiag
  use mpi
  use omp_lib

  use m_common, only: dp, nbytes, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_omp_common, only: SZ
  use m_omp_sendrecv, only: sendrecv_fields
  use m_omp_exec_dist, only: exec_dist_tds_compact
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_test_utils, only: write_perf_minmax_metrics

  implicit none

  character(len=*), parameter :: backend_label = 'omp'
  real(dp), parameter :: periodic_bw = 5.0_dp

  real(dp), allocatable :: u(:, :, :), du(:, :, :)
  real(dp), allocatable :: u_recv_s(:, :, :), u_recv_e(:, :, :)
  real(dp), allocatable :: u_send_s(:, :, :), u_send_e(:, :, :)
  real(dp), allocatable :: send_s(:, :, :), send_e(:, :, :)
  real(dp), allocatable :: recv_s(:, :, :), recv_e(:, :, :)

  type(tdsops_t) :: tdsops

  integer :: n, n_groups, n_halo, n_iters, n_warmup, n_glob
  integer :: ndof, nrank, nproc, pprev, pnext
  integer :: ierr
  real(dp) :: dx_per

  call initialise_mpi()
  call configure_benchmark()
  call allocate_fields()
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
    n_groups = 64*64/SZ
    n_iters = 1000
    n_warmup = 10
    n_halo = 4
    ndof = n*n_groups*SZ
    dx_per = 2*pi/n_glob
  end subroutine configure_benchmark

  subroutine allocate_fields()
    allocate (u(SZ, n, n_groups), du(SZ, n, n_groups))
    allocate (u_send_s(SZ, n_halo, n_groups), u_send_e(SZ, n_halo, n_groups))
    allocate (u_recv_s(SZ, n_halo, n_groups), u_recv_e(SZ, n_halo, n_groups))
    allocate (send_s(SZ, 1, n_groups), send_e(SZ, 1, n_groups))
    allocate (recv_s(SZ, 1, n_groups), recv_e(SZ, 1, n_groups))
  end subroutine allocate_fields

  subroutine run_case(case_name, delta, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: delta
    real(dp), intent(in) :: consumed_bw

    real(dp) :: tstart, tend
    real(dp) :: achieved_bw, achieved_bw_min, achieved_bw_max

    call initialise_input(delta)

    tdsops = tdsops_init(n, delta, operation='second-deriv', &
                         scheme='compact6', &
                         bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    call run_iterations(n_warmup)
    call sync_backend()

    call start_timer(tstart)
    call run_iterations(n_iters)
    call sync_backend()
    call stop_timer(tend)

    call collect_perf_stats(tend - tstart, consumed_bw, achieved_bw, &
                            achieved_bw_min, achieved_bw_max)

    if (nrank == 0) then
      call write_perf_minmax_metrics(trim(backend_label)//'_tridiag_'//trim(case_name), &
                                     tend - tstart, &
                                     trim(backend_label)//'_tridiag_bw', &
                                     achieved_bw_min, achieved_bw_max)
    end if
  end subroutine run_case

  subroutine initialise_input(delta)
    real(dp), intent(in) :: delta

    integer :: j

    do j = 1, n
      u(:, j, :) = sin(((j - 1) + nrank*n)*delta)
    end do
  end subroutine initialise_input

  subroutine run_iterations(niters)
    integer, intent(in) :: niters

    integer :: iters, i, j, k

    do iters = 1, niters
      !$omp parallel do
      do k = 1, n_groups
        do j = 1, n_halo
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

  subroutine collect_perf_stats(time, consumed_bw, achieved_bw, achieved_bw_min, &
                                achieved_bw_max)
    real(dp), intent(in) :: time, consumed_bw
    real(dp), intent(out) :: achieved_bw, achieved_bw_min, achieved_bw_max

    achieved_bw = consumed_bw*n_iters*ndof*nbytes/time

    call MPI_Allreduce(achieved_bw, achieved_bw_max, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(achieved_bw, achieved_bw_min, 1, MPI_X3D2_DP, &
                       MPI_MIN, MPI_COMM_WORLD, ierr)
  end subroutine collect_perf_stats

  subroutine sync_backend()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine sync_backend

  subroutine start_timer(t)
    real(dp), intent(out) :: t

    t = omp_get_wtime()
  end subroutine start_timer

  subroutine stop_timer(t)
    real(dp), intent(out) :: t

    t = omp_get_wtime()
  end subroutine stop_timer

  subroutine finalise_benchmark()
    call MPI_Finalize(ierr)
  end subroutine finalise_benchmark

end program perf_omp_tridiag
