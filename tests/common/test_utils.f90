module m_test_utils
  use mpi
  use iso_fortran_env, only: stderr => error_unit
  use m_common, only: dp, nbytes
  implicit none

  private
  public :: initialise_mpi, finalise_test, global_all, global_sum, checkerr, &
            write_perf_metric, write_perf_minmax_metrics, &
            write_perf_summary, write_perf_minmax_summary, &
            write_device_bw_metric

contains

  subroutine initialise_mpi(nrank, nproc, pprev, pnext)
    !! Initialise MPI and return the rank and communicator size. Optionally
    !! return the neighbouring ranks for a 1D periodic decomposition.
    integer, intent(out) :: nrank, nproc
    integer, optional, intent(out) :: pprev, pnext

    integer :: ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (present(pprev)) pprev = modulo(nrank - 1, nproc)
    if (present(pnext)) pnext = modulo(nrank - nproc + 1, nproc)
  end subroutine initialise_mpi

  subroutine global_all(value)
    !! Replace a local logical value with the logical AND across all ranks.
    logical, intent(inout) :: value

    integer :: ierr

    call MPI_Allreduce(MPI_IN_PLACE, value, 1, MPI_LOGICAL, MPI_LAND, &
                       MPI_COMM_WORLD, ierr)
  end subroutine global_all

  subroutine global_sum(value)
    !! Replace a local integer value with the sum across all ranks.
    integer, intent(inout) :: value

    integer :: ierr

    call MPI_Allreduce(MPI_IN_PLACE, value, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
  end subroutine global_sum

  subroutine finalise_test(allpass, nrank, finalize_mpi)
    !! Report the aggregate test result and finalise MPI, aborting via
    !! `error stop` (non-zero exit status) if any check failed. When nrank is
    !! given, the success message is only printed by rank 0. Set finalize_mpi
    !! to .false. for serial tests that never called MPI_Init.
    logical, intent(in) :: allpass
    integer, optional, intent(in) :: nrank
    logical, optional, intent(in) :: finalize_mpi

    integer :: ierr
    logical :: is_root, do_finalize

    is_root = .true.
    if (present(nrank)) is_root = (nrank == 0)

    do_finalize = .true.
    if (present(finalize_mpi)) do_finalize = finalize_mpi

    if (allpass) then
      if (is_root) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if

    if (do_finalize) call MPI_Finalize(ierr)
  end subroutine finalise_test

  subroutine check_norm(norm, tol, label, allpass)
    real(dp), intent(in) :: norm
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    logical, intent(inout) :: allpass

    if (norm > tol) then
      print '(a, a, a, es12.4, a, es12.4, a)', &
        'CHECK ', label, ' ... FAILED (norm=', norm, ' tol=', tol, ')'
      allpass = .false.
    else
      print '(a, a, a, es12.4, a)', &
        'CHECK ', label, ' ... PASSED (norm=', norm, ')'
    end if
  end subroutine check_norm

  subroutine checkerr(u, du, tol, label, allpass)
    real(dp), intent(in) :: u(:, :, :)
    real(dp), intent(in) :: du(:, :, :)
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    logical, intent(inout) :: allpass

    real(dp) :: norm_residual

    norm_residual = sum((u + du)**2)/real(size(u), dp)
    norm_residual = sqrt(norm_residual)

    print *, "Check error:"
    print *, "min:", minval(u + du), "max: ", maxval(u + du)
    print *, "error norm", norm_residual

    call check_norm(norm_residual, tol, label, allpass)
  end subroutine checkerr

  subroutine write_perf_metric(label, time, n_iters, ndof, consumed_bw)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: time
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw

    real(dp) :: achievedBW

    achievedBW = compute_achieved_bw(time, n_iters, ndof, consumed_bw)

    print '(a, a, a, f10.6, a, f10.3, a)', &
      'PERF_METRIC: ', label, ' time=', time, 's bw=', &
      achievedBW/real(2**30, dp), ' GiB/s'
  end subroutine write_perf_metric

  subroutine write_perf_minmax_metrics(time_label, time, bw_label, bw_min, bw_max)
    character(len=*), intent(in) :: time_label
    real(dp), intent(in) :: time
    character(len=*), intent(in) :: bw_label
    real(dp), intent(in) :: bw_min
    real(dp), intent(in) :: bw_max

    print '(a, a, a, f10.6, a)', &
      'PERF_METRIC: ', time_label, ' time=', time, 's'
    print '(a, a, a, f10.3, a)', &
      'PERF_METRIC: ', trim(bw_label)//'_min', ' bw=', &
      bw_min/real(2**30, dp), ' GiB/s'
    print '(a, a, a, f10.3, a)', &
      'PERF_METRIC: ', trim(bw_label)//'_max', ' bw=', &
      bw_max/real(2**30, dp), ' GiB/s'
  end subroutine write_perf_minmax_metrics

  subroutine write_perf_summary(time, n_iters, ndof, consumed_bw, &
                                mem_clock_rt, mem_bus_width)
    real(dp), intent(in) :: time
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw
    integer, optional, intent(in) :: mem_clock_rt
    integer, optional, intent(in) :: mem_bus_width

    real(dp) :: achievedBW, deviceBW, utilisation
    integer :: resolved_mem_clock_rt, resolved_mem_bus_width

    if (present(mem_clock_rt)) then
      resolved_mem_clock_rt = mem_clock_rt
    else
      resolved_mem_clock_rt = 3200000
    end if

    if (present(mem_bus_width)) then
      resolved_mem_bus_width = mem_bus_width
    else
      resolved_mem_bus_width = 64
    end if

    achievedBW = compute_achieved_bw(time, n_iters, ndof, consumed_bw)
    deviceBW = compute_device_bw(resolved_mem_clock_rt, resolved_mem_bus_width)
    utilisation = achievedBW/deviceBW*100

    print *, "Check performance:"
    print '(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'
    print '(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
    print '(a, f5.2)', 'Effective BW util: %', utilisation
  end subroutine write_perf_summary

  subroutine write_perf_minmax_summary(bw_min, bw_max, mem_clock_rt, mem_bus_width)
    real(dp), intent(in) :: bw_min
    real(dp), intent(in) :: bw_max
    integer, intent(in) :: mem_clock_rt
    integer, intent(in) :: mem_bus_width

    real(dp) :: deviceBW, utilisation_min, utilisation_max

    deviceBW = compute_device_bw(mem_clock_rt, mem_bus_width)
    utilisation_min = bw_min/deviceBW*100
    utilisation_max = bw_max/deviceBW*100

    print *, "Check performance:"
    print '(a, f8.3, a)', 'Achieved BW min: ', bw_min/2**30, ' GiB/s'
    print '(a, f8.3, a)', 'Achieved BW max: ', bw_max/2**30, ' GiB/s'
    print '(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
    print '(a, f5.2)', 'Effective BW util min: %', utilisation_min
    print '(a, f5.2)', 'Effective BW util max: %', utilisation_max
  end subroutine write_perf_minmax_summary

  subroutine write_device_bw_metric(mem_clock_rt, mem_bus_width)
    integer, intent(in) :: mem_clock_rt
    integer, intent(in) :: mem_bus_width

    real(dp) :: deviceBW

    deviceBW = compute_device_bw(mem_clock_rt, mem_bus_width)

    print '(a, f10.3, a)', 'PERF_METRIC: device_bw ref=', &
      deviceBW/real(2**30, dp), ' GiB/s'
  end subroutine write_device_bw_metric

  pure real(dp) function compute_achieved_bw(time, n_iters, ndof, &
                                             consumed_bw) result(achievedBW)
    real(dp), intent(in) :: time
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw

    achievedBW = consumed_bw*n_iters*ndof*nbytes/time
  end function compute_achieved_bw

  pure real(dp) function compute_device_bw(mem_clock_rt, mem_bus_width) result(deviceBW)
    integer, intent(in) :: mem_clock_rt
    integer, intent(in) :: mem_bus_width

    deviceBW = 2.0_dp*mem_bus_width/nbytes*mem_clock_rt*(10**3)
  end function compute_device_bw

end module m_test_utils
