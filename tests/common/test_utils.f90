module m_test_utils
  use m_common, only: dp, nbytes
  implicit none

  private
  public :: checkerr, write_perf_metric, write_perf_minmax_metrics, &
            write_perf_summary, write_perf_minmax_summary, &
            write_device_bw_metric

contains

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

  subroutine write_perf_summary(time, n_iters, ndof, consumed_bw, mem_clock_rt, &
                                mem_bus_width)
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

  pure real(dp) function compute_achieved_bw(time, n_iters, ndof, consumed_bw) &
    result(achievedBW)
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
