module m_test_utils
  use m_common, only: dp, nbytes
  implicit none

contains

  subroutine check_norm(norm, tol, label, allpass)
    real(dp), intent(in) :: norm
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    logical, intent(inout) :: allpass

    if (norm > tol) then
      print '(a, a, a, es12.4, a, es12.4)', &
        'CHECK ', label, ' ... FAILED (norm=', norm, ' tol=', tol, ')'
      allpass = .false.
    else
      print '(a, a, a, es12.4)', &
        'CHECK ', label, ' ... PASSED (norm=', norm, ')'
    end if
  end subroutine check_norm

  subroutine report_perf(label, time, n_iters, ndof, consumed_bw)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: time
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw

    real(dp) :: achievedBW

    achievedBW = consumed_bw*n_iters*ndof*nbytes/time

    print '(a, a, a, f10.6, a, f10.3, a)', &
      'PERF_METRIC: ', label, ' time=', time, 's bw=', &
      achievedBW/real(2**30, dp), ' GiB/s'
  end subroutine report_perf

end module m_test_utils
