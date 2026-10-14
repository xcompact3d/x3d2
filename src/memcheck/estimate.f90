module m_memcheck_estimate
  !! Static (Tier 1) and Tier 2 per-GPU memory estimate of x3d2-memcheck.
  use m_common, only: dp, i8
  implicit none

  !> <80% of card memory: FITS. 80-95%: BORDERLINE. >95%: DOES_NOT_FIT.
  real(dp), parameter :: FITS_FRACTION = 0.80_dp
  real(dp), parameter :: DOES_NOT_FIT_FRACTION = 0.95_dp

contains

  pure function to_gib(n_bytes) result(gib)
    !! Bytes to GiB, the one conversion every printed figure goes through
    !! so the digits cannot drift between call sites.
    integer(i8), intent(in) :: n_bytes
    real(dp) :: gib

    gib = real(n_bytes, dp)/1024._dp**3
  end function to_gib

  pure function classify(per_gpu_gib, total_gib) result(verdict)
    !! Shared three-way verdict classification, used by both the static
    !! estimate (estimate_for_ng) and the real --build measurement
    !! (report_measured_table) so the two paths can never silently diverge
    !! on what counts as FITS/BORDERLINE/DOES_NOT_FIT.
    real(dp), intent(in) :: per_gpu_gib, total_gib
    character(len=12) :: verdict

    if (per_gpu_gib < FITS_FRACTION*total_gib) then
      verdict = 'FITS'
    else if (per_gpu_gib <= DOES_NOT_FIT_FRACTION*total_gib) then
      verdict = 'BORDERLINE'
    else
      verdict = 'DOES_NOT_FIT'
    end if
  end function classify

end module m_memcheck_estimate
