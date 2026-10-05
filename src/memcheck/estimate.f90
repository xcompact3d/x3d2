module m_memcheck_estimate
  !! Static (Tier 1) and Tier 2 per-GPU memory estimate of x3d2-memcheck.
  use m_common, only: dp, i8
  use m_memcheck_context, only: memcheck_ctx_t
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

  function ng_unsupported_reason(ctx, ng) result(reason)
    !! Whether GPU count ng is a decomposition this tool/solver actually
    !! supports - blank when it is. Shared by report()'s per-ng table,
    !! report_measured_table()'s per-ng table, and report()'s requested
    !! nproc_dir line, so the three can never silently diverge on what
    !! counts as unsupported. Checked in this order: ng itself must be a
    !! positive count; ng>1 needs BC_y periodic (multi_gpu_supported);
    !! gdims(3) must divide by ng (the physical Z decomposition); and, for
    !! whichever BC actually decomposes the spectral slab along dim2 (000/
    !! 010: cy: 100: cx - see spectral_slab_bytes), that dimension must
    !! also divide by ng, or the solver would truncate or error-stop.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: ng
    character(len=96) :: reason

    reason = ''
    if (ng < 1) then
      reason = 'nproc_dir(3) must be >= 1'
    else if (ng > 1 .and. .not. ctx%multi_gpu_supported) then
      reason = 'not supported: BC_y is non-periodic - &
               &src/poisson_fft.f90 error-stops at nproc>1'
    else if (mod(ctx%gdims(3), ng) /= 0) then
      write (reason, '(a,i0,a)') 'z=', ctx%gdims(3), ' not divisible'
    else if ((ctx%bc_is_000 .or. ctx%bc_is_010) .and. &
             mod(ctx%cdims(2), ng) /= 0) then
      write (reason, '(a,i0,a)') 'y cells=', ctx%cdims(2), ' not divisible &
        &by ng, the solver would truncate the spectral slab'
    else if (ctx%bc_is_100 .and. mod(ctx%cdims(1), ng) /= 0) then
      write (reason, '(a,i0,a)') 'x cells=', ctx%cdims(1), ' not divisible &
        &by ng, src/backend/cuda/poisson_fft.f90 error-stops'
    end if
  end function ng_unsupported_reason

end module m_memcheck_estimate
