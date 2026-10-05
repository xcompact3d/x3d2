module m_memcheck_estimate
  !! Static (Tier 1) and Tier 2 per-GPU memory estimate of x3d2-memcheck.
  use m_common, only: dp, i8, nbytes
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memory_estimate, only: padded_cells, padded_halo_bytes, &
                               spectral_slab_bytes, mirror_buffer_bytes_100, &
                               spectral_extra_bytes_110, &
                               stretched_y_matrix_bytes, gpu_io_staging_bytes
  implicit none

  private
  public :: to_gib, classify, ng_unsupported_reason, estimate_for_ng, &
            fields_plus_halo_bytes_n, n_halo

  !> <80% of card memory: FITS. 80-95%: BORDERLINE. >95%: DOES_NOT_FIT.
  real(dp), parameter :: FITS_FRACTION = 0.80_dp
  real(dp), parameter :: DOES_NOT_FIT_FRACTION = 0.95_dp
  integer, parameter :: n_halo = 4

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

  function fields_plus_halo_bytes_n(ctx, ng, npeak) result(nbytes8)
    !! Same formula as fields_plus_halo_bytes, but for an explicit npeak
    !! rather than the module's static peak_fields - used by the --build
    !! measured table, which has its own real measured field count.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: ng, npeak
    integer(i8) :: nbytes8
    integer :: local_dims(3)

    local_dims = [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)/ng]
    nbytes8 = int(npeak, i8)*padded_cells(local_dims, ctx%device%sz()) &
              *int(nbytes, i8) + &
              padded_halo_bytes(local_dims, ctx%device%sz(), n_halo)
  end function fields_plus_halo_bytes_n

  function fields_plus_halo_bytes(ctx, ng) result(nbytes8)
    !! Exact, ng dependent field+halo term (Tier 1, no GPU/plan needed),
    !! using the static peak_fields_lookup value.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: ng
    integer(i8) :: nbytes8

    nbytes8 = fields_plus_halo_bytes_n(ctx, ng, ctx%peak_fields)
  end function fields_plus_halo_bytes

  function spectral_plus_mirror_bytes(ctx, ng) result(nbytes8)
    !! waves_dev (all cases, src/backend/cuda/poisson_fft.f90:322-324) is
    !! exactly spectral_slab_bytes. Two BC-specific extra terms:
    !!   100, ng>1: mirror/exchange buffers (:339-347), zero at ng<=1.
    !!   110 (never ng>1 - see multi_gpu_supported): a SECOND complex array
    !!     of the same spectral shape, c_dev (:416-418), plus a real,
    !!     globally-shaped transposed workspace r_dev_110(nz_glob,nx_glob,
    !!     ny_glob) (:414) - i.e. the same total element count as cdims,
    !!     just permuted. Missing this term underestimated the 110 case by
    !!     ~14% (0.98 vs 1.134 GiB measured) during this feature's design.
    !!   010, non-uniform y-stretching (never ng>1 - 010 error-stops at
    !!     nproc>1): the a_*_dev Poisson coefficient matrices, see
    !!     stretched_y_matrix_bytes.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: ng
    integer(i8) :: nbytes8

    nbytes8 = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, ng)
    if (ctx%bc_is_100) then
      nbytes8 = nbytes8 + mirror_buffer_bytes_100(ctx%cdims, ng)
    else if (ctx%bc_is_110) then
      nbytes8 = nbytes8 + spectral_extra_bytes_110(ctx%cdims, ng)
    else if (ctx%bc_is_010) then
      nbytes8 = nbytes8 + stretched_y_matrix_bytes(ctx%bc_is_010, &
                                                ctx%domain_cfg%stretching(2), &
                                      ctx%solver_cfg%lowmem_fft, ctx%cdims, ng)
    end if
  end function spectral_plus_mirror_bytes

  subroutine estimate_for_ng(ctx, ng, per_gpu_gib, exact, verdict, &
                             workspace_gib, overhead_gib, io_gib)
    !! Tier 1 floor first; only calls into Tier 2 (a real, throwaway GPU
    !! FFT plan) when Tier 1 alone cannot already answer DOES_NOT_FIT, and
    !! never under --static. workspace_gib/overhead_gib (optional) are the
    !! same two-term breakdown the --build measured table reports:
    !! workspace = exact fields+halo+spectral+mirror (Tier 1, no GPU
    !! needed); overhead = everything from the Tier 2 GPU query, the
    !! hardware constants (worksize/xtdesc/heap/context), and the
    !! GPU-aware IO staging term (also broken out on its own via the
    !! optional io_gib, in GiB, for callers that need to track it apart
    !! from context/FFT overhead - e.g. run_tier3()'s CHECK line).
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: ng
    real(dp), intent(out) :: per_gpu_gib
    logical, intent(out) :: exact
    character(len=12), intent(out) :: verdict
    real(dp), intent(out), optional :: workspace_gib, overhead_gib, io_gib

    integer(i8) :: base_bytes, spec_bytes, floor_bytes, overhead_bytes, &
                   io_bytes
    real(dp) :: floor_gib, io_gib_local

    base_bytes = fields_plus_halo_bytes(ctx, ng)
    spec_bytes = spectral_plus_mirror_bytes(ctx, ng)
    io_bytes = gpu_io_staging_bytes([ctx%gdims(1), ctx%gdims(2), &
                                     ctx%gdims(3)/ng], &
                                    ctx%checkpoint_cfg, &
                                    ctx%gpu_io_device_write)
    io_gib_local = to_gib(io_bytes)
    if (present(io_gib)) io_gib = io_gib_local
    if (present(workspace_gib)) &
      workspace_gib = to_gib(base_bytes + spec_bytes)
    ! Floor: assume cuFFTMp is used (the realistic case for every BC this
    ! solver ever attempts it for - see m_cuda_memory_estimate's 110 guard)
    ! since that gives the larger, safer lower bound.
    floor_bytes = base_bytes + spec_bytes + io_bytes + &
                  ctx%device%overhead_floor_bytes(ng, ctx%bc_is_110)
    floor_gib = to_gib(floor_bytes)

    if (trim(ctx%run_mode) == 'STATIC' .or. &
        floor_gib > DOES_NOT_FIT_FRACTION*ctx%card_gib) then
      per_gpu_gib = floor_gib
      exact = .false.
      if (present(overhead_gib)) &
        overhead_gib = to_gib(ctx%device%overhead_floor_bytes( &
                              ng, ctx%bc_is_110)) + io_gib_local
      if (trim(ctx%run_mode) == 'STATIC') then
        ! --static: Tier 1 only, never queries the GPU FFT plan. floor_gib
        ! IS the estimate here - not just a DOES_NOT_FIT early-return floor
        ! - so classify it with the full three-way verdict.
        verdict = classify(per_gpu_gib, ctx%card_gib)
      else
        verdict = 'DOES_NOT_FIT'
      end if
      return
    end if

    overhead_bytes = ctx%device%overhead_query_bytes( &
                     ng, ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, &
                     ctx%irank == 0)
    per_gpu_gib = to_gib(base_bytes + spec_bytes + overhead_bytes + io_bytes)
    exact = .true.
    if (present(overhead_gib)) &
      overhead_gib = to_gib(overhead_bytes) + io_gib_local

    verdict = classify(per_gpu_gib, ctx%card_gib)
  end subroutine estimate_for_ng

end module m_memcheck_estimate
