module m_memcheck_estimate
  !! Formula estimate and cuFFTMp probe per-GPU memory estimate of
  !! x3d2-memcheck.
  use m_common, only: dp, i8, nbytes
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memory_estimate, only: padded_cells, padded_halo_bytes, &
                               spectral_slab_bytes, mirror_buffer_bytes_100, &
                               spectral_extra_bytes_110, &
                               stretched_y_matrix_bytes, gpu_io_staging_bytes
  implicit none

  private
  public :: to_gib, classify, ng_unsupported_reason, &
            requested_unsupported_reason, estimate_for_ng, &
            fields_plus_halo_bytes_n, n_halo, basis_measured, &
            basis_extrapolated, basis_formula, basis_name

  !> <80% of card memory: FITS. 80-95%: BORDERLINE. >95%: DOES_NOT_FIT.
  real(dp), parameter :: FITS_FRACTION = 0.80_dp
  real(dp), parameter :: DOES_NOT_FIT_FRACTION = 0.95_dp
  integer, parameter :: n_halo = 4

  !> How the overhead of an estimate was obtained. measured: the ng=1 GPU
  !> FFT plan query. extrapolated: that query scaled to ng>1. formula: the
  !> formula estimate only, no plan query ran or it returned a failure.
  integer, parameter :: basis_measured = 1, basis_extrapolated = 2, &
                        basis_formula = 3
  character(len=12), parameter :: basis_name(3) = [character(len=12) :: &
                                                   'measured', &
                                                   'extrapolated', &
                                                   'formula-only']

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
    !! positive count; ng>1 needs BC_y periodic (multi_gpu_supported); for
    !! 000 and 100, ng>1 also needs the distributed FFT (cuFFTMp), checked
    !! once the FFT path is known (a probe or the real build has run);
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
    else if (ng > 1 .and. (ctx%bc_is_000 .or. ctx%bc_is_100) .and. &
             ctx%device%fft_path_known() .and. &
             (.not. ctx%device%uses_distributed_fft())) then
      ! Plain cuFFT cannot decompose 100 at all, and for 000 it runs a
      ! local per-rank transform with no exchange
      ! (src/backend/cuda/poisson_fft.f90:721-740), so the result would be
      ! wrong.
      reason = 'not supported: the FFT fell back to plain cuFFT, no &
               &multi-GPU transform'
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

  function requested_unsupported_reason(ctx) result(reason)
    !! Whether the nproc_dir the input file asks for is supported - blank
    !! when it is. ng_unsupported_reason judges the Z count nproc_dir(3);
    !! on top of that the CUDA backend only supports nproc_dir = [1,1,ng].
    !! Shared by report()'s requested line and run_build_mode, so the initial
    !! report and the --build verdict cannot disagree on it.
    type(memcheck_ctx_t), intent(in) :: ctx
    character(len=96) :: reason

    reason = ng_unsupported_reason(ctx, ctx%domain_cfg%nproc_dir(3))
    if (ctx%domain_cfg%nproc_dir(1) /= 1 .or. &
        ctx%domain_cfg%nproc_dir(2) /= 1) &
      reason = 'only nproc_dir = [1,1,ng] is supported by the CUDA backend'
  end function requested_unsupported_reason

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
    !! Exact, ng dependent field+halo term (formula estimate, no GPU/plan
    !! needed), using the static peak_fields_lookup value.
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
      nbytes8 = nbytes8 + stretched_y_matrix_bytes( &
                ctx%bc_is_010, ctx%domain_cfg%stretching(2), &
                ctx%solver_cfg%lowmem_fft, ctx%cdims, ng)
    end if
  end function spectral_plus_mirror_bytes

  subroutine estimate_for_ng(ctx, ng, per_gpu_gib, basis, verdict, &
                             workspace_gib, overhead_gib, io_gib, used_gib, &
                             npeak)
    !! Formula floor first; only calls into the cuFFTMp probe (a real,
    !! throwaway GPU FFT plan) when the floor alone cannot already answer
    !! DOES_NOT_FIT, and
    !! never under --static. workspace_gib/overhead_gib (optional) are the
    !! same two-term breakdown the --build measured table reports:
    !! workspace = exact fields+halo+spectral+mirror (formula estimate, no GPU
    !! needed); overhead = everything from the cuFFTMp probe, the
    !! measured values (worksize/xtdesc/heap/context), and the
    !! GPU-aware IO staging term (also broken out on its own via the
    !! optional io_gib, in GiB, for callers that need to track it apart
    !! from context/FFT overhead - e.g. run_build_mode()'s CHECK line).
    !! used_gib and npeak (optional) are the ng=1 measurement of the real
    !! build and its measured peak_fields. The estimate is then that
    !! measurement scaled to ng GPUs by the same terms (spectral slab, mirror
    !! buffers, the NVSHMEM heap over ng when the real solver's FFT is
    !! distributed, IO staging), and workspace_gib holds the fields+halo
    !! alone.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: ng
    real(dp), intent(out) :: per_gpu_gib
    integer, intent(out) :: basis
    character(len=12), intent(out) :: verdict
    real(dp), intent(out), optional :: workspace_gib, overhead_gib, io_gib
    real(dp), intent(in), optional :: used_gib
    integer, intent(in), optional :: npeak

    integer(i8) :: base_bytes, spec_bytes, floor_bytes, overhead_bytes, &
                   io_bytes, total_b, free_b, probe_bytes
    real(dp) :: floor_gib, io_gib_local, residual_gib
    logical :: probe_skipped
    integer :: peak

    peak = ctx%peak_fields
    if (present(npeak)) peak = npeak
    base_bytes = fields_plus_halo_bytes_n(ctx, ng, peak)
    spec_bytes = spectral_plus_mirror_bytes(ctx, ng)
    io_bytes = gpu_io_staging_bytes([ctx%gdims(1), ctx%gdims(2), &
                                     ctx%gdims(3)/ng], &
                                    ctx%checkpoint_cfg, &
                                    ctx%gpu_io_device_write)
    io_gib_local = to_gib(io_bytes)
    if (present(io_gib)) io_gib = io_gib_local
    if (present(workspace_gib)) &
      workspace_gib = to_gib(base_bytes + spec_bytes)
    if (present(used_gib)) then
      per_gpu_gib = used_gib + to_gib(base_bytes + spec_bytes + io_bytes - &
                                      fields_plus_halo_bytes_n(ctx, 1, peak) &
                                      - spectral_plus_mirror_bytes(ctx, 1))
      if (ctx%device%fft_path_known() .and. &
          ctx%device%uses_distributed_fft()) &
        per_gpu_gib = per_gpu_gib + &
                      to_gib(ctx%device%overhead_floor_bytes( &
                             ng, ctx%bc_is_110) - &
                             ctx%device%overhead_floor_bytes( &
                             1, ctx%bc_is_110))
      basis = merge(basis_measured, basis_extrapolated, ng == 1)
      if (present(workspace_gib)) workspace_gib = to_gib(base_bytes)
      if (present(overhead_gib)) &
        overhead_gib = per_gpu_gib - to_gib(base_bytes)
      verdict = classify(per_gpu_gib, ctx%card_gib)
      return
    end if
    ! Floor: assume cuFFTMp is used (the realistic case for every BC this
    ! solver ever attempts it for - see the 110 guard of
    ! fft_workspace_bytes_query) since that gives the larger, safer lower
    ! bound.
    floor_bytes = base_bytes + spec_bytes + io_bytes + &
                  ctx%device%overhead_floor_bytes(ng, ctx%bc_is_110)
    floor_gib = to_gib(floor_bytes)

    ! The cuFFTMp probe builds the FFT plans on the whole global grid, one
    ! GPU, whatever ng is. Skip it when free memory cannot hold that:
    ! seven spectral copies plus the fixed overhead. cuFFTMp aborts inside
    ! NVSHMEM on a plan that does not fit and returns no status. On a 40 GB
    ! card the heap took 2.5 to 5 copies, depending on the grid size, and
    ! the xtdesc buffer one. The 984^3 probe aborted with a guess of five
    ! copies, the 1056^3 probe with four and the 1280^3 probe with two.
    probe_skipped = .false.
    if (trim(ctx%run_mode) /= 'STATIC' .and. &
        floor_gib <= DOES_NOT_FIT_FRACTION*ctx%card_gib) then
      if (.not. ctx%device%query_residual_gib(residual_gib)) then
        call ctx%device%mem_info(total_b, free_b)
        probe_bytes = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, &
                                          ctx%cdims, 1)
        probe_skipped = free_b < 7_i8*probe_bytes + &
                        ctx%device%overhead_floor_bytes(1, ctx%bc_is_110)
      end if
    end if

    if (trim(ctx%run_mode) == 'STATIC' .or. probe_skipped .or. &
        floor_gib > DOES_NOT_FIT_FRACTION*ctx%card_gib) then
      per_gpu_gib = floor_gib
      basis = basis_formula
      if (present(overhead_gib)) &
        overhead_gib = to_gib(ctx%device%overhead_floor_bytes( &
                              ng, ctx%bc_is_110)) + io_gib_local
      if (trim(ctx%run_mode) == 'STATIC' .or. probe_skipped) then
        ! --static: formula estimate only, no cuFFTMp probe. floor_gib
        ! IS the estimate here - not just a DOES_NOT_FIT early-return floor
        ! - so classify it with the full three-way verdict. A skipped probe
        ! is the same case.
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
    ! A probe that returned a failure status leaves the formula. A probe at
    ! ng>1 is the single-GPU probe scaled by 1/ng, not a measurement.
    if (ctx%device%query_residual_gib(residual_gib)) then
      basis = merge(basis_measured, basis_extrapolated, ng == 1)
    else
      basis = basis_formula
    end if
    if (present(overhead_gib)) &
      overhead_gib = to_gib(overhead_bytes) + io_gib_local

    verdict = classify(per_gpu_gib, ctx%card_gib)
  end subroutine estimate_for_ng

end module m_memcheck_estimate
