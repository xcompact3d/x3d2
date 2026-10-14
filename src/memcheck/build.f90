module m_memcheck_build
  !! Tier 3 of x3d2-memcheck: the real case build, its measurement and the
  !! measured report.
  use m_common, only: dp, i8
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memcheck_report, only: n_gpu_list, print_rule, print_table_header, &
                               print_table_row
  use m_memcheck_estimate, only: to_gib, classify, ng_unsupported_reason, &
                                 fields_plus_halo_bytes_n
  use m_cuda_memory_estimate, only: context_floor_bytes
  use m_memory_estimate, only: spectral_slab_bytes, mirror_buffer_bytes_100, &
                               peak_fields_lookup, gpu_io_staging_bytes
  implicit none

contains

  subroutine check_peak_fields_table(ctx, measured)
    !! Drift detector for m_memory_estimate's peak_fields_lookup: compares
    !! its static, no build estimate against the real measured peak_fields
    !! (allocator%next_id) for the config this run just built, and warns
    !! (does not fail the tool) on any mismatch. This is what lets the
    !! static table grow/self-correct as the solver evolves, instead of
    !! silently going stale the next time someone adds a new persistent
    !! allocation to a code path the table already covers.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: measured
    integer :: predicted

    if (ctx%irank /= 0) return

    predicted = peak_fields_lookup(trim(ctx%domain_cfg%flow_case_name), &
                                   ctx%solver_cfg%n_species, &
                                   trim(ctx%les_cfg%model) /= 'none', &
                                   ctx%solver_cfg%ibm_on, &
                                   ctx%output_vorticity, &
                                   ctx%output_qcriterion, &
                                   ctx%solver_cfg%lowmem_transeq)
    if (predicted /= measured) then
      print '(a,i0,a,i0,a)', &
        'WARNING: peak_fields_lookup predicted ', predicted, &
        ' but the real build measured ', measured, &
        ' - m_memory_estimate''s table is stale for this config and &
        &should be recalibrated (see src/memory_estimate.f90).'
    end if
  end subroutine check_peak_fields_table

  subroutine report_measured_table(ctx, measured_peak_fields, used_gib, &
                                   workspace_gib_measured, requested_verdict)
    !! Print the real, per-GPU-measured table. "workspace" at ng=1 is the
    !! exact fields+halo term with the MEASURED peak_fields; "overhead" is
    !! used_gib - workspace_gib_measured (context + FFT/spectral, whatever
    !! it actually was). For ng>1, apply the same exact ng-scaling
    !! corrections estimate_for_ng already encodes on top of that ng=1
    !! overhead: the spectral array shrinking (relative to the ng=1 slab
    !! the measurement already includes), the 100 case's multi rank only
    !! mirror buffers, and (cuFFTMp cases) the NVSHMEM heap's 1/ng shrink
    !! (confirmed by a real mpirun -n 2 xcompact run on
    !! examples/TGV/input.x3d, 2026-09-10). multi_gpu_supported here uses
    !! the REAL use_cufftmp from the build just performed, more
    !! authoritative than the static table's Tier 2 plan-query signal.
    !! The GPU-aware IO staging term (gpu_io_staging_bytes) is added the
    !! same way, since the real build performs no snapshot/checkpoint write
    !! and so never measures it directly.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: measured_peak_fields
    real(dp), intent(in) :: used_gib, workspace_gib_measured
    character(len=12), intent(out) :: requested_verdict

    integer :: k, ng, local_dims(3)
    real(dp) :: overhead, per_gpu, overhead_term, mirror_gib, w_local
    integer(i8) :: spec_bytes_1
    logical :: multi_gpu_supported_measured
    character(len=12) :: verdict
    character(len=96) :: reason

    spec_bytes_1 = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, &
                                       ctx%cdims, 1)
    overhead = used_gib - workspace_gib_measured
    requested_verdict = 'DOES_NOT_FIT'

    multi_gpu_supported_measured = .true.
    if (ctx%bc_is_010 .or. ctx%bc_is_110) then
      ! BC_y non-periodic: no nproc>1 support for this BC
      ! (src/poisson_fft.f90:178-180,196-198).
      multi_gpu_supported_measured = .false.
    else if ((ctx%bc_is_100 .or. ctx%bc_is_000) .and. &
             ctx%use_cufftmp_known .and. (.not. ctx%use_cufftmp)) then
      ! Covers the 100 case (needs cuFFTMp to decompose at all) and the
      ! fully-periodic 000 case (the plain-cuFFT fallback performs a purely
      ! local per-rank transform with no cross-rank exchange -
      ! src/backend/cuda/poisson_fft.f90:721-740 - so at nproc>1 it would
      ! run without erroring but compute invalid results). cuFFTMp
      ! unavailable here means the build just performed fell back to plain
      ! cuFFT.
      multi_gpu_supported_measured = .false.
    end if

    call print_rule('-')
    print '(a,i0,a,i0,a,f0.2,a,f0.2,a)', 'Real build (1 GPU, ', &
      ctx%measured_n_substeps, ' substep(s)): peak_fields measured ', &
      measured_peak_fields, ', workspace ', workspace_gib_measured, &
      ' GiB, used ', used_gib, ' GiB'
    ! Single machine-parseable line for --extensive sweeps (see
    ! scripts/memcheck_extensive_sweep.sh), so a wrapper script comparing
    ! several substep counts can grep one line per run instead of parsing
    ! the whole table.
    if (ctx%extensive_substeps > 0) &
      print '(a,i0,a,i0,a,f0.3,a,f0.3)', 'EXTENSIVE result: substeps=', &
      ctx%measured_n_substeps, ' peak_fields=', measured_peak_fields, &
      ' used_gib=', used_gib, ' pct_card=', 100._dp*used_gib/ctx%card_gib
    call print_table_header()

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      ! ng_unsupported_reason covers the base checks shared with report()'s
      ! table; multi_gpu_supported_measured on top of that additionally
      ! excludes ng>1 when the real build just fell back from cuFFTMp to
      ! plain cuFFT for a BC that needs it, a signal only this real-build
      ! path has (see its own derivation above).
      reason = ng_unsupported_reason(ctx, ng)
      if (ng > 1 .and. ctx%multi_gpu_supported .and. &
          .not. multi_gpu_supported_measured) &
        reason = 'not supported for this BC/environment'
      if (len_trim(reason) > 0) then
        print '(a,i0,a,a,a)', ' ', ng, '     (skipped: ', trim(reason), ')'
        cycle
      end if
      local_dims = [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)/ng]
      w_local = to_gib(fields_plus_halo_bytes_n(ctx, ng, measured_peak_fields))

      mirror_gib = 0._dp
      if (ctx%bc_is_100) &
        mirror_gib = to_gib(mirror_buffer_bytes_100(ctx%cdims, ng))

      overhead_term = overhead + &
                      to_gib(spectral_slab_bytes(ctx%bc_is_100, &
                                ctx%bc_is_110, ctx%cdims, ng) - spec_bytes_1) &
                      + mirror_gib
      if (ctx%use_cufftmp_known .and. ctx%use_cufftmp) &
        overhead_term = overhead_term + &
                        to_gib(context_floor_bytes(ng, .true.) - &
                               context_floor_bytes(1, .true.))
      overhead_term = overhead_term + &
                      to_gib(gpu_io_staging_bytes(local_dims, &
                                                  ctx%checkpoint_cfg, &
                                                  ctx%gpu_io_device_write))

      per_gpu = w_local + overhead_term
      verdict = classify(per_gpu, ctx%card_gib)
      call print_table_row(ng, local_dims, w_local, overhead_term, per_gpu, &
                           ctx%card_gib, verdict)
      if (ng == ctx%domain_cfg%nproc_dir(3)) requested_verdict = verdict
    end do
    call print_rule('-')
  end subroutine report_measured_table

end module m_memcheck_build
