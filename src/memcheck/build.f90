module m_memcheck_build
  !! Tier 3 of x3d2-memcheck: the real case build, its measurement and the
  !! measured report.
  use m_common, only: dp, i8, VERT
  use m_mesh, only: mesh_t
  use m_postprocess, only: compute_derived_fields, compute_pressure_vert
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_case_channel, only: case_channel_t
  use m_case_cylinder, only: case_cylinder_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t
  use m_field, only: flist_t
  use m_memcheck_scratch, only: ibm_mask_filename, validate_build_inputs, &
                                make_scratch, leave_build_scratch
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memcheck_report, only: n_gpu_list, print_rule, print_table_header, &
                               print_table_row
  use m_memcheck_estimate, only: to_gib, classify, ng_unsupported_reason, &
                                 fields_plus_halo_bytes_n, n_halo
  use m_memory_estimate, only: spectral_slab_bytes, mirror_buffer_bytes_100, &
                               peak_fields_lookup, gpu_io_staging_bytes, &
                               output_field_active, padded_halo_bytes
  implicit none

  private
  public :: run_tier3

  !> --build: do not attempt a real build whose fields only workspace alone
  !> exceeds this fraction of the FREE card memory (card_free_gib, not the
  !> card's full capacity) - leaves headroom for context + FFT scratch so
  !> the build never OOMs. Same value and role this tool's real build probe
  !> has always used.
  real(dp), parameter :: BUILD_FRACTION = 0.55_dp
  !> --build's CHECK line: how far the static estimate may be from the real
  !> measurement before it is flagged MISMATCH rather than OK. Justified by
  !> this tool's own calibration record: TGV matched to 1.5%; a real
  !> calibration bug (an LES peak_fields delta wrong by 8) showed up as an
  !> ~8% byte mismatch, well outside this tolerance - and is separately,
  !> exactly caught by the drift check below regardless of byte tolerance.
  real(dp), parameter :: CHECK_TOLERANCE = 0.05_dp

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
             ctx%device%fft_path_known() .and. &
             (.not. ctx%device%uses_distributed_fft())) then
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
                      to_gib(spectral_slab_bytes( &
                             ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, ng) &
                             - spec_bytes_1) &
                      + mirror_gib
      ! The floor takes bc_is_110 where the pre-refactor code hard-coded
      ! uses_cufftmp=.true.: equal here because the 110 case never uses the
      ! distributed FFT and its ng>1 rows are skipped above, so this term
      ! is never evaluated for it.
      if (ctx%device%fft_path_known() .and. &
          ctx%device%uses_distributed_fft()) &
        overhead_term = overhead_term + &
                        to_gib(ctx%device%overhead_floor_bytes( &
                               ng, ctx%bc_is_110) - &
                               ctx%device%overhead_floor_bytes( &
                               1, ctx%bc_is_110))
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

  subroutine build_mesh(ctx, dims_in, mesh, dims, ibm_missing)
    !! Mesh of the real build at dims_in plus the ibm_on mask pre-check;
    !! ibm_missing=.true. tells build_and_measure to stop.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: dims_in(3)
    type(mesh_t), intent(out) :: mesh
    integer, intent(out) :: dims(3)
    logical, intent(out) :: ibm_missing

    character(len=16) :: ibm_file
    logical :: ibm_file_exists

    ibm_missing = .false.
    mesh = mesh_t(dims_in, [1, 1, 1], ctx%domain_cfg%L_global, &
                  ctx%domain_cfg%BC_x, ctx%domain_cfg%BC_y, &
                  ctx%domain_cfg%BC_z, ctx%domain_cfg%stretching, &
                  ctx%domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)

    ! ibm_on triggers reading an external ibm_<BC-suffix>.bp mask file
    ! inside solver init (src/module/ibm.f90), which MPI_Aborts if that
    ! file is missing. That is correct for production xcompact, but not
    ! for a tool that should degrade gracefully - check for it here, using
    ! the same suffix construction as src/module/ibm.f90:70-75 (shared via
    ! ibm_mask_filename, also used by make_scratch to mirror the
    ! mask into the build scratch directory), and bail out before
    ! triggering the case/solver construction that would abort.
    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(mesh%grid%periodic_BC(1), &
                                   mesh%grid%periodic_BC(2), &
                                   mesh%grid%periodic_BC(3))
      inquire (file=ibm_file, exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        ibm_missing = .true.
        return
      end if
    end if
  end subroutine build_mesh

  subroutine make_flow_case(ctx, backend, mesh, host_allocator, flow_case)
    !! Build the flow case named in the input through the same case
    !! dispatch select case xcompact.f90 uses.
    type(memcheck_ctx_t), intent(in) :: ctx
    class(base_backend_t), pointer, intent(in) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    class(base_case_t), allocatable, intent(out) :: flow_case

    select case (trim(ctx%domain_cfg%flow_case_name))
    case ('channel')
      allocate (case_channel_t :: flow_case)
      flow_case = case_channel_t(backend, mesh, host_allocator)
    case ('cylinder')
      allocate (case_cylinder_t :: flow_case)
      flow_case = case_cylinder_t(backend, mesh, host_allocator)
    case ('generic')
      allocate (case_generic_t :: flow_case)
      flow_case = case_generic_t(backend, mesh, host_allocator)
    case ('tgv')
      allocate (case_tgv_t :: flow_case)
      flow_case = case_tgv_t(backend, mesh, host_allocator)
    case default
      error stop 'Undefined flow_case.'
    end select
  end subroutine make_flow_case

  subroutine drive_case(ctx, flow_case)
    !! Drive the built flow case the way run() does so the allocator reaches
    !! its high-water mark; records measured_n_substeps and the output flags.
    type(memcheck_ctx_t), intent(inout) :: ctx
    class(base_case_t), intent(inout) :: flow_case

    type(flist_t), allocatable :: curr(:), deriv(:)
    integer :: i, n_substeps, i_substep

    allocate (curr(flow_case%solver%time_integrator%nvars))
    allocate (deriv(flow_case%solver%time_integrator%nvars))
    curr(1)%ptr => flow_case%solver%u
    curr(2)%ptr => flow_case%solver%v
    curr(3)%ptr => flow_case%solver%w
    do i = 1, flow_case%solver%nspecies
      curr(3 + i)%ptr => flow_case%solver%species(i)%ptr
    end do

    ! Mirrors run()'s pre-loop call to the case-specific postprocess hook.
    call flow_case%postprocess(0, 0._dp)

    ! One real sub-stage drives the allocator to the same high-water mark a
    ! full time step would reach; field values are irrelevant to memory.
    ! --extensive <n> overrides the count (default 1) so the peak can be
    ! re-checked against more sub-stages, e.g. to confirm no later-only
    ! allocation path changes it. This only re-checks the allocator's own
    ! pool (structurally leak-free - next_id is monotonic, get_block/
    ! release_block just recycle a free list, see src/allocator.f90): it
    ! never repeats the per-iteration I/O paths (compute_pressure_vert,
    ! compute_derived_fields, io_mgr%update_stats, snapshot/checkpoint
    ! writes) that a real xcompact run would hit on every iteration and
    ! where an actual leak is far more likely to live.
    n_substeps = 1
    if (ctx%extensive_substeps > 0) n_substeps = ctx%extensive_substeps
    do i_substep = 1, n_substeps
      call flow_case%substep(curr, deriv, i_substep)
    end do
    ctx%measured_n_substeps = n_substeps

    ! Mirrors run()'s per-iteration postprocessing calls (base_case.f90,
    ! immediately after the sub-stage loop): keep_pressure and
    ! vorticity/Q-criterion allocations are gated HERE, not inside
    ! postprocess() above, so they must be driven separately to be
    ! captured - postprocess(0,.) alone does not reach them.
    if (flow_case%solver%keep_pressure) &
      call compute_pressure_vert(flow_case%solver)
    ctx%output_vorticity = output_field_active( &
                           flow_case%io_mgr%snapshot_mgr%config, 'vorticity')
    ctx%output_qcriterion = output_field_active( &
                            flow_case%io_mgr%snapshot_mgr%config, 'qcriterion')
    if (ctx%output_vorticity .or. ctx%output_qcriterion) &
      call compute_derived_fields(flow_case%solver, ctx%output_vorticity, &
                                  ctx%output_qcriterion)
  end subroutine drive_case

  subroutine build_and_measure(ctx, dims_in, npeak, dev_used, ibm_missing)
    !! Build the real flow case at dims_in on this single GPU (via the same
    !! case-dispatch select case xcompact.f90 uses), run one (or, under
    !! --extensive <n>, n) substep(s) to drive the allocator to its
    !! work-field high-water mark, and return that mark (npeak) plus the
    !! absolute device memory in use (dev_used).
    !!
    !! --extensive re-checks this same high-water mark against more RK
    !! sub-stages (drive_case says what that does and does not cover, in
    !! particular that it cannot catch a leak in the per-iteration I/O
    !! paths). The GPU-aware ADIOS2 staging buffer is therefore added
    !! analytically (gpu_io_staging_bytes) rather than measured here. See
    !! scripts/memcheck_extensive_sweep.sh
    !! (this tool's own sweep, allocator-path regression guard only) vs.
    !! scripts/memcheck_extensive_xcompact.sh (a real xcompact run polled
    !! externally via nvidia-smi, which CAN catch an I/O-path leak).
    !!
    !! Driving the real flow_case_t (case_init, one postprocess(0,.) call,
    !! one substep) rather than a bespoke transeq/step/pressure_correction
    !! sequence means case-specific persistent allocations are captured the
    !! same way a real xcompact run captures them: channel/cylinder's
    !! lazily-allocated BC ghost blocks (define_BC_channel/cylinder), and
    !! anything gated by the input's I/O config (keep_pressure,
    !! vorticity/Q-criterion derived fields via the per-iteration
    !! compute_pressure_vert/compute_derived_fields calls drive_case
    !! mirrors - these are NOT part of postprocess()).
    !!
    !! The real case constructors re-read get_argument(1) (this process's
    !! own CLI argument) rather than taking domain_cfg as a parameter; this
    !! already resolves to the same input file read_config() parsed above -
    !! see parse_args' comment on why the input path must stay positional
    !! argument 1.
    !!
    !! case_init also runs the input's own restart logic (io_mgr%is_restart)
    !! before this subroutine gets control: if a checkpoint file matching
    !! the input's restart path already exists in the CWD, this will
    !! measure a restarted state instead of the input's initial conditions.
    !! This mirrors production behaviour, not a tool-specific choice.
    type(memcheck_ctx_t), target, intent(inout) :: ctx
    integer, intent(in) :: dims_in(3)
    integer, intent(out) :: npeak
    real(dp), intent(out) :: dev_used
    logical, intent(out) :: ibm_missing

    type(mesh_t), target :: mesh
    class(allocator_t), pointer :: allocator
    class(base_backend_t), pointer :: backend
    class(base_case_t), allocatable :: flow_case
    integer :: dims(3)
    type(allocator_t), target :: host_allocator
    integer(i8) :: free_b, total_b

    npeak = 0
    dev_used = 0._dp
    ctx%output_vorticity = .false.
    ctx%output_qcriterion = .false.

    call build_mesh(ctx, dims_in, mesh, dims, ibm_missing)
    if (ibm_missing) return

    call ctx%device%make_backend(mesh, dims, allocator, backend)
    host_allocator = allocator_t(dims, ctx%device%sz())

    call make_flow_case(ctx, backend, mesh, host_allocator, flow_case)

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_measured_table's safety gating on
    ! the 100 case, since there is no static "is cuFFTMp available" query.
    call ctx%device%capture_solver_info(flow_case)

    call drive_case(ctx, flow_case)

    npeak = allocator%next_id
    call ctx%device%mem_info(total_b, free_b)
    dev_used = to_gib(total_b - free_b)
  end subroutine build_and_measure

  subroutine run_tier3(ctx)
    !! Tier 3: a REAL case build + one substep on 1 GPU, measured with
    !! cudaMemGetInfo - the last-resort ground truth, now run in-process
    !! (see the top-of-file note on why this removes the old nested mpirun
    !! limitation). Skips the build entirely (the estimate above stands)
    !! if the fields only workspace alone already exceeds BUILD_FRACTION of
    !! what is currently FREE on the card (card_free_gib, not the card's
    !! full capacity - see its declaration), matching this tool's real
    !! build probe, or if ibm_on=T and the matching mask file is not
    !! present in the working directory. The real build itself runs inside
    !! a throwaway x3d2-memcheck-build.<pid> scratch directory
    !! (make_scratch/leave_build_scratch; make_scratch's docstring says why,
    !! when it skips the build, and what an unexpected error stop after its
    !! chdir leaves behind). validate_build_inputs pre-validates the known
    !! causes up front - a missing input file, an unsupported flow case, or
    !! (ibm_on=T) a missing mask file - before touching the filesystem.
    type(memcheck_ctx_t), intent(inout) :: ctx
    real(dp) :: used_gib, workspace_gib_measured, ws_guess
    real(dp) :: pct_error, estimate_ng1_excl_io_gib, residual_gib
    character(len=8) :: check_word
    integer :: measured_peak_fields
    logical :: ibm_missing, build_scratch_ok
    character(len=12) :: measured_verdict

    ! --build gate: skip the real build (after printing why) when the fields
    ! only workspace alone already exceeds BUILD_FRACTION of the FREE card
    ! memory.
    ws_guess = to_gib(fields_plus_halo_bytes_n(ctx, 1, ctx%peak_fields))
    ! fields_plus_halo_bytes_n includes the halo term; the historical
    ! BUILD_FRACTION guard was calibrated against fields alone, so subtract
    ! it back out here rather than changing the (already-validated)
    ! threshold itself.
    ws_guess = ws_guess - &
               to_gib(padded_halo_bytes(ctx%gdims, ctx%device%sz(), n_halo))
    if (ws_guess >= BUILD_FRACTION*ctx%card_free_gib) then
      call print_rule('-')
      print '(a,f0.2,a,f0.1,a,f0.2,a)', 'Real build skipped: fields only &
        &workspace ', ws_guess, ' GiB exceeds ', 100._dp*BUILD_FRACTION, &
        '% of the FREE card memory (', ctx%card_free_gib, ' GiB free); the &
        &estimate above stands.'
      return
    end if

    ! Printed under --build only, right before the real build: the only mode
    ! where a second cuFFTMp plan lifecycle follows in the same process, so
    ! a non-zero residual would mean the real build's own measurement is
    ! inflated by whatever Tier 2 left behind.
    if (ctx%device%query_residual_gib(residual_gib)) then
      call print_rule('-')
      print '(a,f0.2,a)', 'Residual after plan query: ', &
        residual_gib, ' GiB (should be ~0 if cufftDestroy &
        &released it; if not, the real build below may double-reserve the &
        &NVSHMEM heap).'
    end if

    call validate_build_inputs(ctx, build_scratch_ok)
    if (.not. build_scratch_ok) return
    call make_scratch(ctx, build_scratch_ok)
    if (.not. build_scratch_ok) return

    call build_and_measure(ctx, [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)], &
                           measured_peak_fields, used_gib, ibm_missing)
    call leave_build_scratch(ctx)
    if (ibm_missing) then
      call print_rule('-')
      print '(a)', 'Real build skipped: ibm_on=T but the matching &
        &ibm_<BC-suffix>.bp mask file was not found in the working &
        &directory; the estimate above stands.'
      return
    end if

    call check_peak_fields_table(ctx, measured_peak_fields)

    workspace_gib_measured = &
      to_gib(fields_plus_halo_bytes_n(ctx, 1, measured_peak_fields))

    call report_measured_table(ctx, measured_peak_fields, used_gib, &
                               workspace_gib_measured, measured_verdict)
    ctx%final_verdict = measured_verdict

    ! --build's CHECK line: the static ng=1 estimate against the measured
    ! device memory in use, OK within CHECK_TOLERANCE and MISMATCH beyond it.
    ! The real build above performs no snapshot/checkpoint write, so the
    ! GPU-aware IO staging term (analytical only, never measured here) is
    ! excluded from both the estimate compared and the printed CHECK line.
    estimate_ng1_excl_io_gib = ctx%estimate_ng1_gib - ctx%io_staging_ng1_gib
    pct_error = 100._dp*(estimate_ng1_excl_io_gib - used_gib)/used_gib
    if (abs(estimate_ng1_excl_io_gib - used_gib) <= &
        CHECK_TOLERANCE*used_gib) then
      check_word = 'OK'
    else
      check_word = 'MISMATCH'
    end if
    print '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a,a,a)', 'CHECK ng=1: estimated ', &
      estimate_ng1_excl_io_gib, ' GiB, measured ', used_gib, ' GiB (', &
      pct_error, '%) - ', trim(check_word), ' (tolerance 5%)'
    if (ctx%io_staging_ng1_gib > 0._dp) &
      print '(a,f0.1,a)', '  (GPU-aware IO staging ', &
        ctx%io_staging_ng1_gib*1024._dp, ' MiB excluded from CHECK: the real &
        &build performs no snapshot/checkpoint write.)'
  end subroutine run_tier3

end module m_memcheck_build
