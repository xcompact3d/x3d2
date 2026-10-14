module m_memcheck_build
  !! Tier 3 of x3d2-memcheck: the real case build, its measurement and the
  !! measured report.
  use m_common, only: dp, i8, VERT
  use cudafor, only: cudaMemGetInfo, cuda_count_kind
  use m_mesh, only: mesh_t
  use m_cuda_common, only: SZ
  use m_postprocess, only: compute_derived_fields, compute_pressure_vert
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_case_channel, only: case_channel_t
  use m_case_cylinder, only: case_cylinder_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t
  use m_field, only: flist_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t
  use m_memcheck_scratch, only: ibm_mask_filename
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memcheck_report, only: n_gpu_list, print_rule, print_table_header, &
                               print_table_row
  use m_memcheck_estimate, only: to_gib, classify, ng_unsupported_reason, &
                                 fields_plus_halo_bytes_n
  use m_cuda_memory_estimate, only: context_floor_bytes
  use m_memory_estimate, only: spectral_slab_bytes, mirror_buffer_bytes_100, &
                               peak_fields_lookup, gpu_io_staging_bytes, &
                               output_field_active
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

  subroutine make_cuda_backend(mesh, dims, cuda_allocator, host_allocator, &
                               cuda_backend, allocator, backend)
    !! Build the device and host allocators and the CUDA backend; the
    !! caller owns all three so the pointers cannot outlive them.
    type(mesh_t), target, intent(inout) :: mesh
    integer, intent(in) :: dims(3)
    type(cuda_allocator_t), target, intent(inout) :: cuda_allocator
    type(allocator_t), target, intent(inout) :: host_allocator
    type(cuda_backend_t), target, intent(inout) :: cuda_backend
    class(allocator_t), pointer, intent(out) :: allocator
    class(base_backend_t), pointer, intent(out) :: backend

    cuda_allocator = cuda_allocator_t(dims, SZ)
    allocator => cuda_allocator
    host_allocator = allocator_t(dims, SZ)
    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend
  end subroutine make_cuda_backend

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
    !! sub-stages, but only re-exercises the allocator's own pool - it does
    !! NOT repeat the per-iteration I/O paths (compute_pressure_vert,
    !! compute_derived_fields, io_mgr%update_stats, snapshot/checkpoint
    !! writes) that a real xcompact run hits every iteration, so it cannot
    !! catch a leak in those paths. The GPU-aware ADIOS2 staging buffer is
    !! therefore added analytically (gpu_io_staging_bytes) rather than
    !! measured here. See scripts/memcheck_extensive_sweep.sh
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
    !! compute_pressure_vert/compute_derived_fields calls mirrored below,
    !! immediately after substep - these are NOT part of postprocess()).
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
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: dims_in(3)
    integer, intent(out) :: npeak
    real(dp), intent(out) :: dev_used
    logical, intent(out) :: ibm_missing

    type(mesh_t), target :: mesh
    class(allocator_t), pointer :: allocator
    class(base_backend_t), pointer :: backend
    class(base_case_t), allocatable :: flow_case
    integer :: dims(3)
    type(cuda_allocator_t), target :: cuda_allocator
    type(cuda_backend_t), target :: cuda_backend
    type(allocator_t), target :: host_allocator
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    npeak = 0
    dev_used = 0._dp
    ctx%output_vorticity = .false.
    ctx%output_qcriterion = .false.

    call build_mesh(ctx, dims_in, mesh, dims, ibm_missing)
    if (ibm_missing) return

    call make_cuda_backend(mesh, dims, cuda_allocator, host_allocator, &
                           cuda_backend, allocator, backend)

    call make_flow_case(ctx, backend, mesh, host_allocator, flow_case)

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_measured_table's safety gating on
    ! the 100 case, since there is no static "is cuFFTMp available" query.
    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      ctx%use_cufftmp = pf%use_cufftmp
      ctx%use_cufftmp_known = .true.
    end select

    call drive_case(ctx, flow_case)

    npeak = allocator%next_id
    ierr = cudaMemGetInfo(free_b, total_b)
    dev_used = to_gib(int(total_b - free_b, i8))
  end subroutine build_and_measure

end module m_memcheck_build
