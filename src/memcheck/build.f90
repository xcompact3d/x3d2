module m_memcheck_build
  !! The real build of x3d2-memcheck (--build): the real case build, its
  !! measurement and the measured report.
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
                               print_table_row, print_requested_line
  use m_memcheck_estimate, only: to_gib, num_gpus_unsupported_reason, &
                                 requested_unsupported_reason, &
                                 estimate_for_num_gpus, &
                                 fields_plus_halo_bytes_n, n_halo, basis_name
  use m_memory_estimate, only: peak_fields_lookup, output_field_active, &
                               padded_halo_bytes
  implicit none

  private
  public :: run_build_mode

  !> --build skips the real build when the fields alone take more than this
  !> fraction of the FREE card memory (card_free_gib), which leaves headroom
  !> for the context and the FFT scratch.
  real(dp), parameter :: BUILD_FRACTION = 0.55_dp
  !> The distance the static estimate may be from the measurement before the
  !> CHECK line of --build says MISMATCH.
  real(dp), parameter :: CHECK_TOLERANCE = 0.05_dp

contains

  subroutine check_peak_fields_table(ctx, measured)
    !! Warns when peak_fields_lookup differs from the measured peak_fields.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: measured
    integer :: predicted

    if (ctx%irank /= 0) return

    predicted = peak_fields_lookup(trim(ctx%domain_cfg%flow_case_name), &
                                   ctx%solver_cfg%time_intg, &
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
    !! Prints the measured table: the num_gpus=1 measurement scaled to every
    !! GPU count by estimate_for_num_gpus. A supported requested num_gpus
    !! outside n_gpu_list gets its own row. requested_verdict is the verdict
    !! of the requested num_gpus.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: measured_peak_fields
    real(dp), intent(in) :: used_gib, workspace_gib_measured
    character(len=12), intent(out) :: requested_verdict

    integer :: k, num_gpus, basis
    real(dp) :: per_gpu_gib, workspace_gib, overhead_gib
    logical :: any_extrapolated
    character(len=12) :: verdict
    character(len=96) :: reason

    requested_verdict = 'DOES_NOT_FIT'

    call print_rule('-')
    print '(a,i0,a,i0,a,f0.2,a,f0.2,a)', 'Real build (1 GPU, ', &
      ctx%measured_n_substeps, ' substep(s)): peak_fields measured ', &
      measured_peak_fields, ', workspace ', workspace_gib_measured, &
      ' GiB, used ', used_gib, ' GiB'
    ! One line for scripts that compare --extensive runs.
    if (ctx%extensive_substeps > 0) &
      print '(a,i0,a,i0,a,f0.3,a,f0.3)', 'EXTENSIVE result: substeps=', &
      ctx%measured_n_substeps, ' peak_fields=', measured_peak_fields, &
      ' used_gib=', used_gib, ' pct_card=', 100._dp*used_gib/ctx%card_gib
    call print_table_header()

    any_extrapolated = .false.
    do k = 1, size(n_gpu_list) + 1
      num_gpus = ctx%domain_cfg%nproc_dir(3)
      if (k <= size(n_gpu_list)) num_gpus = n_gpu_list(k)
      if (k > size(n_gpu_list) .and. &
          (any(n_gpu_list == num_gpus) .or. &
           len_trim(requested_unsupported_reason(ctx)) > 0)) exit
      reason = num_gpus_unsupported_reason(ctx, num_gpus)
      if (len_trim(reason) > 0) then
        print '(a,i0,a,a,a)', ' ', num_gpus, '     (skipped: ', &
          trim(reason), ')'
        cycle
      end if
      call estimate_for_num_gpus(ctx, num_gpus, per_gpu_gib, basis, verdict, &
                                 workspace_gib, overhead_gib, &
                                 used_gib=used_gib, npeak=measured_peak_fields)
      call print_table_row(num_gpus, [ctx%gdims(1), ctx%gdims(2), &
                                      ctx%gdims(3)/num_gpus], workspace_gib, &
                           overhead_gib, per_gpu_gib, ctx%card_gib, &
                           basis_name(basis), verdict)
      if (num_gpus == ctx%domain_cfg%nproc_dir(3)) requested_verdict = verdict
      if (num_gpus > 1) any_extrapolated = .true.
    end do
    call print_rule('-')
    if (any_extrapolated) &
      print '(a)', '(num_gpus>1 rows extrapolate the single-GPU measurement &
        &by 1/num_gpus.)'
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

    ! A missing mask file makes the solver init abort: stop before it.
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

    ! One sub-stage reaches the allocator high-water mark of a time step.
    ! --extensive <n> drives n sub-stages, which re-checks the allocator
    ! only, not the per-iteration IO paths.
    n_substeps = 1
    if (ctx%extensive_substeps > 0) n_substeps = ctx%extensive_substeps
    do i_substep = 1, n_substeps
      call flow_case%substep(curr, deriv, i_substep)
    end do
    ctx%measured_n_substeps = n_substeps

    ! The per-iteration postprocessing of run(): keep_pressure and the
    ! derived fields are allocated here, not in postprocess().
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
    !! Builds the real flow case at dims_in on this GPU, drives one substep
    !! (n under --extensive) and returns the allocator mark npeak and the
    !! device memory in use dev_used. The real case constructors re-read
    !! argument 1 as the input file, and case_init runs the restart logic of
    !! the input, as in production.
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

    ! The solver construction settled the FFT path: record it.
    call ctx%device%capture_solver_info(flow_case)

    call drive_case(ctx, flow_case)

    npeak = allocator%next_id
    call ctx%device%mem_info(total_b, free_b)
    dev_used = to_gib(total_b - free_b)
  end subroutine build_and_measure

  subroutine run_build_mode(ctx)
    !! The real build: a case build and one substep on 1 GPU in this process,
    !! measured with cudaMemGetInfo. The estimate stands when the fields
    !! alone take more than BUILD_FRACTION of the free card memory or the
    !! inputs are not valid. The build runs in a scratch directory because a
    !! case build writes monitoring.csv in the working directory.
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
    ! inflated by whatever the cuFFTMp probe left behind.
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
    ! An unsupported nproc_dir stays UNSUPPORTED. The real build can make a
    ! request unsupported (cuFFTMp fell back): print that line now.
    if (len_trim(requested_unsupported_reason(ctx)) == 0) then
      ctx%final_verdict = measured_verdict
    else if (trim(ctx%final_verdict) /= 'UNSUPPORTED') then
      call print_requested_line(ctx)
    end if

    ! The CHECK line compares the num_gpus=1 estimate without the IO staging
    ! term, which the real build never measures, with the measured memory.
    estimate_ng1_excl_io_gib = ctx%estimate_ng1_gib - ctx%io_staging_ng1_gib
    pct_error = 100._dp*(estimate_ng1_excl_io_gib - used_gib)/used_gib
    if (abs(estimate_ng1_excl_io_gib - used_gib) <= &
        CHECK_TOLERANCE*used_gib) then
      check_word = 'OK'
    else
      check_word = 'MISMATCH'
    end if
    print '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a,a,a)', &
      'CHECK num_gpus=1: estimated ', &
      estimate_ng1_excl_io_gib, ' GiB, measured ', used_gib, ' GiB (', &
      pct_error, '%) - ', trim(check_word), ' (tolerance 5%)'
    if (ctx%io_staging_ng1_gib > 0._dp) &
      print '(a,f0.1,a)', '  (GPU-aware IO staging ', &
        ctx%io_staging_ng1_gib*1024._dp, ' MiB excluded from CHECK: the real &
        &build performs no snapshot/checkpoint write.)'
  end subroutine run_build_mode

end module m_memcheck_build
