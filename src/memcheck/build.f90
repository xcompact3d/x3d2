module m_memcheck_build
  !! The real build of x3d2-memcheck (--build): the real case build, its
  !! measurement and the measured report.
  use m_common, only: dp, i8, VERT
  use m_logger, only: log_info, log_warn
  use m_mesh, only: mesh_t
  use m_postprocess, only: compute_derived_fields, compute_pressure_vert
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_case_factory, only: allocate_flow_case, is_supported_flow_case
  use m_field, only: flist_t
  use m_monitoring, only: monitoring_default_enabled
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memcheck_report, only: n_gpu_list, print_rule, print_table_header, &
                               print_table_row
  use m_memcheck_estimate, only: to_gib, num_gpus_unsupported_reason, &
                                 requested_unsupported_reason, &
                                 estimate_for_num_gpus, &
                                 fields_plus_halo_bytes_n, basis_name
  use m_memory_estimate, only: peak_fields_lookup, output_field_active
  implicit none

  private
  public :: run_build_mode

  !> The distance the static estimate may be from the measurement before the
  !> CHECK line of --build says MISMATCH.
  real(dp), parameter :: CHECK_TOLERANCE = 0.05_dp

contains

  pure function ibm_mask_filename(px, py, pz) result(fname)
    !! The ibm_<BC-suffix>.bp mask filename for the periodic flags, as in
    !! src/module/ibm.f90.
    logical, intent(in) :: px, py, pz
    character(len=16) :: fname
    character(len=3) :: bc_suffix

    bc_suffix(1:1) = '0'
    if (.not. px) bc_suffix(1:1) = '1'
    bc_suffix(2:2) = '0'
    if (.not. py) bc_suffix(2:2) = '1'
    bc_suffix(3:3) = '0'
    if (.not. pz) bc_suffix(3:3) = '1'
    fname = "ibm_"//bc_suffix//".bp"
  end function ibm_mask_filename

  subroutine validate_build_inputs(ctx, ok)
    !! Checks the known causes of a failed real build: a missing input file,
    !! a flow case without a dispatch, a missing ibm mask. ok is false and
    !! the skip line is printed on a failure.
    type(memcheck_ctx_t), intent(in) :: ctx
    logical, intent(out) :: ok

    character(len=16) :: ibm_file
    logical :: ibm_file_exists, input_exists

    ok = .false.
    inquire (file=trim(ctx%input_path), exist=input_exists)
    if (.not. input_exists) then
      call log_info('Real build skipped: input file not found: '// &
                    trim(ctx%input_path))
      return
    end if

    if (.not. is_supported_flow_case(ctx%domain_cfg%flow_case_name)) then
      call log_info("Real build skipped: flow case '"// &
                    trim(ctx%domain_cfg%flow_case_name)// &
                    "' has no dispatch in x3d2-memcheck")
      return
    end if

    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(ctx%periodic_x, ctx%periodic_y, &
                                   ctx%periodic_z)
      inquire (file=trim(ibm_file), exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        call log_info('Real build skipped: ibm_on=T but the matching &
          &ibm_<BC-suffix>.bp mask file was not found in the &
          &working directory; the estimate above stands.')
        return
      end if
    end if
    ok = .true.
  end subroutine validate_build_inputs

  subroutine check_peak_fields_table(ctx, measured)
    !! Warns when peak_fields_lookup differs from the measured peak_fields.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: measured
    integer :: predicted
    character(len=256) :: line_buf

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
      write (line_buf, '(a,i0,a,i0,a)') &
        'peak_fields_lookup predicted ', predicted, &
        ' but the real build measured ', measured, &
        ' - m_memory_estimate''s table is stale for this config and &
        &should be recalibrated (see src/memory_estimate.f90).'
      call log_warn(trim(line_buf))
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
    character(len=256) :: line_buf

    requested_verdict = 'DOES_NOT_FIT'

    call print_rule('-')
    write (line_buf, '(a,i0,a,i0,a,f0.2,a,f0.2,a)') &
      'Real build (1 GPU, ', &
      ctx%measured_n_substeps, ' substep(s)): peak_fields measured ', &
      measured_peak_fields, ', workspace ', workspace_gib_measured, &
      ' GiB, used ', used_gib, ' GiB'
    call log_info(trim(line_buf))
    ! One line for scripts that compare --extensive runs.
    if (ctx%extensive_substeps > 0) then
      write (line_buf, '(a,i0,a,i0,a,f0.3,a,f0.3)') &
        'EXTENSIVE result: substeps=', &
        ctx%measured_n_substeps, ' peak_fields=', measured_peak_fields, &
        ' used_gib=', used_gib, ' pct_card=', 100._dp*used_gib/ctx%card_gib
      call log_info(trim(line_buf))
    end if
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
        write (line_buf, '(a,i0,a,a,a)') &
          ' ', num_gpus, '     (skipped: ', &
          trim(reason), ')'
        call log_info(trim(line_buf))
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
      call log_info('(num_gpus>1 rows extrapolate the single-GPU measurement &
        &by 1/num_gpus.)')
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

  subroutine build_and_measure(ctx, dims_in, npeak, free_end_bytes, &
                               ibm_missing)
    !! Builds the real flow case at dims_in on this GPU, drives one substep
    !! (n under --extensive) and returns the allocator mark npeak and the
    !! free device memory free_end_bytes at that point. The real case
    !! constructors re-read argument 1 as the input file, and case_init
    !! runs the restart logic of the input, as in production.
    type(memcheck_ctx_t), target, intent(inout) :: ctx
    integer, intent(in) :: dims_in(3)
    integer, intent(out) :: npeak
    integer(i8), intent(out) :: free_end_bytes
    logical, intent(out) :: ibm_missing

    type(mesh_t), target :: mesh
    class(allocator_t), pointer :: allocator
    class(base_backend_t), pointer :: backend
    class(base_case_t), allocatable :: flow_case
    integer :: dims(3)
    type(allocator_t), target :: host_allocator
    integer(i8) :: total_b

    npeak = 0
    free_end_bytes = 0_i8
    ctx%output_vorticity = .false.
    ctx%output_qcriterion = .false.

    call build_mesh(ctx, dims_in, mesh, dims, ibm_missing)
    if (ibm_missing) return

    call ctx%device%make_backend(mesh, dims, allocator, backend)
    host_allocator = allocator_t(dims, ctx%device%sz())

    ! The case build writes no monitoring.csv in the working directory.
    monitoring_default_enabled = .false.
    call allocate_flow_case(ctx%domain_cfg%flow_case_name, backend, mesh, &
                            host_allocator, flow_case)

    ! The solver construction settled the FFT path: record it.
    call ctx%device%capture_solver_info(flow_case)

    call drive_case(ctx, flow_case)

    npeak = allocator%next_id
    call ctx%device%mem_info(total_b, free_end_bytes)
  end subroutine build_and_measure

  subroutine run_real_build(ctx, measured_peak_fields, used_gib, ibm_missing)
    !! Runs the real build in this process, after the FFT probe has destroyed
    !! its plans. used_gib is the card memory the build holds: the context at
    !! the start plus the fall of the free memory since then.
    type(memcheck_ctx_t), target, intent(inout) :: ctx
    integer, intent(out) :: measured_peak_fields
    real(dp), intent(out) :: used_gib
    logical, intent(out) :: ibm_missing

    integer(i8) :: free_end_bytes

    call build_and_measure(ctx, [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)], &
                           measured_peak_fields, free_end_bytes, ibm_missing)
    used_gib = 0._dp
    if (.not. ibm_missing) &
      used_gib = to_gib(ctx%device%footprint_bytes(free_end_bytes))
  end subroutine run_real_build

  subroutine run_build_mode(ctx)
    !! The real build: a case build and one substep on 1 GPU in this process,
    !! measured with cudaMemGetInfo. The estimate stands when the
    !! num_gpus=1 row is DOES_NOT_FIT or the inputs are not valid.
    type(memcheck_ctx_t), intent(inout) :: ctx
    real(dp) :: used_gib, workspace_gib_measured
    real(dp) :: pct_error, estimate_ng1_excl_io_gib
    character(len=8) :: check_word
    integer :: measured_peak_fields
    logical :: ibm_missing, inputs_ok
    character(len=12) :: measured_verdict
    character(len=256) :: line_buf

    if (trim(ctx%ng1_verdict) == 'DOES_NOT_FIT') then
      call print_rule('-')
      call log_info('Real build skipped: the num_gpus=1 estimate is &
        &DOES_NOT_FIT; the estimate above stands.')
      return
    end if

    call validate_build_inputs(ctx, inputs_ok)
    if (.not. inputs_ok) return

    call run_real_build(ctx, measured_peak_fields, used_gib, ibm_missing)
    if (ibm_missing) then
      call print_rule('-')
      call log_info('Real build skipped: ibm_on=T but the matching &
        &ibm_<BC-suffix>.bp mask file was not found in the working &
        &directory; the estimate above stands.')
      return
    end if

    call check_peak_fields_table(ctx, measured_peak_fields)

    workspace_gib_measured = &
      to_gib(fields_plus_halo_bytes_n(ctx, 1, measured_peak_fields))

    call report_measured_table(ctx, measured_peak_fields, used_gib, &
                               workspace_gib_measured, measured_verdict)
    ! An unsupported nproc_dir stays UNSUPPORTED.
    if (len_trim(requested_unsupported_reason(ctx)) == 0) &
      ctx%final_verdict = measured_verdict

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
    write (line_buf, '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a,a,a)') &
      'CHECK num_gpus=1: estimated ', &
      estimate_ng1_excl_io_gib, ' GiB, measured ', used_gib, ' GiB (', &
      pct_error, '%) - ', trim(check_word), ' (tolerance 5%)'
    call log_info(trim(line_buf))
    if (ctx%io_staging_ng1_gib > 0._dp) then
      write (line_buf, '(a,f0.1,a)') &
        '  (GPU-aware IO staging ', &
        ctx%io_staging_ng1_gib*1024._dp, ' MiB excluded from CHECK: the real &
        &build performs no snapshot/checkpoint write.)'
      call log_info(trim(line_buf))
    end if
  end subroutine run_build_mode

end module m_memcheck_build
