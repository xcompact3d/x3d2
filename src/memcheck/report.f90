module m_memcheck_report
  !! Printed report of x3d2-memcheck: header, tables and verdict lines.
  use m_common, only: dp, i8, is_sp
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memory_estimate, only: gpu_io_staging_bytes
  use m_memcheck_estimate, only: to_gib, ng_unsupported_reason, &
                                 estimate_for_ng
  implicit none

  private
  public :: report, n_gpu_list, print_rule, print_table_header, &
            print_table_row

  !> Candidate GPU counts to scan for "smallest ng that fits". The CUDA
  !> backend only supports a Z-only pencil decomposition (nproc_dir=
  !> [1,1,ng]) today.
  integer, parameter :: n_gpu_list(4) = [1, 2, 4, 8]

contains

  subroutine print_rule(fill)
    !! One 60 character separator line: '=' around the report header, '-'
    !! between its sections.
    character(len=1), intent(in) :: fill

    print '(a)', repeat(fill, 60)
  end subroutine print_rule

  subroutine print_table_header()
    character(len=14) :: hdr_grid

    ! hdr_grid is a fixed-length variable (not a literal passed straight to
    ! a14) so the assignment's right-hand-side blank-padding left-justifies
    ! it, matching adjustl(local_grid_str) in print_table_row below - a
    ! literal would instead be right-justified by the a14 edit descriptor,
    ! misaligning the column under the data rows.
    hdr_grid = 'local grid'
    print '(a,a4,a,a14,a,a9,a,a8,a,a8,a,a6,a,a)', ' ', 'GPUs', ' | ', &
      hdr_grid, ' | ', 'workspace', ' | ', 'overhead', ' | ', 'per-GPU', &
      ' | ', '%', ' | ', 'verdict'
  end subroutine print_table_header

  subroutine print_table_row(ng, local_dims, workspace_gib, overhead_gib, &
                             per_gpu_gib, pct_denominator_gib, verdict)
    integer, intent(in) :: ng, local_dims(3)
    real(dp), intent(in) :: workspace_gib, overhead_gib, per_gpu_gib, &
                            pct_denominator_gib
    character(*), intent(in) :: verdict

    character(len=14) :: local_grid_str

    local_grid_str = ''
    write (local_grid_str, '(i0,a,i0,a,i0)') local_dims(1), 'x', &
      local_dims(2), 'x', local_dims(3)
    print '(a,i4,a,a,a,f9.2,a,f8.2,a,f8.2,a,f5.1,a,a)', ' ', ng, ' | ', &
      adjustl(local_grid_str), ' | ', workspace_gib, ' | ', overhead_gib, &
      ' | ', per_gpu_gib, ' | ', 100._dp*per_gpu_gib/pct_denominator_gib, &
      '% | ', trim(verdict)
  end subroutine print_table_row

  subroutine print_io_staging_line(ctx)
    !! The GPU-aware IO staging line of the report header: the term at ng=1
    !! and what causes it, or why there is none.
    type(memcheck_ctx_t), intent(in) :: ctx
    real(dp) :: mib
    logical :: unit_stride, snapshot_active, checkpoint_active
    character(len=64) :: what, io_none_reason
    integer(i8) :: io_ng1_bytes

    ! GPU-aware IO staging: computed directly at ng=1 (local_dims == gdims
    ! there) so it is available here, ahead of the per-ng table below.
    io_ng1_bytes = gpu_io_staging_bytes(ctx%gdims, ctx%checkpoint_cfg, &
                                        ctx%gpu_io_device_write)
    unit_stride = all(ctx%checkpoint_cfg%output_stride == 1)
    snapshot_active = ctx%checkpoint_cfg%snapshot_freq > 0 .and. unit_stride
    checkpoint_active = ctx%checkpoint_cfg%checkpoint_freq > 0
    if (io_ng1_bytes > 0_i8) then
      mib = to_gib(io_ng1_bytes)*1024._dp
      if (snapshot_active .and. checkpoint_active) then
        what = 'snapshot at unit stride + checkpoint'
      else if (checkpoint_active) then
        what = 'checkpoint'
      else
        what = 'snapshot at unit stride'
      end if
      if (.not. checkpoint_active .and. ctx%checkpoint_cfg%snapshot_sp .and. &
          .not. is_sp .and. snapshot_active) what = trim(what)//' (sp)'
      print '(a,f0.1,a,a,a,a,a)', 'GPU-aware IO staging: ', mib, &
        ' MiB/GPU at ng=1 (write mode ', trim(ctx%gpu_io_mode_name), ', ', &
        trim(what), ')'
    else
      if (.not. ctx%gpu_io_device_write) then
        io_none_reason = trim(ctx%gpu_io_reason)
      else if (ctx%checkpoint_cfg%snapshot_freq > 0 .and. .not. unit_stride &
               .and. ctx%checkpoint_cfg%checkpoint_freq == 0) then
        io_none_reason = 'snapshot striding falls back to host path'
      else
        io_none_reason = 'no unit-stride snapshot and no checkpoint enabled'
      end if
      print '(a,a,a)', 'GPU-aware IO staging: none (', trim(io_none_reason), &
        ')'
    end if
  end subroutine print_io_staging_line

  subroutine print_requested_line(ctx)
    !! The verdict line for the nproc_dir the input file asks for (or why
    !! that decomposition is not supported); sets final_verdict.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer :: requested_ng
    real(dp) :: requested_gib
    logical :: exact
    character(len=12) :: requested_verdict
    character(len=96) :: reason

    requested_ng = ctx%domain_cfg%nproc_dir(3)
    reason = ng_unsupported_reason(ctx, requested_ng)
    if (ctx%domain_cfg%nproc_dir(1) /= 1 .or. &
        ctx%domain_cfg%nproc_dir(2) /= 1) &
      reason = 'only nproc_dir = [1,1,ng] is supported by the CUDA backend'
    if (len_trim(reason) > 0) then
      print '(a,i0,a,i0,a,i0,a,a,a)', 'Requested nproc_dir [', &
        ctx%domain_cfg%nproc_dir(1), ',', ctx%domain_cfg%nproc_dir(2), ',', &
        ctx%domain_cfg%nproc_dir(3), ']: not supported (', trim(reason), ')'
      ctx%final_verdict = 'UNSUPPORTED'
    else
      call estimate_for_ng(ctx, requested_ng, requested_gib, exact, &
                           requested_verdict)
      ctx%final_verdict = requested_verdict
      print '(a,i0,a,f0.2,a,f0.1,a,a)', 'Requested nproc_dir gives ng=', &
        requested_ng, ': ', requested_gib, ' GiB/GPU (', &
        100._dp*requested_gib/ctx%card_gib, '% of card) - ', &
        trim(ctx%final_verdict)
      if (.not. exact) then
        if (trim(ctx%run_mode) == 'STATIC') then
          print '(a)', '  (Tier 1 estimate only: --static skips the GPU &
            &FFT plan query.)'
        else
          print '(a)', '  (Tier 1 floor only: this input is well over &
            &the limit, no GPU plan probe was attempted.)'
        end if
      end if
    end if
  end subroutine print_requested_line

  subroutine print_ng_table(ctx, smallest_fits)
    !! The estimate table, one row per scanned GPU count (the ng=1 row also
    !! feeds run_tier3's CHECK line); smallest_fits is 0 if none FITS.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(out) :: smallest_fits

    integer :: k, ng, local_dims(3)
    real(dp) :: per_gpu_gib, workspace_gib, overhead_gib, io_gib
    logical :: exact
    character(len=12) :: verdict
    character(len=96) :: reason

    call print_table_header()
    smallest_fits = 0
    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      reason = ng_unsupported_reason(ctx, ng)
      if (len_trim(reason) > 0) then
        print '(a,i0,a,a,a)', ' ', ng, '     (skipped: ', trim(reason), ')'
        cycle
      end if
      local_dims = [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)/ng]
      call estimate_for_ng(ctx, ng, per_gpu_gib, exact, verdict, &
                           workspace_gib, overhead_gib, io_gib)
      if (ng == 1) then
        ctx%estimate_ng1_gib = per_gpu_gib
        ctx%io_staging_ng1_gib = io_gib
      end if
      call print_table_row(ng, local_dims, workspace_gib, overhead_gib, &
                           per_gpu_gib, ctx%card_gib, verdict)
      if (smallest_fits == 0 .and. trim(verdict) == 'FITS') smallest_fits = ng
    end do
    call print_rule('-')
  end subroutine print_ng_table

  subroutine report(ctx)
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer :: smallest_fits

    call print_rule('=')
    print '(a,i0,a,i0,a,i0,a)', 'Input grid: ', ctx%gdims(1), 'x', &
      ctx%gdims(2), 'x', ctx%gdims(3)
    print '(a,f0.2,a,f0.2,a)', 'Card memory: ', ctx%card_gib, ' GiB (', &
      ctx%card_free_gib, ' GiB free now)'
    if (ctx%card_gib - ctx%card_free_gib > 0.5_dp) &
      print '(a,f0.2,a)', 'WARNING: ', ctx%card_gib - ctx%card_free_gib, &
        ' GiB of &
        &this card is in use by other processes; verdicts are against the &
        &full card, and the real build only proceeds if it fits in what is &
        &free.'
    print '(a,i0)', 'peak_fields (static estimate): ', ctx%peak_fields
    call print_io_staging_line(ctx)
    select case (trim(ctx%run_mode))
    case ('STATIC')
      print '(a)', 'Mode: static (Tier 1 only, no FFT plan query)'
    case ('BUILD')
      print '(a)', 'Mode: build (Tier 1 + Tier 2, then a real Tier 3 &
        &build/measure)'
    case default
      print '(a)', 'Mode: default (Tier 1 + Tier 2 FFT plan query)'
    end select
    call print_rule('=')
    call print_requested_line(ctx)

    call print_rule('-')
    call print_ng_table(ctx, smallest_fits)
    print '(a)', 'Notes: workspace + overhead = per-GPU. % is per-GPU &
      &against this card''s total memory.'

    if (smallest_fits > 0) then
      print '(a,i0,a)', 'Smallest GPU count that fits: ', smallest_fits, '.'
    else
      print '(a)', 'No scanned GPU count fits with headroom to spare.'
    end if

    if (trim(ctx%final_verdict) == 'BORDERLINE' .and. &
        trim(ctx%run_mode) /= 'BUILD') &
      print '(a)', 'To confirm with a real measurement, re-run with --build.'
  end subroutine report

end module m_memcheck_report
