module m_memcheck_report
  !! Printed report of x3d2-memcheck: header, tables and verdict lines.
  use m_common, only: dp, i8, is_sp
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memory_estimate, only: gpu_io_staging_bytes
  use m_memcheck_estimate, only: to_gib, ng_unsupported_reason, &
                                 estimate_for_ng
  implicit none

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
                           ctx%final_verdict)
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

end module m_memcheck_report
