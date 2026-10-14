module m_memcheck_report
  !! Printed report of x3d2-memcheck: header, tables and verdict lines.
  use m_common, only: dp
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

end module m_memcheck_report
