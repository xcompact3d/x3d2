module m_memcheck_context
  !! Run state of x3d2-memcheck, passed explicitly to every routine.
  implicit none

  type :: memcheck_ctx_t
    character(len=256) :: input_path
    !> STATIC = --static (Tier 1 only), DEFAULT = no flag (Tier 1+2, today's
    !> behaviour), BUILD = --build (adds Tier 3).
    character(len=7) :: run_mode
    !> Set by parse_args() from --extensive <n> (0 = not requested, the
    !> normal --build path of exactly 1 substep). build_and_measure() reads
    !> this to decide how many RK sub-stages to drive before measuring;
    !> report_measured_table()'s banner line and EXTENSIVE result line both
    !> read measured_n_substeps afterwards to state the real count.
    integer :: extensive_substeps = 0
  end type memcheck_ctx_t
end module m_memcheck_context
