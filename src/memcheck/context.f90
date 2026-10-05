module m_memcheck_context
  !! Run state of x3d2-memcheck, passed explicitly to every routine.
  use m_common, only: dp
  use m_config, only: domain_config_t, solver_config_t, les_config_t, &
                      checkpoint_config_t
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
    type(domain_config_t) :: domain_cfg
    type(solver_config_t) :: solver_cfg
    type(les_config_t) :: les_cfg
    type(checkpoint_config_t) :: checkpoint_cfg
    integer :: gdims(3)
    integer :: cdims(3)
    logical :: periodic_x, periodic_y, periodic_z
    logical :: bc_is_000, bc_is_010, bc_is_100, bc_is_110
    logical :: multi_gpu_supported
    integer :: peak_fields
    real(dp) :: card_gib
    !> Free memory on the card right now (query_card_gib), i.e. card_gib minus
    !> whatever other processes already hold - used by report()'s header
    !> WARNING line and by run_tier3's BUILD_FRACTION gate, both of which must
    !> judge headroom against what is actually free, not the card's full
    !> capacity.
    real(dp) :: card_free_gib = 0._dp
    !> Set by resolve_gpu_io_mode() - whether this run's snapshot/checkpoint
    !> writes would stage through a device buffer (true) or fall back to a
    !> host copy (false) before reaching disk.
    logical :: gpu_io_device_write = .false.
    !> Resolved write mode name (mirrors runtime_gpu_write_mode_name in
    !> src/io/adios2/io.f90): 'auto', 'gpu', or 'host'.
    character(len=16) :: gpu_io_mode_name = 'auto'
    !> Human-readable reason gpu_io_device_write is false, used by report()'s
    !> GPU-aware IO staging line when there is no term to report.
    character(len=64) :: gpu_io_reason = ''
  end type memcheck_ctx_t
end module m_memcheck_context
