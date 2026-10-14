module m_memcheck_context
  !! Run state of x3d2-memcheck, passed explicitly to every routine.
  use m_common, only: dp, i8
  use m_config, only: domain_config_t, solver_config_t, les_config_t, &
                      checkpoint_config_t
  use m_mesh, only: periodic_dir
  use m_memory_estimate, only: cell_dims
  use m_memcheck_device, only: memcheck_device_t
  implicit none

  type :: memcheck_ctx_t
    !> The compute backend, allocated by the program before anything else.
    class(memcheck_device_t), allocatable :: device
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
    !> Set by report()'s ng=1 row, used by run_tier3()'s CHECK line to
    !> compare the static estimate against the real --build measurement.
    real(dp) :: estimate_ng1_gib
    !> The requested configuration's verdict - set by report() from the
    !> estimate, and (--build only) overwritten by run_tier3() with the more
    !> authoritative measured verdict if a real build actually ran. Used
    !> after MPI_Finalize to set the process exit code: 0=FITS, 1=BORDERLINE,
    !> 2=DOES_NOT_FIT, so a calling script can check $? instead of scraping
    !> stdout.
    character(len=12) :: final_verdict
    !> GPU-aware IO staging term (GiB) at ng=1, captured by report()'s table
    !> loop and consumed by run_tier3()'s CHECK line, which must exclude it -
    !> the real --build never performs a snapshot/checkpoint write.
    real(dp) :: io_staging_ng1_gib = 0._dp
    !> orig_dir/scratch_dir: set by validate_build_inputs and make_scratch,
    !> used by leave_build_scratch to chdir back and remove the scratch
    !> directory - see their docstrings and run_tier3's docstring for why a
    !> real --build runs inside a throwaway subdirectory rather than the
    !> invoking directory.
    character(len=4096) :: orig_dir = '', scratch_dir = ''
    !> Queried once, at this process's own rank count (always 1 - see
    !> ensure_fft_query). These are real, live cudaMemGetInfo measurements
    !> for a SINGLE-RANK communicator; there is no way to measure the true
    !> ng>1 cost without actually running under mpirun -n ng, which this
    !> single-process estimator does not do. estimate_for_ng instead
    !> extrapolates ng>1 by the 1/ng law confirmed by a real mpirun -n 2
    !> xcompact run on examples/TGV/input.x3d (2026-09-10: 3339 MiB/GPU
    !> measured vs 3.21 GiB predicted this way) - ng1_worksize_bytes (plain
    !> cuFFT path only) and ng1_heap_bytes (cuFFTMp's NVSHMEM heap) both
    !> scale this way; ng1_xtdesc_bytes (cuFFTMp's distributed data buffer)
    !> does not - it is added only at ng=1, since it was measured to cost
    !> 0 extra at ng=2 (drawn from heap space the plan already reserved).
    integer(i8) :: ng1_worksize_bytes, ng1_heap_bytes, ng1_xtdesc_bytes
    logical :: ng1_used_cufftmp
    !> Set by ensure_fft_query, the cudaMemGetInfo delta across the whole
    !> throwaway plan query (create+destroy) - if cufftDestroy fully released
    !> the NVSHMEM heap this reserved, this should read ~0. Printed under
    !> --build only, right before the real build, since that is the only
    !> mode where a second cuFFTMp plan lifecycle follows in the same
    !> process and a non-zero residual would mean the real build's own
    !> measurement is inflated by whatever Tier 2 left behind.
    real(dp) :: ng1_query_residual_gib
    logical :: ng1_query_ran = .false.
    !> Whether a --build's real build actually used cuFFTMp (only known once
    !> the real solver is constructed) - the only real "cuFFTMp available"
    !> signal in this codebase, more authoritative than Tier 2's throwaway
    !> plan query (ng1_used_cufftmp above) for gating the measured table's
    !> ng>1 rows.
    logical :: use_cufftmp = .false., use_cufftmp_known = .false.
    !> Set inside build_and_measure, once the real config driven output
    !> gating is known - needed after it returns to cross-check the measured
    !> peak_fields against m_memory_estimate's static peak_fields_lookup (see
    !> check_peak_fields_table).
    logical :: output_vorticity = .false., output_qcriterion = .false.
    integer :: measured_n_substeps
    integer :: irank
  end type memcheck_ctx_t

contains

  subroutine parse_args(ctx)
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer :: i, nargs, iostat_n
    character(len=256) :: arg
    logical :: static_flag, build_flag

    nargs = command_argument_count()
    if (nargs < 1) error stop 'usage: x3d2-memcheck <input.x3d> &
      &[--static | --build [--extensive <n_substeps>]]'
    ! The input file must stay positional argument 1: the real case
    ! constructors this tool calls under --build (case_channel_init,
    ! case_cylinder_init, solver_init) independently re-read
    ! get_argument(1) themselves rather than taking the path as a
    ! parameter, so flags must only ever appear after it.
    call get_command_argument(1, ctx%input_path)

    static_flag = .false.
    build_flag = .false.
    i = 2
    do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--static')
        static_flag = .true.
      case ('--build')
        build_flag = .true.
      case ('--extensive')
        if (i == nargs) error stop 'x3d2-memcheck: --extensive needs a &
          &positive integer substep count, e.g. --extensive 1000'
        i = i + 1
        call get_command_argument(i, arg)
        read (arg, *, iostat=iostat_n) ctx%extensive_substeps
        if (iostat_n /= 0 .or. ctx%extensive_substeps < 1) &
          error stop 'x3d2-memcheck: --extensive needs a positive &
            &integer substep count'
        build_flag = .true.
      case default
        error stop 'x3d2-memcheck: unknown flag '//trim(arg)// &
          ' (expected --static, --build, or --extensive <n>)'
      end select
      i = i + 1
    end do
    if (static_flag .and. build_flag) &
      error stop 'x3d2-memcheck: --static and --build/--extensive are &
        &mutually exclusive (--static stops at Tier 1).'

    if (static_flag) then
      ctx%run_mode = 'STATIC'
    else if (build_flag) then
      ctx%run_mode = 'BUILD'
    else
      ctx%run_mode = 'DEFAULT'
    end if
  end subroutine parse_args

  subroutine read_config(ctx)
    type(memcheck_ctx_t), intent(inout) :: ctx

    call ctx%domain_cfg%read(nml_file=trim(ctx%input_path))
    call ctx%solver_cfg%read(nml_file=trim(ctx%input_path))
    call ctx%les_cfg%read(nml_file=trim(ctx%input_path))
    call ctx%checkpoint_cfg%read(nml_file=trim(ctx%input_path))
    ctx%gdims = ctx%domain_cfg%dims_global
  end subroutine read_config

  subroutine classify_bc(ctx)
    !! Mirrors src/backend/cuda/poisson_fft.f90:227-239's four-way split.
    type(memcheck_ctx_t), intent(inout) :: ctx

    ctx%periodic_x = periodic_dir(ctx%domain_cfg%BC_x)
    ctx%periodic_y = periodic_dir(ctx%domain_cfg%BC_y)
    ctx%periodic_z = periodic_dir(ctx%domain_cfg%BC_z)
    ctx%bc_is_000 = ctx%periodic_x .and. ctx%periodic_y .and. ctx%periodic_z
    ctx%bc_is_010 = ctx%periodic_x .and. (.not. ctx%periodic_y) .and. &
                    ctx%periodic_z
    ctx%bc_is_100 = (.not. ctx%periodic_x) .and. ctx%periodic_y .and. &
                    ctx%periodic_z
    ctx%bc_is_110 = (.not. ctx%periodic_x) .and. &
                    (.not. ctx%periodic_y) .and. ctx%periodic_z
    ctx%cdims = cell_dims(ctx%gdims, [ctx%periodic_x, ctx%periodic_y, &
                                      ctx%periodic_z])
  end subroutine classify_bc

#ifdef X3D2_ADIOS2_CUDA
  pure function to_lower(text) result(lowered)
    !! Small private ASCII lower-caser for resolve_gpu_io_mode - not shared
    !! with src/io/adios2/io.f90's to_lower_ascii, which is private to that
    !! module (see resolve_gpu_io_mode's comment).
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i, code

    lowered = text
    do i = 1, len(text)
      code = iachar(lowered(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lowered(i:i) = achar(code + 32)
      end if
    end do
  end function to_lower
#endif

  subroutine resolve_gpu_io_mode(ctx)
    !! Resolve whether this run's snapshot/checkpoint writes would stage
    !! through a device buffer. Mirrors init_runtime_options in
    !! src/io/adios2/io.f90:296-313, which is private to the ADIOS2 writer
    !! and belongs to PR 277, so it is not shared here - keep the two in
    !! sync by hand if the env var contract changes.
    type(memcheck_ctx_t), intent(inout) :: ctx
#ifdef X3D2_ADIOS2_CUDA
    character(len=64) :: raw_value, mode_value
    integer :: status, value_length

    call get_environment_variable('X3D2_ADIOS2_GPU_WRITE_MODE', raw_value, &
                                  length=value_length, status=status)
    if (status /= 0 .or. value_length == 0) then
      mode_value = 'auto'
    else
      mode_value = to_lower(adjustl(raw_value(1:min(value_length, &
                                                    len(raw_value)))))
    end if

    select case (trim(mode_value))
    case ('host', 'd2h', 'staged')
      ctx%gpu_io_device_write = .false.
      ctx%gpu_io_mode_name = 'host'
      ctx%gpu_io_reason = 'write mode host'
    case default
      ctx%gpu_io_device_write = .true.
      ctx%gpu_io_mode_name = trim(mode_value)
    end select
#else
    ctx%gpu_io_device_write = .false.
    ctx%gpu_io_mode_name = 'host'
    ctx%gpu_io_reason = 'build lacks X3D2_ADIOS2_CUDA'
#endif
  end subroutine resolve_gpu_io_mode

end module m_memcheck_context
