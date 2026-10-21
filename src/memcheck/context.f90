module m_memcheck_context
  !! Run state of x3d2-memcheck, passed explicitly to every routine.
  use m_common, only: dp
  use m_config, only: domain_config_t, solver_config_t, les_config_t, &
                      checkpoint_config_t
  use m_logger, only: abort_run
  use m_mesh, only: classify_boundaries, global_cell_dims, FFT_BC_000, &
                    FFT_BC_010, FFT_BC_100, FFT_BC_110
  use m_cuda_memcheck_device, only: cuda_memcheck_device_t
  implicit none

  private
  public :: memcheck_ctx_t, parse_args, read_config, classify_bc, &
            resolve_gpu_io_mode

  type :: memcheck_ctx_t
    type(cuda_memcheck_device_t) :: device
    character(len=256) :: input_path
    !> STATIC (--static), DEFAULT (no flag) or BUILD (--build).
    character(len=7) :: run_mode
    !> The substep count of --extensive, 0 when it is not given.
    integer :: extensive_substeps = 0
    type(domain_config_t) :: domain_cfg
    type(solver_config_t) :: solver_cfg
    type(les_config_t) :: les_cfg
    type(checkpoint_config_t) :: checkpoint_cfg
    integer :: gdims(3)
    integer :: cdims(3)
    logical :: periodic_x, periodic_y, periodic_z
    integer :: bc_code
    logical :: bc_is_000, bc_is_010, bc_is_100, bc_is_110
    logical :: multi_gpu_supported
    integer :: peak_fields
    real(dp) :: card_gib
    !> Free memory on the card at the start.
    real(dp) :: card_free_gib = 0._dp
    !> Whether snapshot and checkpoint writes stage through a device buffer,
    !> the write mode name, and the reason when there is no staging.
    logical :: gpu_io_device_write = .false.
    character(len=16) :: gpu_io_mode_name = 'auto'
    character(len=64) :: gpu_io_reason = ''
    !> The num_gpus=1 row of the estimate: the CHECK line of --build
    !> compares it with the measurement. The verdict is blank for an
    !> unsupported row.
    real(dp) :: estimate_ng1_gib
    character(len=12) :: ng1_verdict = ''
    !> The FFT probe plan: 'full' when the estimate uses the full-grid
    !> probe's result, else 'small' (also when the full-grid probe failed).
    !> probe_reason says why the full-grid probe has no result, and
    !> probe_note_shown is set once estimate_for_num_gpus has printed it.
    character(len=5) :: probe_plan = 'small'
    character(len=96) :: probe_reason = ''
    logical :: probe_note_shown = .false.
    !> The verdict of the requested nproc_dir, which sets the exit code.
    character(len=12) :: final_verdict
    !> The IO staging term (GiB) at num_gpus=1.
    real(dp) :: io_staging_ng1_gib = 0._dp
    !> The output fields and the substep count the real build measured.
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
    if (nargs < 1) call abort_run('usage: x3d2-memcheck <input.x3d> &
      &[--static | --build [--extensive <n_substeps>]]', collective=.true.)
    ! The input file stays argument 1: the real case constructors re-read it.
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
        if (i == nargs) call abort_run('x3d2-memcheck: --extensive needs a &
          &positive integer substep count, e.g. --extensive 1000', &
          collective=.true.)
        i = i + 1
        call get_command_argument(i, arg)
        read (arg, *, iostat=iostat_n) ctx%extensive_substeps
        if (iostat_n /= 0 .or. ctx%extensive_substeps < 1) &
          call abort_run('x3d2-memcheck: --extensive needs a positive &
            &integer substep count', collective=.true.)
        build_flag = .true.
      case default
        call abort_run('x3d2-memcheck: unknown flag '//trim(arg)// &
                       ' (expected --static, --build, or --extensive <n>)', &
                       collective=.true.)
      end select
      i = i + 1
    end do
    if (static_flag .and. build_flag) &
      call abort_run('x3d2-memcheck: --static and --build/--extensive are &
        &mutually exclusive (--static stops at the formula estimate).', &
        collective=.true.)

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
    !! The BC parse and four-way split of classify_boundaries in m_mesh.
    type(memcheck_ctx_t), intent(inout) :: ctx

    logical :: periodic(3)

    call classify_boundaries(ctx%domain_cfg%BC_x, ctx%domain_cfg%BC_y, &
                             ctx%domain_cfg%BC_z, periodic, ctx%bc_code)
    ctx%periodic_x = periodic(1)
    ctx%periodic_y = periodic(2)
    ctx%periodic_z = periodic(3)
    ctx%bc_is_000 = ctx%bc_code == FFT_BC_000
    ctx%bc_is_010 = ctx%bc_code == FFT_BC_010
    ctx%bc_is_100 = ctx%bc_code == FFT_BC_100
    ctx%bc_is_110 = ctx%bc_code == FFT_BC_110
    ctx%cdims = global_cell_dims(ctx%gdims, periodic)
  end subroutine classify_bc

  subroutine resolve_gpu_io_mode(ctx)
#ifdef X3D2_ADIOS2_CUDA
    use m_io_backend, only: gpu_write_mode_from_env, gpu_write_mode_names, &
                            gpu_write_mode_force_host
#endif
    !! Resolves whether snapshot and checkpoint writes stage through a
    !! device buffer, with the environment reader the solver uses for the
    !! write mode.
    type(memcheck_ctx_t), intent(inout) :: ctx
#ifdef X3D2_ADIOS2_CUDA
    integer :: mode

    mode = gpu_write_mode_from_env()

    ctx%gpu_io_mode_name = trim(gpu_write_mode_names(mode))
    if (mode == gpu_write_mode_force_host) then
      ctx%gpu_io_device_write = .false.
      ctx%gpu_io_reason = 'write mode host'
    else
      ctx%gpu_io_device_write = .true.
    end if
#else
    ctx%gpu_io_device_write = .false.
    ctx%gpu_io_mode_name = 'host'
    ctx%gpu_io_reason = 'build lacks X3D2_ADIOS2_CUDA'
#endif
  end subroutine resolve_gpu_io_mode

end module m_memcheck_context
