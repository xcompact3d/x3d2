program x3d2_memcheck
  !! Standalone GPU memory footprint estimator AND real build probe for an
  !! x3d2 .x3d input file, in one tool:
  !!   (no flag)  Tier 1 (pure formula, no GPU build) + Tier 2 (a throwaway
  !!              cuFFT/cuFFTMp plan query, needs a GPU but builds nothing
  !!              else) unless Tier 1 already answers DOES_NOT_FIT.
  !!   --static   Tier 1 only - never touches cuFFT, instant, deterministic.
  !!   --build    the above, then Tier 3 - a REAL case build + one substep
  !!              on 1 GPU, measured via cudaMemGetInfo, with a
  !!              self verifying drift check against the static
  !!              peak_fields_lookup table and a CHECK line comparing the
  !!              static estimate to the real measurement. Runs in the SAME
  !!              process (unlike the two binary design this replaced): no
  !!              nested mpirun limitation (mpirun refuses a second mpirun
  !!              from a process it already launched, which used to mean a
  !!              separate test_memory_estimate_cuda_1 binary had to be run
  !!              by hand to confirm a BORDERLINE result).
  !!              The real build itself runs inside a throwaway
  !!              x3d2-memcheck-build.<pid> subdirectory of the invoking
  !!              directory, because a case build initialises monitoring
  !!              (writes monitoring.csv) and calls postprocess(0), which
  !!              clobbered run directories on 2026-09-15; a relative
  !!              input path is mirrored in by symlink, one containing
  !!              '..' must be passed as absolute instead; the scratch
  !!              directory is removed once the build finishes.
  !!
  !! Also estimates the GPU-aware ADIOS2 I/O staging buffer (one extra
  !! unpadded local field held on the device while a snapshot/checkpoint
  !! field is packed - src/io/adios2/io.f90) whenever a device-side write is
  !! in play, added analytically via m_memory_estimate's
  !! gpu_io_staging_bytes rather than measured, since --build's real case
  !! never performs a snapshot/checkpoint write. Which mode is in play
  !! follows the X3D2_ADIOS2_GPU_WRITE_MODE environment variable, resolved
  !! the same way src/io/adios2/io.f90's own runtime option does (see
  !! resolve_gpu_io_mode below).
  !!
  !! Usage: x3d2-memcheck <input.x3d> [--static | --build]
  use mpi
  use iso_c_binding, only: c_char, c_int, c_size_t, c_ptr, c_null_char
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice
  use m_common, only: dp, i8, nbytes, VERT, is_sp
  use m_config, only: domain_config_t, solver_config_t, les_config_t, &
                      checkpoint_config_t
  use m_mesh, only: mesh_t, periodic_dir
  use m_cuda_common, only: SZ
  use m_memory_estimate, only: padded_cells, cell_dims, &
                               spectral_slab_bytes, mirror_buffer_bytes_100, &
                               spectral_extra_bytes_110, &
                               stretched_y_matrix_bytes, output_field_active, &
                               peak_fields_lookup, padded_halo_bytes, &
                               gpu_io_staging_bytes
  use m_cuda_memory_estimate, only: fft_workspace_bytes_query, &
                                    context_floor_bytes, check_status
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

  implicit none

  !> <80% of card memory: FITS. 80-95%: BORDERLINE. >95%: DOES_NOT_FIT.
  real(dp), parameter :: FITS_FRACTION = 0.80_dp
  real(dp), parameter :: DOES_NOT_FIT_FRACTION = 0.95_dp
  !> --build: do not attempt a real build whose fields only workspace alone
  !> exceeds this fraction of the FREE card memory (card_free_gib, not the
  !> card's full capacity) - leaves headroom for context + FFT scratch so
  !> the build never OOMs. Same value and role this tool's real build probe
  !> has always used.
  real(dp), parameter :: BUILD_FRACTION = 0.55_dp
  !> --build's CHECK line: how far the static estimate may be from the real
  !> measurement before it is flagged MISMATCH rather than OK. Justified by
  !> this tool's own calibration record: TGV matched to 1.5%; a real
  !> calibration bug (an LES peak_fields delta wrong by 8) showed up as an
  !> ~8% byte mismatch, well outside this tolerance - and is separately,
  !> exactly caught by the drift check below regardless of byte tolerance.
  real(dp), parameter :: CHECK_TOLERANCE = 0.05_dp
  integer, parameter :: n_halo = 4
  !> Candidate GPU counts to scan for "smallest ng that fits". The CUDA
  !> backend only supports a Z-only pencil decomposition (nproc_dir=
  !> [1,1,ng]) today.
  integer, parameter :: n_gpu_list(4) = [1, 2, 4, 8]

  !> What a --build measurement fixes for the per-ng table rows: the
  !> measured peak_fields, the measured per-GPU overhead at ng=1 (used GiB
  !> minus the exact workspace), and whether ng>1 is usable given the FFT
  !> path the real build took.
  type :: measured_basis_t
    integer :: peak_fields
    real(dp) :: overhead_gib
    logical :: multi_gpu_ok
  end type measured_basis_t

  character(len=256) :: input_path
  !> STATIC = --static (Tier 1 only), DEFAULT = no flag (Tier 1+2, today's
  !> behaviour), BUILD = --build (adds Tier 3).
  character(len=7) :: run_mode
  integer :: ierr, irank, nproc, ndevs
  type(domain_config_t) :: domain_cfg
  type(solver_config_t) :: solver_cfg
  type(les_config_t) :: les_cfg
  type(checkpoint_config_t) :: checkpoint_cfg
  integer :: gdims(3), cdims(3)
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
  !> Set by parse_args() from --extensive <n> (0 = not requested, the
  !> normal --build path of exactly 1 substep). build_and_measure() reads
  !> this to decide how many RK sub-stages to drive before measuring;
  !> report_measured_table()'s banner line and EXTENSIVE result line both
  !> read measured_n_substeps afterwards to state the real count.
  integer :: extensive_substeps = 0
  integer :: measured_n_substeps
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

  interface
    function c_chdir(path) bind(C, name='chdir') result(rc)
      import :: c_char, c_int
      character(kind=c_char), dimension(*), intent(in) :: path
      integer(c_int) :: rc
    end function c_chdir
    function c_getcwd(buf, size) bind(C, name='getcwd') result(ptr)
      import :: c_char, c_size_t, c_ptr
      character(kind=c_char), dimension(*), intent(inout) :: buf
      integer(c_size_t), value :: size
      type(c_ptr) :: ptr
    end function c_getcwd
    function c_getpid() bind(C, name='getpid') result(pid)
      import :: c_int
      integer(c_int) :: pid
    end function c_getpid
  end interface

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 1) error stop 'x3d2-memcheck: run single-rank (mpirun -n 1).'

  ierr = cudaGetDeviceCount(ndevs)
  if (ierr /= 0 .or. ndevs < 1) then
    print '(a,i0,a,i0,a)', 'x3d2-memcheck: no usable CUDA device &
      &(cudaGetDeviceCount code ', ierr, ', devices ', ndevs, ')'
    call no_device_exit()
  end if
  ierr = cudaSetDevice(0)
  if (ierr /= 0) then
    print '(a,i0,a)', 'x3d2-memcheck: no usable CUDA device (cudaSetDevice &
      &failed, code ', ierr, ')'
    call no_device_exit()
  end if

  call parse_args()
  call read_config()
  call resolve_gpu_io_mode()
  call classify_bc()
  call query_card_gib()

  peak_fields = peak_fields_lookup(trim(domain_cfg%flow_case_name), &
                                   solver_cfg%n_species, &
                                   trim(les_cfg%model) /= 'none', &
                                   solver_cfg%ibm_on, &
                                   output_field_active(checkpoint_cfg, &
                                                       'vorticity'), &
                                   output_field_active(checkpoint_cfg, &
                                                       'qcriterion'), &
                                   solver_cfg%lowmem_transeq)

  multi_gpu_supported = .not. (bc_is_010 .or. bc_is_110)

  call report()
  if (trim(run_mode) == 'BUILD') call run_tier3()

  call MPI_Finalize(ierr)

  ! Exit code mirrors the final verdict: 0=FITS, 1=BORDERLINE,
  ! 2=DOES_NOT_FIT (also covers UNSUPPORTED - an unsupported nproc_dir
  ! request, via the default case below, since neither fits nor merely
  ! borderline is the right signal for it).
  select case (trim(final_verdict))
  case ('FITS')
    call exit(0)
  case ('BORDERLINE')
    call exit(1)
  case default
    call exit(2)
  end select

contains

  subroutine no_device_exit()
    !! Common tail of both no-device exits above (which print their own
    !! reason first): shut MPI down cleanly, then error stop.
    integer :: ierr_finalize

    call MPI_Finalize(ierr_finalize)
    error stop 'x3d2-memcheck: no usable CUDA device'
  end subroutine no_device_exit

  subroutine parse_args()
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
    call get_command_argument(1, input_path)

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
        read (arg, *, iostat=iostat_n) extensive_substeps
        if (iostat_n /= 0 .or. extensive_substeps < 1) &
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
      run_mode = 'STATIC'
    else if (build_flag) then
      run_mode = 'BUILD'
    else
      run_mode = 'DEFAULT'
    end if
  end subroutine parse_args

  subroutine read_config()
    call domain_cfg%read(nml_file=trim(input_path))
    call solver_cfg%read(nml_file=trim(input_path))
    call les_cfg%read(nml_file=trim(input_path))
    call checkpoint_cfg%read(nml_file=trim(input_path))
    gdims = domain_cfg%dims_global
  end subroutine read_config

  subroutine resolve_gpu_io_mode()
    !! Resolve whether this run's snapshot/checkpoint writes would stage
    !! through a device buffer. Mirrors init_runtime_options in
    !! src/io/adios2/io.f90:296-313, which is private to the ADIOS2 writer
    !! and belongs to PR 277, so it is not shared here - keep the two in
    !! sync by hand if the env var contract changes.
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
      gpu_io_device_write = .false.
      gpu_io_mode_name = 'host'
      gpu_io_reason = 'write mode host'
    case default
      gpu_io_device_write = .true.
      gpu_io_mode_name = trim(mode_value)
    end select
#else
    gpu_io_device_write = .false.
    gpu_io_mode_name = 'host'
    gpu_io_reason = 'build lacks X3D2_ADIOS2_CUDA'
#endif
  end subroutine resolve_gpu_io_mode

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

  subroutine classify_bc()
    !! Mirrors src/backend/cuda/poisson_fft.f90:227-239's four-way split.
    periodic_x = periodic_dir(domain_cfg%BC_x)
    periodic_y = periodic_dir(domain_cfg%BC_y)
    periodic_z = periodic_dir(domain_cfg%BC_z)
    bc_is_000 = periodic_x .and. periodic_y .and. periodic_z
    bc_is_010 = periodic_x .and. (.not. periodic_y) .and. periodic_z
    bc_is_100 = (.not. periodic_x) .and. periodic_y .and. periodic_z
    bc_is_110 = (.not. periodic_x) .and. (.not. periodic_y) .and. periodic_z
    cdims = cell_dims(gdims, [periodic_x, periodic_y, periodic_z])
  end subroutine classify_bc

  function current_dir() result(path)
    !! Working directory via libc getcwd, trimmed at the first embedded NUL
    !! (the C string terminator - Fortran character variables are not
    !! themselves NUL-terminated, so the raw buffer must be cut there before
    !! use). Used by validate_build_inputs to record the invoking directory.
    character(len=4096) :: path
    character(kind=c_char, len=4096) :: buf
    type(c_ptr) :: ptr
    integer :: nul_pos

    ptr = c_getcwd(buf, 4096_c_size_t)
    nul_pos = index(buf, c_null_char)
    if (nul_pos > 0) then
      path = buf(1:nul_pos - 1)
    else
      path = buf
    end if
  end function current_dir

  pure function ibm_mask_filename(px, py, pz) result(fname)
    !! Build the ibm_<BC-suffix>.bp mask filename from the three periodic_BC
    !! flags, e.g. mesh%grid%periodic_BC(1:3) in build_and_measure or the
    !! module periodic_x/y/z classify_bc sets - both carry the same
    !! information since both derive from domain_cfg%BC_x/y/z via
    !! periodic_dir. Shared by build_and_measure's ibm_on pre-check and
    !! make_scratch's mask mirroring so the suffix logic can never
    !! diverge between the two (src/module/ibm.f90:70-75 is the third,
    !! authoritative copy this mirrors).
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

  pure function has_dotdot_component(path) result(found)
    !! True if any '/'-separated component of path is exactly '..' - used
    !! by make_scratch to reject a relative input path it cannot
    !! safely mirror by symlink.
    character(len=*), intent(in) :: path
    logical :: found
    integer :: start, slash_pos

    found = .false.
    start = 1
    do
      slash_pos = index(path(start:), '/')
      if (slash_pos == 0) then
        if (path(start:) == '..') found = .true.
        exit
      end if
      if (path(start:start + slash_pos - 2) == '..') found = .true.
      start = start + slash_pos
    end do
  end function has_dotdot_component

  subroutine query_card_gib()
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo (card total)')
    card_gib = to_gib(int(total_b, i8))
    card_free_gib = to_gib(int(free_b, i8))
  end subroutine query_card_gib

  pure function classify(per_gpu_gib, total_gib) result(verdict)
    !! Shared three-way verdict classification, used by both the static
    !! estimate (estimate_for_ng) and the real --build measurement
    !! (report_measured_table) so the two paths can never silently diverge
    !! on what counts as FITS/BORDERLINE/DOES_NOT_FIT.
    real(dp), intent(in) :: per_gpu_gib, total_gib
    character(len=12) :: verdict

    if (per_gpu_gib < FITS_FRACTION*total_gib) then
      verdict = 'FITS'
    else if (per_gpu_gib <= DOES_NOT_FIT_FRACTION*total_gib) then
      verdict = 'BORDERLINE'
    else
      verdict = 'DOES_NOT_FIT'
    end if
  end function classify

  pure function to_gib(n_bytes) result(gib)
    !! Bytes to GiB, the one conversion every figure this tool prints goes
    !! through, so the expression (and so the printed digits) cannot drift
    !! between call sites.
    integer(i8), intent(in) :: n_bytes
    real(dp) :: gib

    gib = real(n_bytes, dp)/1024._dp**3
  end function to_gib

  function ng_unsupported_reason(ng) result(reason)
    !! Whether GPU count ng is a decomposition this tool/solver actually
    !! supports - blank when it is. Shared by report()'s per-ng table,
    !! report_measured_table()'s per-ng table, and report()'s requested
    !! nproc_dir line, so the three can never silently diverge on what
    !! counts as unsupported. Checked in this order: ng itself must be a
    !! positive count; ng>1 needs BC_y periodic (multi_gpu_supported);
    !! gdims(3) must divide by ng (the physical Z decomposition); and, for
    !! whichever BC actually decomposes the spectral slab along dim2 (000/
    !! 010: cy: 100: cx - see spectral_slab_bytes), that dimension must
    !! also divide by ng, or the solver would truncate or error-stop.
    integer, intent(in) :: ng
    character(len=96) :: reason

    reason = ''
    if (ng < 1) then
      reason = 'nproc_dir(3) must be >= 1'
    else if (ng > 1 .and. .not. multi_gpu_supported) then
      reason = 'not supported: BC_y is non-periodic - &
               &src/poisson_fft.f90 error-stops at nproc>1'
    else if (mod(gdims(3), ng) /= 0) then
      write (reason, '(a,i0,a)') 'z=', gdims(3), ' not divisible'
    else if ((bc_is_000 .or. bc_is_010) .and. mod(cdims(2), ng) /= 0) then
      write (reason, '(a,i0,a)') 'y cells=', cdims(2), ' not divisible &
        &by ng, the solver would truncate the spectral slab'
    else if (bc_is_100 .and. mod(cdims(1), ng) /= 0) then
      write (reason, '(a,i0,a)') 'x cells=', cdims(1), ' not divisible &
        &by ng, src/backend/cuda/poisson_fft.f90 error-stops'
    end if
  end function ng_unsupported_reason

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

  function fields_plus_halo_bytes_n(ng, npeak) result(nbytes8)
    !! Same formula as fields_plus_halo_bytes, but for an explicit npeak
    !! rather than the module's static peak_fields - used by the --build
    !! measured table, which has its own real measured field count.
    integer, intent(in) :: ng, npeak
    integer(i8) :: nbytes8
    integer :: local_dims(3)

    local_dims = [gdims(1), gdims(2), gdims(3)/ng]
    nbytes8 = int(npeak, i8)*padded_cells(local_dims, SZ) &
              *int(nbytes, i8) + padded_halo_bytes(local_dims, SZ, n_halo)
  end function fields_plus_halo_bytes_n

  function fields_plus_halo_bytes(ng) result(nbytes8)
    !! Exact, ng dependent field+halo term (Tier 1, no GPU/plan needed),
    !! using the static peak_fields_lookup value.
    integer, intent(in) :: ng
    integer(i8) :: nbytes8

    nbytes8 = fields_plus_halo_bytes_n(ng, peak_fields)
  end function fields_plus_halo_bytes

  function spectral_plus_mirror_bytes(ng) result(nbytes8)
    !! waves_dev (all cases, src/backend/cuda/poisson_fft.f90:322-324) is
    !! exactly spectral_slab_bytes. Two BC-specific extra terms:
    !!   100, ng>1: mirror/exchange buffers (:339-347), zero at ng<=1.
    !!   110 (never ng>1 - see multi_gpu_supported): a SECOND complex array
    !!     of the same spectral shape, c_dev (:416-418), plus a real,
    !!     globally-shaped transposed workspace r_dev_110(nz_glob,nx_glob,
    !!     ny_glob) (:414) - i.e. the same total element count as cdims,
    !!     just permuted. Missing this term underestimated the 110 case by
    !!     ~14% (0.98 vs 1.134 GiB measured) during this feature's design.
    !!   010, non-uniform y-stretching (never ng>1 - 010 error-stops at
    !!     nproc>1): the a_*_dev Poisson coefficient matrices, see
    !!     stretched_y_matrix_bytes.
    integer, intent(in) :: ng
    integer(i8) :: nbytes8

    nbytes8 = spectral_slab_bytes(bc_is_100, bc_is_110, cdims, ng)
    if (bc_is_100) then
      nbytes8 = nbytes8 + mirror_buffer_bytes_100(cdims, ng)
    else if (bc_is_110) then
      nbytes8 = nbytes8 + spectral_extra_bytes_110(cdims, ng)
    else if (bc_is_010) then
      nbytes8 = nbytes8 + stretched_y_matrix_bytes(bc_is_010, &
                          domain_cfg%stretching(2), solver_cfg%lowmem_fft, &
                          cdims, ng)
    end if
  end function spectral_plus_mirror_bytes

  subroutine ensure_fft_query()
    !! Runs the real (single-rank) Tier 2 plan query at most once per
    !! process and caches the result in ng1_worksize_bytes/ng1_heap_bytes/
    !! ng1_xtdesc_bytes/ng1_used_cufftmp, plus the cudaMemGetInfo delta
    !! across the whole call in ng1_query_residual_gib (see its
    !! declaration).
    logical, save :: done = .false.
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b

    if (done) return
    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before FFT probe)')
    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   irank == 0, ng1_worksize_bytes, &
                                   ng1_heap_bytes, ng1_xtdesc_bytes, &
                                   ng1_used_cufftmp)
    ierr = cudaMemGetInfo(free_after, total_b)
    call check_status(ierr, 'cudaMemGetInfo (after FFT probe)')
    ng1_query_residual_gib = to_gib(int(free_before, i8) - &
                                    int(free_after, i8))
    ng1_query_ran = .true.
    done = .true.
  end subroutine ensure_fft_query

  subroutine estimate_for_ng(ng, per_gpu_gib, exact, verdict, workspace_gib, &
                             overhead_gib, io_gib)
    !! Tier 1 floor first; only calls into Tier 2 (a real, throwaway GPU
    !! FFT plan) when Tier 1 alone cannot already answer DOES_NOT_FIT, and
    !! never under --static. workspace_gib/overhead_gib (optional) are the
    !! same two-term breakdown the --build measured table reports:
    !! workspace = exact fields+halo+spectral+mirror (Tier 1, no GPU
    !! needed); overhead = everything from the Tier 2 GPU query, the
    !! hardware constants (worksize/xtdesc/heap/context), and the
    !! GPU-aware IO staging term (also broken out on its own via the
    !! optional io_gib, in GiB, for callers that need to track it apart
    !! from context/FFT overhead - e.g. run_tier3()'s CHECK line, which
    !! must exclude it since the real --build never performs a
    !! snapshot/checkpoint write).
    integer, intent(in) :: ng
    real(dp), intent(out) :: per_gpu_gib
    logical, intent(out) :: exact
    character(len=12), intent(out) :: verdict
    real(dp), intent(out), optional :: workspace_gib, overhead_gib, io_gib

    integer(i8) :: base_bytes, spec_bytes, worksize_bytes, floor_bytes, &
                  xtdesc_bytes, heap_bytes, context_bytes, io_bytes
    logical :: used_cufftmp
    real(dp) :: floor_gib, io_gib_local

    base_bytes = fields_plus_halo_bytes(ng)
    spec_bytes = spectral_plus_mirror_bytes(ng)
    io_bytes = gpu_io_staging_bytes([gdims(1), gdims(2), gdims(3)/ng], &
                                    checkpoint_cfg, gpu_io_device_write)
    io_gib_local = to_gib(io_bytes)
    if (present(io_gib)) io_gib = io_gib_local
    if (present(workspace_gib)) &
      workspace_gib = to_gib(base_bytes + spec_bytes)
    ! Floor: assume cuFFTMp is used (the realistic case for every BC this
    ! solver ever attempts it for - see m_cuda_memory_estimate's 110 guard)
    ! since that gives the larger, safer lower bound.
    floor_bytes = base_bytes + spec_bytes + io_bytes + &
                  context_floor_bytes(ng, .not. bc_is_110)
    floor_gib = to_gib(floor_bytes)

    if (trim(run_mode) == 'STATIC' .or. &
        floor_gib > DOES_NOT_FIT_FRACTION*card_gib) then
      per_gpu_gib = floor_gib
      exact = .false.
      if (present(overhead_gib)) &
        overhead_gib = to_gib(context_floor_bytes(ng, &
                                                  .not. bc_is_110)) &
                       + io_gib_local
      if (trim(run_mode) == 'STATIC') then
        ! --static: Tier 1 only, never queries the GPU FFT plan. floor_gib
        ! IS the estimate here - not just a DOES_NOT_FIT early-return floor
        ! - so classify it with the full three-way verdict.
        verdict = classify(per_gpu_gib, card_gib)
      else
        verdict = 'DOES_NOT_FIT'
      end if
      return
    end if

    call ensure_fft_query()
    used_cufftmp = ng1_used_cufftmp
    ! context_bytes: the CUDA-context-alone baseline (no ng argument
    ! matters here - context_floor_bytes only adds its hardcoded NVSHMEM
    ! heap when uses_cufftmp=.true., so passing .false. always yields just
    ! the context). The cuFFTMp heap itself now comes from the live
    ! ng1_heap_bytes measurement below instead of that hardcoded constant.
    context_bytes = context_floor_bytes(ng, .false.)
    if (used_cufftmp) then
      ! Live-measured (fft_workspace_bytes_query): worksize is carved out
      ! of the NVSHMEM heap for cuFFTMp, not a separate allocation - do
      ! not add it on top of heap_bytes (see that function's docstring).
      worksize_bytes = 0_i8
      heap_bytes = int(real(ng1_heap_bytes, dp)/real(ng, dp), i8)
      ! xtdesc cost was measured to be 0 extra at ng=2 (drawn from heap
      ! space the plan already reserved) - add the ng=1 value only there.
      xtdesc_bytes = merge(ng1_xtdesc_bytes, 0_i8, ng == 1)
    else
      ! See ng1_worksize_bytes' declaration: this is an extrapolation for
      ! ng>1, not an independent per-ng measurement (though in practice
      ! this branch only ever runs at ng=1 - plain cuFFT is only reached
      ! by the 110 BC or a cuFFTMp fallback, neither of which supports
      ! nproc>1 in this solver).
      worksize_bytes = int(real(ng1_worksize_bytes, dp)/real(ng, dp), i8)
      heap_bytes = 0_i8
      xtdesc_bytes = 0_i8
    end if
    per_gpu_gib = to_gib(base_bytes + spec_bytes + worksize_bytes + &
                         xtdesc_bytes + heap_bytes + context_bytes + &
                         io_bytes)
    exact = .true.
    if (present(overhead_gib)) &
      overhead_gib = to_gib(worksize_bytes + xtdesc_bytes + heap_bytes + &
                            context_bytes) + io_gib_local

    verdict = classify(per_gpu_gib, card_gib)
  end subroutine estimate_for_ng

  subroutine print_io_staging_line()
    !! The GPU-aware IO staging line of the report header: the term at ng=1
    !! with the snapshot/checkpoint kind that causes it, or why there is
    !! none.
    real(dp) :: mib
    logical :: unit_stride, snapshot_active, checkpoint_active
    character(len=64) :: what, io_none_reason
    integer(i8) :: io_ng1_bytes

    ! GPU-aware IO staging: computed directly at ng=1 (local_dims == gdims
    ! there) so it is available here, ahead of the per-ng table below.
    io_ng1_bytes = gpu_io_staging_bytes(gdims, checkpoint_cfg, &
                                        gpu_io_device_write)
    unit_stride = all(checkpoint_cfg%output_stride == 1)
    snapshot_active = checkpoint_cfg%snapshot_freq > 0 .and. unit_stride
    checkpoint_active = checkpoint_cfg%checkpoint_freq > 0
    if (io_ng1_bytes > 0_i8) then
      mib = to_gib(io_ng1_bytes)*1024._dp
      if (snapshot_active .and. checkpoint_active) then
        what = 'snapshot at unit stride + checkpoint'
      else if (checkpoint_active) then
        what = 'checkpoint'
      else
        what = 'snapshot at unit stride'
      end if
      if (.not. checkpoint_active .and. checkpoint_cfg%snapshot_sp .and. &
          .not. is_sp .and. snapshot_active) what = trim(what)//' (sp)'
      print '(a,f0.1,a,a,a,a,a)', 'GPU-aware IO staging: ', mib, &
        ' MiB/GPU at ng=1 (write mode ', trim(gpu_io_mode_name), ', ', &
        trim(what), ')'
    else
      if (.not. gpu_io_device_write) then
        io_none_reason = trim(gpu_io_reason)
      else if (checkpoint_cfg%snapshot_freq > 0 .and. .not. unit_stride &
              .and. checkpoint_cfg%checkpoint_freq == 0) then
        io_none_reason = 'snapshot striding falls back to host path'
      else
        io_none_reason = 'no unit-stride snapshot and no checkpoint enabled'
      end if
      print '(a,a,a)', 'GPU-aware IO staging: none (', trim(io_none_reason), &
        ')'
    end if
  end subroutine print_io_staging_line

  subroutine print_requested_line()
    !! The verdict line for the nproc_dir the input file asks for (or why
    !! that decomposition is not supported); sets final_verdict.
    integer :: requested_ng
    real(dp) :: requested_gib
    logical :: exact
    character(len=96) :: reason

    requested_ng = domain_cfg%nproc_dir(3)
    reason = ng_unsupported_reason(requested_ng)
    if (domain_cfg%nproc_dir(1) /= 1 .or. domain_cfg%nproc_dir(2) /= 1) &
      reason = 'only nproc_dir = [1,1,ng] is supported by the CUDA backend'
    if (len_trim(reason) > 0) then
      print '(a,i0,a,i0,a,i0,a,a,a)', 'Requested nproc_dir [', &
        domain_cfg%nproc_dir(1), ',', domain_cfg%nproc_dir(2), ',', &
        domain_cfg%nproc_dir(3), ']: not supported (', trim(reason), ')'
      final_verdict = 'UNSUPPORTED'
    else
      call estimate_for_ng(requested_ng, requested_gib, exact, final_verdict)
      print '(a,i0,a,f0.2,a,f0.1,a,a)', 'Requested nproc_dir gives ng=', &
        requested_ng, ': ', requested_gib, ' GiB/GPU (', &
        100._dp*requested_gib/card_gib, '% of card) - ', trim(final_verdict)
      if (.not. exact) then
        if (trim(run_mode) == 'STATIC') then
          print '(a)', '  (Tier 1 estimate only: --static skips the GPU &
            &FFT plan query.)'
        else
          print '(a)', '  (Tier 1 floor only: this input is well over &
            &the limit, no GPU plan probe was attempted.)'
        end if
      end if
    end if
  end subroutine print_requested_line

  subroutine print_ng_table(smallest_fits, basis, requested_verdict)
    !! Table header and one row per scanned GPU count (an unsupported count
    !! prints why it is skipped), then a closing rule. Shared by report()
    !! and report_measured_table(): rows come from estimate_for_ng, or from
    !! measured_row when basis (the --build measurement) is given.
    !! smallest_fits (optional) is the first scanned count that FITS, 0 if
    !! none. requested_verdict (optional) is set to the verdict of the
    !! requested nproc_dir(3) row when that row is printed. The estimate
    !! ng=1 row also records the estimate and its IO staging term for
    !! run_tier3's CHECK line.
    integer, intent(out), optional :: smallest_fits
    type(measured_basis_t), intent(in), optional :: basis
    character(len=12), intent(inout), optional :: requested_verdict

    integer :: k, ng, local_dims(3), first_fit
    real(dp) :: per_gpu_gib, workspace_gib, overhead_gib, io_gib
    logical :: exact
    character(len=12) :: verdict
    character(len=96) :: reason

    call print_table_header()
    first_fit = 0
    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      ! ng_unsupported_reason covers the base checks (ng<1, the static
      ! BC_y multi-gpu rule, z/cy/cx divisibility) shared by both tables;
      ! basis%multi_gpu_ok on top of that additionally excludes ng>1 when
      ! the real build just fell back from cuFFTMp to plain cuFFT for a BC
      ! that needs it, a signal only the real-build path has (see
      ! report_measured_table).
      reason = ng_unsupported_reason(ng)
      if (present(basis)) then
        if (ng > 1 .and. multi_gpu_supported .and. &
            .not. basis%multi_gpu_ok) &
          reason = 'not supported for this BC/environment'
      end if
      if (len_trim(reason) > 0) then
        print '(a,i0,a,a,a)', ' ', ng, '     (skipped: ', trim(reason), ')'
        cycle
      end if
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]
      if (present(basis)) then
        call measured_row(ng, basis, workspace_gib, overhead_gib)
        per_gpu_gib = workspace_gib + overhead_gib
        verdict = classify(per_gpu_gib, card_gib)
      else
        call estimate_for_ng(ng, per_gpu_gib, exact, verdict, workspace_gib, &
                             overhead_gib, io_gib)
        if (ng == 1) then
          estimate_ng1_gib = per_gpu_gib
          io_staging_ng1_gib = io_gib
        end if
      end if
      call print_table_row(ng, local_dims, workspace_gib, overhead_gib, &
                           per_gpu_gib, card_gib, verdict)
      if (first_fit == 0 .and. trim(verdict) == 'FITS') first_fit = ng
      if (present(requested_verdict) .and. ng == domain_cfg%nproc_dir(3)) &
        requested_verdict = verdict
    end do
    call print_rule('-')
    if (present(smallest_fits)) smallest_fits = first_fit
  end subroutine print_ng_table

  subroutine report()
    integer :: smallest_fits

    call print_rule('=')
    print '(a,i0,a,i0,a,i0,a)', 'Input grid: ', gdims(1), 'x', gdims(2), &
      'x', gdims(3)
    print '(a,f0.2,a,f0.2,a)', 'Card memory: ', card_gib, ' GiB (', &
      card_free_gib, ' GiB free now)'
    if (card_gib - card_free_gib > 0.5_dp) &
      print '(a,f0.2,a)', 'WARNING: ', card_gib - card_free_gib, ' GiB of &
        &this card is in use by other processes; verdicts are against the &
        &full card, and the real build only proceeds if it fits in what is &
        &free.'
    print '(a,i0)', 'peak_fields (static estimate): ', peak_fields
    call print_io_staging_line()
    select case (trim(run_mode))
    case ('STATIC')
      print '(a)', 'Mode: static (Tier 1 only, no FFT plan query)'
    case ('BUILD')
      print '(a)', 'Mode: build (Tier 1 + Tier 2, then a real Tier 3 &
        &build/measure)'
    case default
      print '(a)', 'Mode: default (Tier 1 + Tier 2 FFT plan query)'
    end select
    call print_rule('=')
    call print_requested_line()

    call print_rule('-')
    call print_ng_table(smallest_fits)
    print '(a)', 'Notes: workspace + overhead = per-GPU. % is per-GPU &
      &against this card''s total memory.'

    if (smallest_fits > 0) then
      print '(a,i0,a)', 'Smallest GPU count that fits: ', smallest_fits, '.'
    else
      print '(a)', 'No scanned GPU count fits with headroom to spare.'
    end if

    if (trim(final_verdict) == 'BORDERLINE' .and. trim(run_mode) /= 'BUILD') &
      print '(a)', 'To confirm with a real measurement, re-run with --build.'
  end subroutine report

  function scratch_name(parent) result(path)
    !! <parent>/x3d2-memcheck-build.<pid>, the scratch directory of this
    !! process. make_scratch creates it and leave_build_scratch
    !! re-derives it before removing anything, so both must get the name
    !! from here.
    character(len=*), intent(in) :: parent
    character(len=4096) :: path
    integer(c_int) :: pid
    character(len=32) :: pid_str

    pid = c_getpid()
    write (pid_str, '(i0)') pid
    path = trim(parent)//'/x3d2-memcheck-build.'//trim(pid_str)
  end function scratch_name

  subroutine run_sh(cmd, status)
    !! Run cmd through the shell. status (optional) receives its exit
    !! status, 0 on success; leave it out for a best effort cleanup command.
    character(len=*), intent(in) :: cmd
    integer, intent(out), optional :: status

    call execute_command_line(cmd, exitstat=status)
  end subroutine run_sh

  subroutine skip_build(msg, cleanup, ok)
    !! Common tail of every "Real build skipped" exit of the scratch setup:
    !! print msg, remove the scratch directory made so far if cleanup is
    !! set, and flag ok=.false. so run_tier3 stops (the estimate above
    !! stands).
    character(len=*), intent(in) :: msg
    logical, intent(in) :: cleanup
    logical, intent(out) :: ok

    print '(a)', msg
    if (cleanup) call run_sh("rm -rf '"//trim(scratch_dir)//"'")
    ok = .false.
  end subroutine skip_build

  subroutine validate_build_inputs(ok)
    !! Validates the known failure causes of a real build up front, before
    !! touching the filesystem: the input file must exist,
    !! domain_cfg%flow_case_name must be one this tool (and
    !! build_and_measure's own select case) can dispatch, and (ibm_on=T) the
    !! matching ibm_<BC-suffix>.bp mask file must already be present in the
    !! invoking directory. Records the invoking directory in orig_dir for
    !! make_scratch. Sets ok=.false. (and the estimate above stands, like
    !! the other run_tier3 skips) on any failure.
    logical, intent(out) :: ok

    character(len=16) :: ibm_file
    logical :: ibm_file_exists, input_exists, flow_case_supported

    ok = .true.
    orig_dir = current_dir()

    inquire (file=trim(input_path), exist=input_exists)
    if (.not. input_exists) then
      call skip_build('Real build skipped: input file not found: '// &
                      trim(input_path), .false., ok)
      return
    end if

    select case (trim(domain_cfg%flow_case_name))
    case ('tgv', 'generic', 'channel', 'cylinder')
      flow_case_supported = .true.
    case default
      flow_case_supported = .false.
    end select
    if (.not. flow_case_supported) then
      call skip_build("Real build skipped: flow case '"// &
                      trim(domain_cfg%flow_case_name)// &
                      "' has no dispatch in x3d2-memcheck", .false., ok)
      return
    end if

    if (solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(periodic_x, periodic_y, periodic_z)
      inquire (file=trim(orig_dir)//'/'//trim(ibm_file), &
              exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        call skip_build('Real build skipped: ibm_on=T but the matching &
                        &ibm_<BC-suffix>.bp mask file was not found in the &
                        &working directory; the estimate above stands.', &
                        .false., ok)
        return
      end if
    end if
  end subroutine validate_build_inputs

  subroutine make_scratch(ok)
    !! Real build runs inside a throwaway x3d2-memcheck-build.<pid>
    !! subdirectory of the invoking directory, because a case build
    !! initialises monitoring (writes monitoring.csv) and calls
    !! postprocess(0), which clobbered run directories on 2026-09-15. A
    !! relative input path is mirrored into the scratch directory by
    !! symlink; a relative path containing a '..' component (or a single
    !! quote, which the shell quoting below cannot handle) cannot be
    !! mirrored this way and is rejected - pass an absolute path instead.
    !! Sets ok=.false. (and the estimate above stands, like the other
    !! run_tier3 skips) on any failure; leave_build_scratch removes the
    !! scratch directory once the build finishes. validate_build_inputs
    !! must have run first. If an unexpected error stop happens after the
    !! chdir below regardless - a failure mode the pre-validation does not
    !! cover - the scratch directory is left behind under the invoking
    !! directory; the next run with the same pid in the same directory
    !! removes it as a stale leftover before creating its own.
    logical, intent(out) :: ok

    integer(c_int) :: rc
    integer :: st, slash_pos
    character(len=16) :: ibm_file
    logical :: input_exists

    ok = .true.
    scratch_dir = scratch_name(orig_dir)

    ! A stale scratch directory from an earlier run's unexpected error stop
    ! after the chdir below (see this subroutine's docstring) would make a
    ! plain mkdir fail - remove it first, if present.
    call run_sh("test -d '"//trim(scratch_dir)//"'", st)
    if (st == 0) then
      call run_sh("rm -rf '"//trim(scratch_dir)//"'", st)
      if (st /= 0) then
        call skip_build('Real build skipped: could not remove stale &
                        &scratch directory '//trim(scratch_dir), .false., ok)
        return
      end if
      print '(a,a)', 'Removed stale scratch directory ', trim(scratch_dir)
    end if

    call run_sh("mkdir '"//trim(scratch_dir)//"'", st)
    if (st /= 0) then
      call skip_build('Real build skipped: could not create scratch &
                      &directory '//trim(scratch_dir), .false., ok)
      return
    end if

    if (input_path(1:1) /= '/') then
      if (has_dotdot_component(trim(input_path)) .or. &
          index(trim(input_path), "'") > 0) then
        call skip_build("Real build skipped: relative input path with '..' &
                        &cannot be mirrored; pass an absolute path", &
                        .true., ok)
        return
      end if
      slash_pos = index(trim(input_path), '/', back=.true.)
      if (slash_pos > 0) then
        call run_sh("mkdir -p '"//trim(scratch_dir)//'/'// &
                    trim(input_path(1:slash_pos - 1))//"'", st)
        if (st /= 0) then
          call skip_build('Real build skipped: could not prepare scratch &
                          &directory (mkdir of the input''s parent failed)', &
                          .true., ok)
          return
        end if
      end if
      call run_sh("ln -s '"//trim(orig_dir)//'/'// &
                  trim(input_path)//"' '"//trim(scratch_dir)// &
                  '/'//trim(input_path)//"'", st)
      if (st /= 0) then
        call skip_build('Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
      inquire (file=trim(scratch_dir)//'/'//trim(input_path), &
              exist=input_exists)
      if (.not. input_exists) then
        call skip_build('Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
    end if

    if (solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(periodic_x, periodic_y, periodic_z)
      call run_sh("ln -s '"//trim(orig_dir)//'/'// &
                  trim(ibm_file)//"' '"//trim(scratch_dir)// &
                  '/'//trim(ibm_file)//"'", st)
      if (st /= 0) then
        call skip_build('Real build skipped: could not prepare scratch &
                        &directory (ibm mask symlink failed)', .true., ok)
        return
      end if
    end if

    rc = c_chdir(trim(scratch_dir)//c_null_char)
    if (rc /= 0) then
      call skip_build('Real build skipped: could not chdir into scratch &
                      &directory '//trim(scratch_dir), .true., ok)
      return
    end if

    print '(a,a,a)', 'Real build scratch directory: ', trim(scratch_dir), &
      ' (removed after the build)'
  end subroutine make_scratch

  subroutine leave_build_scratch()
    !! Restore the invoking directory and remove the scratch directory
    !! make_scratch created. Called after build_and_measure returns,
    !! on both the normal path and the ibm_missing early return.
    integer(c_int) :: rc
    character(len=4096) :: expected_scratch_dir

    rc = c_chdir(trim(orig_dir)//c_null_char)
    if (rc /= 0) &
      error stop 'x3d2-memcheck: could not chdir back to the invoking &
        &directory after the real build; state is unknown, not removing &
        &the scratch directory.'

    ! Never remove anything other than the exact scratch directory
    ! make_scratch created and chdir'd into.
    expected_scratch_dir = scratch_name(orig_dir)
    if (trim(scratch_dir) == trim(expected_scratch_dir)) &
      call run_sh("rm -rf '"//trim(scratch_dir)//"'")

    orig_dir = ''
    scratch_dir = ''
  end subroutine leave_build_scratch

  subroutine print_check_line(used_gib)
    !! --build's CHECK line: the static ng=1 estimate (without the IO
    !! staging term) against the measured device memory in use, OK within
    !! CHECK_TOLERANCE and MISMATCH beyond it.
    real(dp), intent(in) :: used_gib

    real(dp) :: pct_error, estimate_ng1_excl_io_gib
    character(len=8) :: check_word

    ! The real build above performs no snapshot/checkpoint write, so the
    ! GPU-aware IO staging term (analytical only, never measured here) is
    ! excluded from both the estimate compared and the printed CHECK line.
    estimate_ng1_excl_io_gib = estimate_ng1_gib - io_staging_ng1_gib
    pct_error = 100._dp*(estimate_ng1_excl_io_gib - used_gib)/used_gib
    if (abs(estimate_ng1_excl_io_gib - used_gib) <= &
        CHECK_TOLERANCE*used_gib) then
      check_word = 'OK'
    else
      check_word = 'MISMATCH'
    end if
    print '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a,a,a)', 'CHECK ng=1: estimated ', &
      estimate_ng1_excl_io_gib, ' GiB, measured ', used_gib, ' GiB (', &
      pct_error, '%) - ', trim(check_word), ' (tolerance 5%)'
    if (io_staging_ng1_gib > 0._dp) &
      print '(a,f0.1,a)', '  (GPU-aware IO staging ', &
        io_staging_ng1_gib*1024._dp, ' MiB excluded from CHECK: the real &
        &build performs no snapshot/checkpoint write.)'
  end subroutine print_check_line

  subroutine run_tier3()
    !! Tier 3: a REAL case build + one substep on 1 GPU, measured with
    !! cudaMemGetInfo - the last-resort ground truth, now run in-process
    !! (see the top-of-file note on why this removes the old nested mpirun
    !! limitation). Skips the build entirely (the estimate above stands)
    !! if the fields only workspace alone already exceeds BUILD_FRACTION of
    !! what is currently FREE on the card (card_free_gib, not the card's
    !! full capacity - see its declaration), matching this tool's real
    !! build probe, or if ibm_on=T and the matching mask file is not
    !! present in the working directory. The real build itself runs inside
    !! a throwaway x3d2-memcheck-build.<pid> scratch directory
    !! (make_scratch/leave_build_scratch above), because it
    !! initialises monitoring (writes monitoring.csv) and calls
    !! postprocess(0), which clobbered run directories on 2026-09-15;
    !! make_scratch also skips the build
    !! (ok=.false.) if the scratch directory cannot be created or entered,
    !! or if input_path is relative with a '..' component it cannot mirror.
    !! validate_build_inputs pre-validates the known causes up front - a
    !! missing input file, an unsupported flow case, or (ibm_on=T) a
    !! missing mask file - before touching the filesystem; see
    !! make_scratch's docstring for what happens if an unexpected error stop
    !! occurs after its chdir regardless (a scratch directory left behind
    !! under the invoking directory, cleaned up by the next run with the
    !! same pid).
    real(dp) :: used_gib, workspace_gib_measured, ws_guess
    integer :: measured_peak_fields
    logical :: ibm_missing, build_scratch_ok
    character(len=12) :: measured_verdict

    ! --build gate: skip the real build (after printing why) when the fields
    ! only workspace alone already exceeds BUILD_FRACTION of the FREE card
    ! memory.
    ws_guess = to_gib(fields_plus_halo_bytes_n(1, peak_fields))
    ! fields_plus_halo_bytes_n includes the halo term; the historical
    ! BUILD_FRACTION guard was calibrated against fields alone, so subtract
    ! it back out here rather than changing the (already-validated)
    ! threshold itself.
    ws_guess = ws_guess - to_gib(padded_halo_bytes(gdims, SZ, n_halo))
    if (ws_guess >= BUILD_FRACTION*card_free_gib) then
      call print_rule('-')
      print '(a,f0.2,a,f0.1,a,f0.2,a)', 'Real build skipped: fields only &
        &workspace ', ws_guess, ' GiB exceeds ', 100._dp*BUILD_FRACTION, &
        '% of the FREE card memory (', card_free_gib, ' GiB free); the &
        &estimate above stands.'
      return
    end if

    if (ng1_query_ran) then
      call print_rule('-')
      print '(a,f0.2,a)', 'Residual after plan query: ', &
        ng1_query_residual_gib, ' GiB (should be ~0 if cufftDestroy &
        &released it; if not, the real build below may double-reserve the &
        &NVSHMEM heap).'
    end if

    call validate_build_inputs(build_scratch_ok)
    if (.not. build_scratch_ok) return
    call make_scratch(build_scratch_ok)
    if (.not. build_scratch_ok) return

    call build_and_measure(gdims, measured_peak_fields, used_gib, &
                           ibm_missing)
    call leave_build_scratch()
    if (ibm_missing) then
      call print_rule('-')
      print '(a)', 'Real build skipped: ibm_on=T but the matching &
        &ibm_<BC-suffix>.bp mask file was not found in the working &
        &directory; the estimate above stands.'
      return
    end if

    call check_peak_fields_table(measured_peak_fields)

    workspace_gib_measured = &
      to_gib(fields_plus_halo_bytes_n(1, measured_peak_fields))

    call report_measured_table(measured_peak_fields, used_gib, &
                               workspace_gib_measured, measured_verdict)
    final_verdict = measured_verdict

    call print_check_line(used_gib)
  end subroutine run_tier3

  subroutine check_peak_fields_table(measured)
    !! Drift detector for m_memory_estimate's peak_fields_lookup: compares
    !! its static, no build estimate against the real measured peak_fields
    !! (allocator%next_id) for the config this run just built, and warns
    !! (does not fail the tool) on any mismatch. This is what lets the
    !! static table grow/self-correct as the solver evolves, instead of
    !! silently going stale the next time someone adds a new persistent
    !! allocation to a code path the table already covers.
    integer, intent(in) :: measured
    integer :: predicted

    if (irank /= 0) return

    predicted = peak_fields_lookup(trim(domain_cfg%flow_case_name), &
                                   solver_cfg%n_species, &
                                   trim(les_cfg%model) /= 'none', &
                                   solver_cfg%ibm_on, output_vorticity, &
                                   output_qcriterion, &
                                   solver_cfg%lowmem_transeq)
    if (predicted /= measured) then
      print '(a,i0,a,i0,a)', &
        'WARNING: peak_fields_lookup predicted ', predicted, &
        ' but the real build measured ', measured, &
        ' - m_memory_estimate''s table is stale for this config and &
        &should be recalibrated (see src/memory_estimate.f90).'
    end if
  end subroutine check_peak_fields_table

  subroutine measured_row(ng, basis, w_local, overhead_term)
    !! Workspace and overhead (GiB) of one per-ng row of the --build
    !! measured table. "workspace" is the exact fields+halo term with the
    !! MEASURED peak_fields; "overhead" starts from the measured ng=1
    !! overhead and applies the same exact ng-scaling corrections
    !! estimate_for_ng encodes (see report_measured_table).
    integer, intent(in) :: ng
    type(measured_basis_t), intent(in) :: basis
    real(dp), intent(out) :: w_local, overhead_term

    real(dp) :: mirror_gib
    integer(i8) :: spec_bytes_1

    spec_bytes_1 = spectral_slab_bytes(bc_is_100, bc_is_110, cdims, 1)
    w_local = to_gib(fields_plus_halo_bytes_n(ng, basis%peak_fields))

    mirror_gib = 0._dp
    if (bc_is_100) &
      mirror_gib = to_gib(mirror_buffer_bytes_100(cdims, ng))

    overhead_term = basis%overhead_gib + &
                    to_gib(spectral_slab_bytes(bc_is_100, bc_is_110, cdims, &
                                               ng) - spec_bytes_1) &
                    + mirror_gib
    if (use_cufftmp_known .and. use_cufftmp) &
      overhead_term = overhead_term + &
                      to_gib(context_floor_bytes(ng, .true.) - &
                             context_floor_bytes(1, .true.))
    overhead_term = overhead_term + &
                    to_gib(gpu_io_staging_bytes([gdims(1), gdims(2), &
                                                 gdims(3)/ng], &
                                                checkpoint_cfg, &
                                                gpu_io_device_write))
  end subroutine measured_row

  subroutine report_measured_table(measured_peak_fields, used_gib, &
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
    integer, intent(in) :: measured_peak_fields
    real(dp), intent(in) :: used_gib, workspace_gib_measured
    character(len=12), intent(out) :: requested_verdict

    logical :: multi_gpu_supported_measured
    type(measured_basis_t) :: basis

    requested_verdict = 'DOES_NOT_FIT'

    multi_gpu_supported_measured = .true.
    if (bc_is_010 .or. bc_is_110) then
      ! BC_y non-periodic: no nproc>1 support for this BC
      ! (src/poisson_fft.f90:178-180,196-198).
      multi_gpu_supported_measured = .false.
    else if ((bc_is_100 .or. bc_is_000) .and. use_cufftmp_known .and. &
            (.not. use_cufftmp)) then
      ! Covers the 100 case (needs cuFFTMp to decompose at all) and the
      ! fully-periodic 000 case (the plain-cuFFT fallback performs a purely
      ! local per-rank transform with no cross-rank exchange -
      ! src/backend/cuda/poisson_fft.f90:721-740 - so at nproc>1 it would
      ! run without erroring but compute invalid results). cuFFTMp
      ! unavailable here means the build just performed fell back to plain
      ! cuFFT.
      multi_gpu_supported_measured = .false.
    end if

    basis = measured_basis_t(measured_peak_fields, &
                             used_gib - workspace_gib_measured, &
                             multi_gpu_supported_measured)

    call print_rule('-')
    print '(a,i0,a,i0,a,f0.2,a,f0.2,a)', 'Real build (1 GPU, ', &
      measured_n_substeps, ' substep(s)): peak_fields measured ', &
      measured_peak_fields, ', workspace ', workspace_gib_measured, &
      ' GiB, used ', used_gib, ' GiB'
    ! Single machine-parseable line for --extensive sweeps (see
    ! scripts/memcheck_extensive_sweep.sh), so a wrapper script comparing
    ! several substep counts can grep one line per run instead of parsing
    ! the whole table.
    if (extensive_substeps > 0) &
      print '(a,i0,a,i0,a,f0.3,a,f0.3)', 'EXTENSIVE result: substeps=', &
        measured_n_substeps, ' peak_fields=', measured_peak_fields, &
        ' used_gib=', used_gib, ' pct_card=', 100._dp*used_gib/card_gib
    call print_ng_table(basis=basis, requested_verdict=requested_verdict)
  end subroutine report_measured_table

  subroutine build_mesh(dims_in, mesh, dims, ibm_missing)
    !! Mesh of the real build at dims_in, its vertex dims, and the ibm_on
    !! mask pre-check. ibm_missing=.true. tells build_and_measure to stop
    !! before building anything.
    integer, intent(in) :: dims_in(3)
    type(mesh_t), intent(out) :: mesh
    integer, intent(out) :: dims(3)
    logical, intent(out) :: ibm_missing

    character(len=16) :: ibm_file
    logical :: ibm_file_exists

    ibm_missing = .false.
    mesh = mesh_t(dims_in, [1, 1, 1], domain_cfg%L_global, &
                  domain_cfg%BC_x, domain_cfg%BC_y, domain_cfg%BC_z, &
                  domain_cfg%stretching, domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)

    ! ibm_on triggers reading an external ibm_<BC-suffix>.bp mask file
    ! inside solver init (src/module/ibm.f90), which MPI_Aborts if that
    ! file is missing. That is correct for production xcompact, but not
    ! for a tool that should degrade gracefully - check for it here, using
    ! the same suffix construction as src/module/ibm.f90:70-75 (shared via
    ! ibm_mask_filename, also used by make_scratch to mirror the
    ! mask into the build scratch directory), and bail out before
    ! triggering the case/solver construction that would abort.
    if (solver_cfg%ibm_on) then
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
    !! Build the device allocator, the host allocator and the CUDA backend
    !! of the real build, and point allocator/backend at them. The three
    !! objects belong to the caller (build_and_measure) because the
    !! pointers, the backend and the case built on it must not outlive
    !! them.
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

  subroutine make_flow_case(backend, mesh, host_allocator, flow_case)
    !! Build the flow case named in the input through the same case
    !! dispatch select case xcompact.f90 uses.
    class(base_backend_t), pointer, intent(in) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    class(base_case_t), allocatable, intent(out) :: flow_case

    select case (trim(domain_cfg%flow_case_name))
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

  subroutine drive_case(flow_case)
    !! Drive the freshly built flow case the way run() does so the
    !! allocator reaches its work-field high-water mark: the pre-loop
    !! postprocess hook, one (or --extensive <n>) substep(s), then the
    !! per-iteration pressure and derived-field calls. Records the substep
    !! count in measured_n_substeps and the derived-field gating in
    !! output_vorticity/output_qcriterion.
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
    if (extensive_substeps > 0) n_substeps = extensive_substeps
    do i_substep = 1, n_substeps
      call flow_case%substep(curr, deriv, i_substep)
    end do
    measured_n_substeps = n_substeps

    ! Mirrors run()'s per-iteration postprocessing calls (base_case.f90,
    ! immediately after the sub-stage loop): keep_pressure and
    ! vorticity/Q-criterion allocations are gated HERE, not inside
    ! postprocess() above, so they must be driven separately to be
    ! captured - postprocess(0,.) alone does not reach them.
    if (flow_case%solver%keep_pressure) &
      call compute_pressure_vert(flow_case%solver)
    output_vorticity = output_field_active(flow_case%io_mgr%snapshot_mgr% &
                                           config, 'vorticity')
    output_qcriterion = output_field_active(flow_case%io_mgr%snapshot_mgr% &
                                            config, 'qcriterion')
    if (output_vorticity .or. output_qcriterion) &
      call compute_derived_fields(flow_case%solver, output_vorticity, &
                                  output_qcriterion)
  end subroutine drive_case

  subroutine build_and_measure(dims_in, npeak, dev_used, ibm_missing)
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
    integer(kind=cuda_count_kind) :: free_b, total_b

    npeak = 0
    dev_used = 0._dp
    output_vorticity = .false.
    output_qcriterion = .false.

    call build_mesh(dims_in, mesh, dims, ibm_missing)
    if (ibm_missing) return

    call make_cuda_backend(mesh, dims, cuda_allocator, host_allocator, &
                           cuda_backend, allocator, backend)

    call make_flow_case(backend, mesh, host_allocator, flow_case)

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_measured_table's safety gating on
    ! the 100 case, since there is no static "is cuFFTMp available" query.
    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      use_cufftmp = pf%use_cufftmp
      use_cufftmp_known = .true.
    end select

    call drive_case(flow_case)

    npeak = allocator%next_id
    ierr = cudaMemGetInfo(free_b, total_b)
    dev_used = to_gib(int(total_b - free_b, i8))
  end subroutine build_and_measure

end program x3d2_memcheck
