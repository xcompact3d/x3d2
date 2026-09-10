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
  !!
  !! Usage: x3d2-memcheck <input.x3d> [--static | --build]
  use mpi
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice
  use m_common, only: dp, i8, nbytes, VERT
  use m_config, only: domain_config_t, solver_config_t, les_config_t, &
                      checkpoint_config_t
  use m_mesh, only: mesh_t, periodic_dir
  use m_cuda_common, only: SZ
  use m_memory_estimate, only: padded_dim, padded_cells, cell_dims, &
                               spectral_slab_bytes, mirror_buffer_bytes_100, &
                               output_field_active, peak_fields_lookup, &
                               halo_bytes
  use m_cuda_memory_estimate, only: fft_workspace_bytes_query, &
                                    context_floor_bytes
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
  !> exceeds this fraction of the card - leaves headroom for context + FFT
  !> scratch so the build never OOMs. Same value and role this tool's real
  !> build probe has always used.
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

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 1) error stop 'x3d2-memcheck: run single-rank (mpirun -n 1).'

  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(0)

  call parse_args()
  call read_config()
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
  ! 2=DOES_NOT_FIT.
  select case (trim(final_verdict))
  case ('FITS')
    call exit(0)
  case ('BORDERLINE')
    call exit(1)
  case default
    call exit(2)
  end select

contains

  subroutine parse_args()
    integer :: i, nargs
    character(len=256) :: arg
    logical :: static_flag, build_flag

    nargs = command_argument_count()
    if (nargs < 1) error stop 'usage: x3d2-memcheck <input.x3d> &
      &[--static | --build]'
    ! The input file must stay positional argument 1: the real case
    ! constructors this tool calls under --build (case_channel_init,
    ! case_cylinder_init, solver_init) independently re-read
    ! get_argument(1) themselves rather than taking the path as a
    ! parameter, so flags must only ever appear after it.
    call get_command_argument(1, input_path)

    static_flag = .false.
    build_flag = .false.
    do i = 2, nargs
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--static')
        static_flag = .true.
      case ('--build')
        build_flag = .true.
      case default
        error stop 'x3d2-memcheck: unknown flag '//trim(arg)// &
          ' (expected --static or --build)'
      end select
    end do
    if (static_flag .and. build_flag) &
      error stop 'x3d2-memcheck: --static and --build are mutually &
        &exclusive (--static stops at Tier 1, --build adds Tier 3).'

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

  subroutine query_card_gib()
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    card_gib = real(total_b, dp)/1024._dp**3
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
    integer :: local_dims(3), padded_dims(3)

    local_dims = [gdims(1), gdims(2), gdims(3)/ng]
    padded_dims = [padded_dim(local_dims(1), SZ), &
                   padded_dim(local_dims(2), SZ), local_dims(3)]
    nbytes8 = int(npeak, i8)*padded_cells(local_dims, SZ) &
             *int(nbytes, i8) + halo_bytes(padded_dims, SZ, n_halo)
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
    integer, intent(in) :: ng
    integer(i8) :: nbytes8

    nbytes8 = spectral_slab_bytes(bc_is_100, bc_is_110, cdims, ng)
    if (bc_is_100) then
      nbytes8 = nbytes8 + mirror_buffer_bytes_100(cdims, ng)
    else if (bc_is_110) then
      nbytes8 = nbytes8 + spectral_slab_bytes(bc_is_100, bc_is_110, cdims, ng) &
               + int(cdims(1), i8)*int(cdims(2), i8)*int(cdims(3), i8) &
                *int(nbytes, i8)
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
    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   irank == 0, ng1_worksize_bytes, &
                                   ng1_heap_bytes, ng1_xtdesc_bytes, &
                                   ng1_used_cufftmp)
    ierr = cudaMemGetInfo(free_after, total_b)
    ng1_query_residual_gib = real(int(free_before, i8) - &
                                  int(free_after, i8), dp)/1024._dp**3
    ng1_query_ran = .true.
    done = .true.
  end subroutine ensure_fft_query

  subroutine estimate_for_ng(ng, per_gpu_gib, exact, verdict, workspace_gib, &
                             overhead_gib)
    !! Tier 1 floor first; only calls into Tier 2 (a real, throwaway GPU
    !! FFT plan) when Tier 1 alone cannot already answer DOES_NOT_FIT, and
    !! never under --static. workspace_gib/overhead_gib (optional) are the
    !! same two-term breakdown the --build measured table reports:
    !! workspace = exact fields+halo+spectral+mirror (Tier 1, no GPU
    !! needed); overhead = everything from the Tier 2 GPU query and the
    !! hardware constants (worksize/xtdesc/heap/context).
    integer, intent(in) :: ng
    real(dp), intent(out) :: per_gpu_gib
    logical, intent(out) :: exact
    character(len=12), intent(out) :: verdict
    real(dp), intent(out), optional :: workspace_gib, overhead_gib

    integer(i8) :: base_bytes, spec_bytes, worksize_bytes, floor_bytes, &
                  xtdesc_bytes, heap_bytes, context_bytes
    logical :: used_cufftmp
    real(dp) :: floor_gib

    base_bytes = fields_plus_halo_bytes(ng)
    spec_bytes = spectral_plus_mirror_bytes(ng)
    if (present(workspace_gib)) &
      workspace_gib = real(base_bytes + spec_bytes, dp)/1024._dp**3
    ! Floor: assume cuFFTMp is used (the realistic case for every BC this
    ! solver ever attempts it for - see m_cuda_memory_estimate's 110 guard)
    ! since that gives the larger, safer lower bound.
    floor_bytes = base_bytes + spec_bytes + &
                  context_floor_bytes(ng, .not. bc_is_110)
    floor_gib = real(floor_bytes, dp)/1024._dp**3

    if (trim(run_mode) == 'STATIC') then
      ! --static: Tier 1 only, never queries the GPU FFT plan. floor_gib
      ! IS the estimate here - not just a DOES_NOT_FIT early-return floor
      ! - so classify it with the full three-way verdict.
      per_gpu_gib = floor_gib
      exact = .false.
      if (present(overhead_gib)) &
        overhead_gib = real(context_floor_bytes(ng, .not. bc_is_110), dp) &
                       /1024._dp**3
      verdict = classify(per_gpu_gib, card_gib)
      return
    end if

    if (floor_gib > DOES_NOT_FIT_FRACTION*card_gib) then
      per_gpu_gib = floor_gib
      exact = .false.
      verdict = 'DOES_NOT_FIT'
      if (present(overhead_gib)) &
        overhead_gib = real(context_floor_bytes(ng, .not. bc_is_110), dp) &
                       /1024._dp**3
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
    per_gpu_gib = real(base_bytes + spec_bytes + worksize_bytes + &
                       xtdesc_bytes + heap_bytes + context_bytes, dp) &
                 /1024._dp**3
    exact = .true.
    if (present(overhead_gib)) &
      overhead_gib = real(worksize_bytes + xtdesc_bytes + heap_bytes + &
                          context_bytes, dp)/1024._dp**3

    verdict = classify(per_gpu_gib, card_gib)
  end subroutine estimate_for_ng

  subroutine report()
    integer :: k, ng, requested_ng, smallest_fits, local_dims(3)
    real(dp) :: per_gpu_gib, requested_gib, workspace_gib, overhead_gib
    logical :: exact
    character(len=12) :: verdict

    print '(a)', '============================================================'
    print '(a,i0,a,i0,a,i0,a)', 'Input grid: ', gdims(1), 'x', gdims(2), &
      'x', gdims(3)
    print '(a,f0.2,a)', 'Card memory: ', card_gib, ' GiB'
    print '(a,i0)', 'peak_fields (static estimate): ', peak_fields
    select case (trim(run_mode))
    case ('STATIC')
      print '(a)', 'Mode: static (Tier 1 only, no FFT plan query)'
    case ('BUILD')
      print '(a)', 'Mode: build (Tier 1 + Tier 2, then a real Tier 3 &
        &build/measure)'
    case default
      print '(a)', 'Mode: default (Tier 1 + Tier 2 FFT plan query)'
    end select
    print '(a)', '============================================================'

    requested_ng = domain_cfg%nproc_dir(3)
    call estimate_for_ng(requested_ng, requested_gib, exact, final_verdict)
    print '(a,i0,a,f0.2,a,f0.1,a,a)', 'Requested nproc_dir gives ng=', &
      requested_ng, ': ', requested_gib, ' GiB/GPU (', &
      100._dp*requested_gib/card_gib, '% of card) - ', trim(final_verdict)
    if (.not. exact) then
      if (trim(run_mode) == 'STATIC') then
        print '(a)', '  (Tier 1 estimate only: --static skips the GPU FFT &
          &plan query.)'
      else
        print '(a)', '  (Tier 1 floor only: this input is well over the &
          &limit, no GPU plan probe was attempted.)'
      end if
    end if

    print '(a)', '------------------------------------------------------------'
    call print_table_header()
    smallest_fits = 0
    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      if (ng > 1 .and. .not. multi_gpu_supported) then
        print '(a,i0,a)', ' ', ng, '     (not supported: BC_y is &
          &non-periodic - src/poisson_fft.f90 error-stops at nproc>1)'
        cycle
      end if
      if (mod(gdims(3), ng) /= 0) then
        print '(a,i0,a,i0,a)', ' ', ng, '     (skipped: z=', gdims(3), &
          ' not divisible)'
        cycle
      end if
      call estimate_for_ng(ng, per_gpu_gib, exact, verdict, workspace_gib, &
                           overhead_gib)
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]
      call print_table_row(ng, local_dims, workspace_gib, overhead_gib, &
                           per_gpu_gib, card_gib, verdict)
      if (ng == 1) estimate_ng1_gib = per_gpu_gib
      if (smallest_fits == 0 .and. trim(verdict) == 'FITS') smallest_fits = ng
    end do
    print '(a)', '------------------------------------------------------------'
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

  subroutine run_tier3()
    !! Tier 3: a REAL case build + one substep on 1 GPU, measured with
    !! cudaMemGetInfo - the last-resort ground truth, now run in-process
    !! (see the top-of-file note on why this removes the old nested mpirun
    !! limitation). Skips the build entirely (the estimate above stands)
    !! if the fields only workspace alone already exceeds BUILD_FRACTION of
    !! the card, exactly as this tool's real build probe always has, or if
    !! ibm_on=T and the matching mask file is not present in the working
    !! directory.
    real(dp) :: ws_guess, used_gib, workspace_gib_measured, pct_error
    integer :: measured_peak_fields
    logical :: ibm_missing
    character(len=12) :: measured_verdict

    ws_guess = real(fields_plus_halo_bytes_n(1, peak_fields), dp) &
              /1024._dp**3
    ! fields_plus_halo_bytes_n includes the halo term; the historical
    ! BUILD_FRACTION guard was calibrated against fields alone, so subtract
    ! it back out here rather than changing the (already-validated)
    ! threshold itself.
    ws_guess = ws_guess - real(halo_bytes([padded_dim(gdims(1), SZ), &
                                           padded_dim(gdims(2), SZ), &
                                           gdims(3)], SZ, n_halo), dp) &
                         /1024._dp**3
    if (ws_guess >= BUILD_FRACTION*card_gib) then
      print '(a)', '------------------------------------------------------------'
      print '(a,f0.2,a,f0.1,a)', 'Real build skipped: fields only workspace &
        &', ws_guess, ' GiB exceeds ', 100._dp*BUILD_FRACTION, &
        '% of the card; the estimate above stands.'
      return
    end if

    if (ng1_query_ran) then
      print '(a)', '------------------------------------------------------------'
      print '(a,f0.2,a)', 'Residual after plan query: ', &
        ng1_query_residual_gib, ' GiB (should be ~0 if cufftDestroy &
        &released it; if not, the real build below may double-reserve the &
        &NVSHMEM heap).'
    end if

    call build_and_measure(gdims, measured_peak_fields, used_gib, &
                           ibm_missing)
    if (ibm_missing) then
      print '(a)', '------------------------------------------------------------'
      print '(a)', 'Real build skipped: ibm_on=T but the matching &
        &ibm_<BC-suffix>.bp mask file was not found in the working &
        &directory; the estimate above stands.'
      return
    end if

    call check_peak_fields_table(measured_peak_fields)

    workspace_gib_measured = real(fields_plus_halo_bytes_n(1, &
                                                            measured_peak_fields), &
                                  dp)/1024._dp**3

    call report_measured_table(measured_peak_fields, used_gib, &
                               workspace_gib_measured, measured_verdict)
    final_verdict = measured_verdict

    pct_error = 100._dp*(estimate_ng1_gib - used_gib)/used_gib
    if (abs(estimate_ng1_gib - used_gib) <= CHECK_TOLERANCE*used_gib) then
      print '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a)', 'CHECK ng=1: estimated ', &
        estimate_ng1_gib, ' GiB, measured ', used_gib, ' GiB (', &
        pct_error, '%) - OK (tolerance 5%)'
    else
      print '(a,f0.2,a,f0.2,a,sp,f0.1,ss,a)', 'CHECK ng=1: estimated ', &
        estimate_ng1_gib, ' GiB, measured ', used_gib, ' GiB (', &
        pct_error, '%) - MISMATCH (tolerance 5%)'
    end if
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
    integer, intent(in) :: measured_peak_fields
    real(dp), intent(in) :: used_gib, workspace_gib_measured
    character(len=12), intent(out) :: requested_verdict

    integer :: k, ng, local_dims(3), requested_ng
    real(dp) :: overhead, per_gpu, overhead_term, mirror_gib, w_local
    integer(i8) :: spec_bytes_1
    logical :: multi_gpu_supported_measured
    character(len=12) :: verdict

    spec_bytes_1 = spectral_slab_bytes(bc_is_100, bc_is_110, cdims, 1)
    overhead = used_gib - workspace_gib_measured
    requested_ng = domain_cfg%nproc_dir(3)
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

    print '(a)', '------------------------------------------------------------'
    print '(a,i0,a,f0.2,a,f0.2,a)', 'Real build (1 GPU, one substep): &
      &peak_fields measured ', measured_peak_fields, ', workspace ', &
      workspace_gib_measured, ' GiB, used ', used_gib, ' GiB'
    call print_table_header()

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      if (ng > 1 .and. .not. multi_gpu_supported_measured) then
        print '(a,i0,a)', ' ', ng, '     (skipped: not supported for this &
          &BC/environment)'
        cycle
      end if
      if (mod(gdims(3), ng) /= 0) then
        print '(a,i0,a,i0,a)', ' ', ng, '     (skipped: z=', gdims(3), &
          ' not divisible)'
        cycle
      end if
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]
      w_local = real(fields_plus_halo_bytes_n(ng, measured_peak_fields), dp) &
               /1024._dp**3

      mirror_gib = 0._dp
      if (bc_is_100) &
        mirror_gib = real(mirror_buffer_bytes_100(cdims, ng), dp) &
                    /1024._dp**3

      overhead_term = overhead + &
                      real(spectral_slab_bytes(bc_is_100, bc_is_110, cdims, &
                                               ng) - spec_bytes_1, dp) &
                      /1024._dp**3 + mirror_gib
      if (use_cufftmp_known .and. use_cufftmp) &
        overhead_term = overhead_term + &
                        real(context_floor_bytes(ng, .true.) - &
                            context_floor_bytes(1, .true.), dp)/1024._dp**3

      per_gpu = w_local + overhead_term
      verdict = classify(per_gpu, card_gib)
      call print_table_row(ng, local_dims, w_local, overhead_term, per_gpu, &
                           card_gib, verdict)
      if (ng == requested_ng) requested_verdict = verdict
    end do
    print '(a)', '------------------------------------------------------------'
  end subroutine report_measured_table

  subroutine build_and_measure(dims_in, npeak, dev_used, ibm_missing)
    !! Build the real flow case at dims_in on this single GPU (via the same
    !! case-dispatch select case xcompact.f90 uses), run one
    !! substep to drive the allocator to its work-field high-water mark,
    !! and return that mark (npeak) plus the absolute device memory in use
    !! (dev_used).
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
    type(flist_t), allocatable :: curr(:), deriv(:)
    integer :: dims(3), i
    type(cuda_allocator_t), target :: cuda_allocator
    type(cuda_backend_t), target :: cuda_backend
    type(allocator_t), target :: host_allocator
    character(len=16) :: ibm_file
    character(len=3) :: bc_suffix
    logical :: ibm_file_exists
    integer(kind=cuda_count_kind) :: free_b, total_b

    ibm_missing = .false.
    npeak = 0
    dev_used = 0._dp
    output_vorticity = .false.
    output_qcriterion = .false.

    mesh = mesh_t(dims_in, [1, 1, 1], domain_cfg%L_global, &
                  domain_cfg%BC_x, domain_cfg%BC_y, domain_cfg%BC_z, &
                  domain_cfg%stretching, domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)

    ! ibm_on triggers reading an external ibm_<BC-suffix>.bp mask file
    ! inside solver init (src/module/ibm.f90), which MPI_Aborts if that
    ! file is missing. That is correct for production xcompact, but not
    ! for a tool that should degrade gracefully - check for it here, using
    ! the same suffix construction as src/module/ibm.f90:70-75, and bail
    ! out before triggering the case/solver construction that would abort.
    if (solver_cfg%ibm_on) then
      bc_suffix(1:1) = '0'
      if (.not. mesh%grid%periodic_BC(1)) bc_suffix(1:1) = '1'
      bc_suffix(2:2) = '0'
      if (.not. mesh%grid%periodic_BC(2)) bc_suffix(2:2) = '1'
      bc_suffix(3:3) = '0'
      if (.not. mesh%grid%periodic_BC(3)) bc_suffix(3:3) = '1'
      ibm_file = "ibm_"//bc_suffix//".bp"
      inquire (file=ibm_file, exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        ibm_missing = .true.
        return
      end if
    end if

    cuda_allocator = cuda_allocator_t(dims, SZ)
    allocator => cuda_allocator
    host_allocator = allocator_t(dims, SZ)
    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend

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

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_measured_table's safety gating on
    ! the 100 case, since there is no static "is cuFFTMp available" query.
    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      use_cufftmp = pf%use_cufftmp
      use_cufftmp_known = .true.
    end select

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
    call flow_case%substep(curr, deriv, 1)

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

    npeak = allocator%next_id
    ierr = cudaMemGetInfo(free_b, total_b)
    dev_used = real(total_b - free_b, dp)/1024._dp**3
  end subroutine build_and_measure

end program x3d2_memcheck
