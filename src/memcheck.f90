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
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice
  use m_common, only: dp, i8
  use m_mesh, only: mesh_t
  use m_cuda_common, only: SZ
  use m_memory_estimate, only: output_field_active, peak_fields_lookup, &
                               padded_halo_bytes
  use m_cuda_memory_estimate, only: check_status
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t
  use m_memcheck_context, only: memcheck_ctx_t, parse_args, read_config, &
                                classify_bc, resolve_gpu_io_mode
  use m_memcheck_estimate, only: to_gib, n_halo, fields_plus_halo_bytes_n
  use m_memcheck_report, only: print_rule, report
  use m_memcheck_scratch, only: validate_build_inputs, make_scratch, &
                                leave_build_scratch
  use m_memcheck_build, only: check_peak_fields_table, report_measured_table, &
                              build_mesh, make_cuda_backend, make_flow_case, &
                              drive_case

  implicit none

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

  type(memcheck_ctx_t) :: ctx
  integer :: ierr, nproc, ndevs

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, ctx%irank, ierr)
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

  call parse_args(ctx)
  call read_config(ctx)
  call resolve_gpu_io_mode(ctx)
  call classify_bc(ctx)
  call query_card_gib()

  ctx%peak_fields = peak_fields_lookup(trim(ctx%domain_cfg%flow_case_name), &
                                   ctx%solver_cfg%n_species, &
                                   trim(ctx%les_cfg%model) /= 'none', &
                                   ctx%solver_cfg%ibm_on, &
                                   output_field_active(ctx%checkpoint_cfg, &
                                                       'vorticity'), &
                                   output_field_active(ctx%checkpoint_cfg, &
                                                       'qcriterion'), &
                                   ctx%solver_cfg%lowmem_transeq)

  ctx%multi_gpu_supported = .not. (ctx%bc_is_010 .or. ctx%bc_is_110)

  call report(ctx)
  if (trim(ctx%run_mode) == 'BUILD') call run_tier3()

  call MPI_Finalize(ierr)

  ! Exit code mirrors the final verdict: 0=FITS, 1=BORDERLINE,
  ! 2=DOES_NOT_FIT (also covers UNSUPPORTED - an unsupported nproc_dir
  ! request, via the default case below, since neither fits nor merely
  ! borderline is the right signal for it).
  select case (trim(ctx%final_verdict))
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

  subroutine query_card_gib()
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo (card total)')
    ctx%card_gib = to_gib(int(total_b, i8))
    ctx%card_free_gib = to_gib(int(free_b, i8))
  end subroutine query_card_gib

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
    real(dp) :: pct_error, estimate_ng1_excl_io_gib
    character(len=8) :: check_word
    integer :: measured_peak_fields
    logical :: ibm_missing, build_scratch_ok
    character(len=12) :: measured_verdict

    ! --build gate: skip the real build (after printing why) when the fields
    ! only workspace alone already exceeds BUILD_FRACTION of the FREE card
    ! memory.
    ws_guess = to_gib(fields_plus_halo_bytes_n(ctx, 1, ctx%peak_fields))
    ! fields_plus_halo_bytes_n includes the halo term; the historical
    ! BUILD_FRACTION guard was calibrated against fields alone, so subtract
    ! it back out here rather than changing the (already-validated)
    ! threshold itself.
    ws_guess = ws_guess - to_gib(padded_halo_bytes(ctx%gdims, SZ, n_halo))
    if (ws_guess >= BUILD_FRACTION*ctx%card_free_gib) then
      call print_rule('-')
      print '(a,f0.2,a,f0.1,a,f0.2,a)', 'Real build skipped: fields only &
        &workspace ', ws_guess, ' GiB exceeds ', 100._dp*BUILD_FRACTION, &
        '% of the FREE card memory (', ctx%card_free_gib, ' GiB free); the &
        &estimate above stands.'
      return
    end if

    if (ctx%ng1_query_ran) then
      call print_rule('-')
      print '(a,f0.2,a)', 'Residual after plan query: ', &
        ctx%ng1_query_residual_gib, ' GiB (should be ~0 if cufftDestroy &
        &released it; if not, the real build below may double-reserve the &
        &NVSHMEM heap).'
    end if

    call validate_build_inputs(ctx, build_scratch_ok)
    if (.not. build_scratch_ok) return
    call make_scratch(ctx, build_scratch_ok)
    if (.not. build_scratch_ok) return

    call build_and_measure(ctx%gdims, measured_peak_fields, used_gib, &
                           ibm_missing)
    call leave_build_scratch(ctx)
    if (ibm_missing) then
      call print_rule('-')
      print '(a)', 'Real build skipped: ibm_on=T but the matching &
        &ibm_<BC-suffix>.bp mask file was not found in the working &
        &directory; the estimate above stands.'
      return
    end if

    call check_peak_fields_table(ctx, measured_peak_fields)

    workspace_gib_measured = &
      to_gib(fields_plus_halo_bytes_n(ctx, 1, measured_peak_fields))

    call report_measured_table(ctx, measured_peak_fields, used_gib, &
                               workspace_gib_measured, measured_verdict)
    ctx%final_verdict = measured_verdict

    ! --build's CHECK line: the static ng=1 estimate against the measured
    ! device memory in use, OK within CHECK_TOLERANCE and MISMATCH beyond it.
    ! The real build above performs no snapshot/checkpoint write, so the
    ! GPU-aware IO staging term (analytical only, never measured here) is
    ! excluded from both the estimate compared and the printed CHECK line.
    estimate_ng1_excl_io_gib = ctx%estimate_ng1_gib - ctx%io_staging_ng1_gib
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
    if (ctx%io_staging_ng1_gib > 0._dp) &
      print '(a,f0.1,a)', '  (GPU-aware IO staging ', &
        ctx%io_staging_ng1_gib*1024._dp, ' MiB excluded from CHECK: the real &
        &build performs no snapshot/checkpoint write.)'
  end subroutine run_tier3

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
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    npeak = 0
    dev_used = 0._dp
    ctx%output_vorticity = .false.
    ctx%output_qcriterion = .false.

    call build_mesh(ctx, dims_in, mesh, dims, ibm_missing)
    if (ibm_missing) return

    call make_cuda_backend(mesh, dims, cuda_allocator, host_allocator, &
                           cuda_backend, allocator, backend)

    call make_flow_case(ctx, backend, mesh, host_allocator, flow_case)

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_measured_table's safety gating on
    ! the 100 case, since there is no static "is cuFFTMp available" query.
    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      ctx%use_cufftmp = pf%use_cufftmp
      ctx%use_cufftmp_known = .true.
    end select

    call drive_case(ctx, flow_case)

    npeak = allocator%next_id
    ierr = cudaMemGetInfo(free_b, total_b)
    dev_used = to_gib(int(total_b - free_b, i8))
  end subroutine build_and_measure

end program x3d2_memcheck
