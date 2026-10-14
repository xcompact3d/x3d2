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
  use iso_c_binding, only: c_int, c_null_char
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice
  use m_common, only: dp, i8, VERT
  use m_mesh, only: mesh_t
  use m_cuda_common, only: SZ
  use m_memory_estimate, only: spectral_slab_bytes, mirror_buffer_bytes_100, &
                               output_field_active, peak_fields_lookup, &
                               padded_halo_bytes, gpu_io_staging_bytes
  use m_cuda_memory_estimate, only: context_floor_bytes, check_status
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
  use m_memcheck_context, only: memcheck_ctx_t, parse_args, read_config, &
                                classify_bc, resolve_gpu_io_mode
  use m_memcheck_estimate, only: to_gib, classify, ng_unsupported_reason, &
                                 n_halo, fields_plus_halo_bytes_n
  use m_memcheck_report, only: print_rule, print_table_header, &
                               print_table_row, n_gpu_list, report
  use m_memcheck_scratch, only: c_chdir, current_dir, has_dotdot_component, &
                                run_sh, scratch_name, skip_build

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

  subroutine query_card_gib()
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo (card total)')
    ctx%card_gib = to_gib(int(total_b, i8))
    ctx%card_free_gib = to_gib(int(free_b, i8))
  end subroutine query_card_gib

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
    ctx%orig_dir = current_dir()

    inquire (file=trim(ctx%input_path), exist=input_exists)
    if (.not. input_exists) then
      call skip_build(ctx, 'Real build skipped: input file not found: '// &
                      trim(ctx%input_path), .false., ok)
      return
    end if

    select case (trim(ctx%domain_cfg%flow_case_name))
    case ('tgv', 'generic', 'channel', 'cylinder')
      flow_case_supported = .true.
    case default
      flow_case_supported = .false.
    end select
    if (.not. flow_case_supported) then
      call skip_build(ctx, "Real build skipped: flow case '"// &
                      trim(ctx%domain_cfg%flow_case_name)// &
                      "' has no dispatch in x3d2-memcheck", .false., ok)
      return
    end if

    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(ctx%periodic_x, ctx%periodic_y, &
                                   ctx%periodic_z)
      inquire (file=trim(ctx%orig_dir)//'/'//trim(ibm_file), &
              exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        call skip_build(ctx, 'Real build skipped: ibm_on=T but the matching &
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
    ctx%scratch_dir = scratch_name(ctx%orig_dir)

    ! A stale scratch directory from an earlier run's unexpected error stop
    ! after the chdir below (see this subroutine's docstring) would make a
    ! plain mkdir fail - remove it first, if present.
    call run_sh("test -d '"//trim(ctx%scratch_dir)//"'", st)
    if (st == 0) then
      call run_sh("rm -rf '"//trim(ctx%scratch_dir)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not remove stale &
                        &scratch directory '//trim(ctx%scratch_dir), &
                        .false., ok)
        return
      end if
      print '(a,a)', 'Removed stale scratch directory ', trim(ctx%scratch_dir)
    end if

    call run_sh("mkdir '"//trim(ctx%scratch_dir)//"'", st)
    if (st /= 0) then
      call skip_build(ctx, 'Real build skipped: could not create scratch &
                      &directory '//trim(ctx%scratch_dir), .false., ok)
      return
    end if

    if (ctx%input_path(1:1) /= '/') then
      if (has_dotdot_component(trim(ctx%input_path)) .or. &
          index(trim(ctx%input_path), "'") > 0) then
        call skip_build(ctx, "Real build skipped: relative input path with &
                        &'..' cannot be mirrored; pass an absolute path", &
                        .true., ok)
        return
      end if
      slash_pos = index(trim(ctx%input_path), '/', back=.true.)
      if (slash_pos > 0) then
        call run_sh("mkdir -p '"//trim(ctx%scratch_dir)//'/'// &
                    trim(ctx%input_path(1:slash_pos - 1))//"'", st)
        if (st /= 0) then
          call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                          &directory (mkdir of the input''s parent failed)', &
                          .true., ok)
          return
        end if
      end if
      call run_sh("ln -s '"//trim(ctx%orig_dir)//'/'// &
                  trim(ctx%input_path)//"' '"//trim(ctx%scratch_dir)// &
                  '/'//trim(ctx%input_path)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
      inquire (file=trim(ctx%scratch_dir)//'/'//trim(ctx%input_path), &
              exist=input_exists)
      if (.not. input_exists) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
    end if

    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(ctx%periodic_x, ctx%periodic_y, &
                                   ctx%periodic_z)
      call run_sh("ln -s '"//trim(ctx%orig_dir)//'/'// &
                  trim(ibm_file)//"' '"//trim(ctx%scratch_dir)// &
                  '/'//trim(ibm_file)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (ibm mask symlink failed)', .true., ok)
        return
      end if
    end if

    rc = c_chdir(trim(ctx%scratch_dir)//c_null_char)
    if (rc /= 0) then
      call skip_build(ctx, 'Real build skipped: could not chdir into scratch &
                      &directory '//trim(ctx%scratch_dir), .true., ok)
      return
    end if

    print '(a,a,a)', 'Real build scratch directory: ', trim(ctx%scratch_dir), &
      ' (removed after the build)'
  end subroutine make_scratch

  subroutine leave_build_scratch()
    !! Restore the invoking directory and remove the scratch directory
    !! make_scratch created. Called after build_and_measure returns,
    !! on both the normal path and the ibm_missing early return.
    integer(c_int) :: rc
    character(len=4096) :: expected_scratch_dir

    rc = c_chdir(trim(ctx%orig_dir)//c_null_char)
    if (rc /= 0) &
      error stop 'x3d2-memcheck: could not chdir back to the invoking &
        &directory after the real build; state is unknown, not removing &
        &the scratch directory.'

    ! Never remove anything other than the exact scratch directory
    ! make_scratch created and chdir'd into.
    expected_scratch_dir = scratch_name(ctx%orig_dir)
    if (trim(ctx%scratch_dir) == trim(expected_scratch_dir)) &
      call run_sh("rm -rf '"//trim(ctx%scratch_dir)//"'")

    ctx%orig_dir = ''
    ctx%scratch_dir = ''
  end subroutine leave_build_scratch

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

    call validate_build_inputs(build_scratch_ok)
    if (.not. build_scratch_ok) return
    call make_scratch(build_scratch_ok)
    if (.not. build_scratch_ok) return

    call build_and_measure(ctx%gdims, measured_peak_fields, used_gib, &
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
      to_gib(fields_plus_halo_bytes_n(ctx, 1, measured_peak_fields))

    call report_measured_table(measured_peak_fields, used_gib, &
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

    if (ctx%irank /= 0) return

    predicted = peak_fields_lookup(trim(ctx%domain_cfg%flow_case_name), &
                                   ctx%solver_cfg%n_species, &
                                   trim(ctx%les_cfg%model) /= 'none', &
                                   ctx%solver_cfg%ibm_on, &
                                   ctx%output_vorticity, &
                                   ctx%output_qcriterion, &
                                   ctx%solver_cfg%lowmem_transeq)
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
    !! The GPU-aware IO staging term (gpu_io_staging_bytes) is added the
    !! same way, since the real build performs no snapshot/checkpoint write
    !! and so never measures it directly.
    integer, intent(in) :: measured_peak_fields
    real(dp), intent(in) :: used_gib, workspace_gib_measured
    character(len=12), intent(out) :: requested_verdict

    integer :: k, ng, local_dims(3)
    real(dp) :: overhead, per_gpu, overhead_term, mirror_gib, w_local
    integer(i8) :: spec_bytes_1
    logical :: multi_gpu_supported_measured
    character(len=12) :: verdict
    character(len=96) :: reason

    spec_bytes_1 = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, &
                                       ctx%cdims, 1)
    overhead = used_gib - workspace_gib_measured
    requested_verdict = 'DOES_NOT_FIT'

    multi_gpu_supported_measured = .true.
    if (ctx%bc_is_010 .or. ctx%bc_is_110) then
      ! BC_y non-periodic: no nproc>1 support for this BC
      ! (src/poisson_fft.f90:178-180,196-198).
      multi_gpu_supported_measured = .false.
    else if ((ctx%bc_is_100 .or. ctx%bc_is_000) .and. &
             ctx%use_cufftmp_known .and. (.not. ctx%use_cufftmp)) then
      ! Covers the 100 case (needs cuFFTMp to decompose at all) and the
      ! fully-periodic 000 case (the plain-cuFFT fallback performs a purely
      ! local per-rank transform with no cross-rank exchange -
      ! src/backend/cuda/poisson_fft.f90:721-740 - so at nproc>1 it would
      ! run without erroring but compute invalid results). cuFFTMp
      ! unavailable here means the build just performed fell back to plain
      ! cuFFT.
      multi_gpu_supported_measured = .false.
    end if

    call print_rule('-')
    print '(a,i0,a,i0,a,f0.2,a,f0.2,a)', 'Real build (1 GPU, ', &
      ctx%measured_n_substeps, ' substep(s)): peak_fields measured ', &
      measured_peak_fields, ', workspace ', workspace_gib_measured, &
      ' GiB, used ', used_gib, ' GiB'
    ! Single machine-parseable line for --extensive sweeps (see
    ! scripts/memcheck_extensive_sweep.sh), so a wrapper script comparing
    ! several substep counts can grep one line per run instead of parsing
    ! the whole table.
    if (ctx%extensive_substeps > 0) &
      print '(a,i0,a,i0,a,f0.3,a,f0.3)', 'EXTENSIVE result: substeps=', &
        ctx%measured_n_substeps, ' peak_fields=', measured_peak_fields, &
        ' used_gib=', used_gib, ' pct_card=', 100._dp*used_gib/ctx%card_gib
    call print_table_header()

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      ! ng_unsupported_reason covers the base checks shared with report()'s
      ! table; multi_gpu_supported_measured on top of that additionally
      ! excludes ng>1 when the real build just fell back from cuFFTMp to
      ! plain cuFFT for a BC that needs it, a signal only this real-build
      ! path has (see its own derivation above).
      reason = ng_unsupported_reason(ctx, ng)
      if (ng > 1 .and. ctx%multi_gpu_supported .and. &
          .not. multi_gpu_supported_measured) &
        reason = 'not supported for this BC/environment'
      if (len_trim(reason) > 0) then
        print '(a,i0,a,a,a)', ' ', ng, '     (skipped: ', trim(reason), ')'
        cycle
      end if
      local_dims = [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)/ng]
      w_local = to_gib(fields_plus_halo_bytes_n(ctx, ng, measured_peak_fields))

      mirror_gib = 0._dp
      if (ctx%bc_is_100) &
        mirror_gib = to_gib(mirror_buffer_bytes_100(ctx%cdims, ng))

      overhead_term = overhead + &
                      to_gib(spectral_slab_bytes(ctx%bc_is_100, &
                                ctx%bc_is_110, ctx%cdims, ng) - spec_bytes_1) &
                      + mirror_gib
      if (ctx%use_cufftmp_known .and. ctx%use_cufftmp) &
        overhead_term = overhead_term + &
                        to_gib(context_floor_bytes(ng, .true.) - &
                               context_floor_bytes(1, .true.))
      overhead_term = overhead_term + &
                      to_gib(gpu_io_staging_bytes(local_dims, &
                                                  ctx%checkpoint_cfg, &
                                                  ctx%gpu_io_device_write))

      per_gpu = w_local + overhead_term
      verdict = classify(per_gpu, ctx%card_gib)
      call print_table_row(ng, local_dims, w_local, overhead_term, per_gpu, &
                           ctx%card_gib, verdict)
      if (ng == ctx%domain_cfg%nproc_dir(3)) requested_verdict = verdict
    end do
    call print_rule('-')
  end subroutine report_measured_table

  subroutine build_mesh(dims_in, mesh, dims, ibm_missing)
    !! Mesh of the real build at dims_in plus the ibm_on mask pre-check;
    !! ibm_missing=.true. tells build_and_measure to stop.
    integer, intent(in) :: dims_in(3)
    type(mesh_t), intent(out) :: mesh
    integer, intent(out) :: dims(3)
    logical, intent(out) :: ibm_missing

    character(len=16) :: ibm_file
    logical :: ibm_file_exists

    ibm_missing = .false.
    mesh = mesh_t(dims_in, [1, 1, 1], ctx%domain_cfg%L_global, &
                  ctx%domain_cfg%BC_x, ctx%domain_cfg%BC_y, &
                  ctx%domain_cfg%BC_z, ctx%domain_cfg%stretching, &
                  ctx%domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)

    ! ibm_on triggers reading an external ibm_<BC-suffix>.bp mask file
    ! inside solver init (src/module/ibm.f90), which MPI_Aborts if that
    ! file is missing. That is correct for production xcompact, but not
    ! for a tool that should degrade gracefully - check for it here, using
    ! the same suffix construction as src/module/ibm.f90:70-75 (shared via
    ! ibm_mask_filename, also used by make_scratch to mirror the
    ! mask into the build scratch directory), and bail out before
    ! triggering the case/solver construction that would abort.
    if (ctx%solver_cfg%ibm_on) then
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
    !! Build the device and host allocators and the CUDA backend; the
    !! caller owns all three so the pointers cannot outlive them.
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

    select case (trim(ctx%domain_cfg%flow_case_name))
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
    !! Drive the built flow case the way run() does so the allocator reaches
    !! its high-water mark; records measured_n_substeps and the output flags.
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
    if (ctx%extensive_substeps > 0) n_substeps = ctx%extensive_substeps
    do i_substep = 1, n_substeps
      call flow_case%substep(curr, deriv, i_substep)
    end do
    ctx%measured_n_substeps = n_substeps

    ! Mirrors run()'s per-iteration postprocessing calls (base_case.f90,
    ! immediately after the sub-stage loop): keep_pressure and
    ! vorticity/Q-criterion allocations are gated HERE, not inside
    ! postprocess() above, so they must be driven separately to be
    ! captured - postprocess(0,.) alone does not reach them.
    if (flow_case%solver%keep_pressure) &
      call compute_pressure_vert(flow_case%solver)
    ctx%output_vorticity = output_field_active( &
                           flow_case%io_mgr%snapshot_mgr%config, 'vorticity')
    ctx%output_qcriterion = output_field_active( &
                            flow_case%io_mgr%snapshot_mgr%config, 'qcriterion')
    if (ctx%output_vorticity .or. ctx%output_qcriterion) &
      call compute_derived_fields(flow_case%solver, ctx%output_vorticity, &
                                  ctx%output_qcriterion)
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
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    npeak = 0
    dev_used = 0._dp
    ctx%output_vorticity = .false.
    ctx%output_qcriterion = .false.

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
      ctx%use_cufftmp = pf%use_cufftmp
      ctx%use_cufftmp_known = .true.
    end select

    call drive_case(flow_case)

    npeak = allocator%next_id
    ierr = cudaMemGetInfo(free_b, total_b)
    dev_used = to_gib(int(total_b - free_b, i8))
  end subroutine build_and_measure

end program x3d2_memcheck
