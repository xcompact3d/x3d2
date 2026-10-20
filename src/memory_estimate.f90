module m_memory_estimate
  !! Static (no-run) GPU memory footprint formulas for the CUDA backend.
  !!
  !! Byte-exact sizes for the field, halo, and cuFFT/cuFFTMp spectral-array
  !! terms, mirroring the allocate statements in src/allocator.f90,
  !! src/backend/cuda/backend.f90, and src/backend/cuda/poisson_fft.f90.
  !! These formulas were originally written and tested inside tests/ as
  !! part of a runtime memory probe, since merged into src/memcheck.f90's
  !! --build Tier 3 path; they live here so both that real build measurement
  !! and the static, no build Tiers 1-2 estimate can share them instead of
  !! duplicating them.
  use m_allocator, only: padded_dim
  use m_common, only: i8, nbytes, is_sp, sp
  use m_config, only: checkpoint_config_t, has_output_field, &
                      output_field_active
  use m_time_integrator, only: nolds_of
  implicit none

  private
  public :: padded_dim, padded_cells, cell_dims, spectral_slab_bytes, &
            mirror_buffer_bytes_100, spectral_extra_bytes_110, &
            stretched_y_matrix_bytes, output_field_active, &
            peak_fields_lookup, halo_bytes, padded_halo_bytes, &
            gpu_io_staging_parts, gpu_io_staging_bytes, &
            predicted_fft_heap_bytes, round_up_bytes, &
            fft_xtdesc_step_bytes, fft_probe_margin_bytes

  !> The NVSHMEM heap of a cuFFTMp plan grows in whole GiB for the workspace.
  !> Measured on an A100 40 GB with 512^3, 640^3 and 1024^3 in double
  !> precision and 704^3 and 1280^3 in single precision: the plan memory is
  !> 2.775 GiB for all five grids, and plan memory plus workspace add up to
  !> the memory of a plan with auto-allocation on to 0.001 GiB. Other cards
  !> can differ.
  integer(i8), parameter :: fft_workspace_step_bytes = 1024_i8**3
  !> The xtdesc buffer grows the heap in steps of about 0.5 GiB.
  integer(i8), parameter :: fft_xtdesc_step_bytes = 512_i8*1024_i8**2
  !> The free memory the full-grid FFT probe keeps beyond the predicted heap.
  !> The probe context and the libraries that stay after cufftDestroy use it.
  integer(i8), parameter :: fft_probe_margin_bytes = 1024_i8**3

contains

  pure function padded_cells(dims, sz) result(n)
    !! Cell count of one padded field at dims (x, y padded to sz, z untouched).
    !! i8 so callers can multiply by field count/nbytes without overflow, as
    !! the allocator's own int32 ngrid does for very large grids.
    integer, intent(in) :: dims(3), sz
    integer(i8) :: n

    n = int(padded_dim(dims(1), sz), i8)*int(padded_dim(dims(2), sz), i8) &
        *int(dims(3), i8)
  end function padded_cells

  pure function cell_dims(vert_dims, periodic) result(c)
    !! Global CELL dims from global VERT dims. Mirrors src/mesh.f90:90-100:
    !! a periodic direction has as many cells as vertices, a non-periodic
    !! direction has one fewer.
    integer, intent(in) :: vert_dims(3)
    logical, intent(in) :: periodic(3)
    integer :: c(3)
    integer :: dir

    do dir = 1, 3
      if (periodic(dir)) then
        c(dir) = vert_dims(dir)
      else
        c(dir) = vert_dims(dir) - 1
      end if
    end do
  end function cell_dims

  pure function spectral_slab_bytes(n_spec) result(nbytes8)
    !! Bytes of ONE complex spectral array (waves_dev, and c_dev on the
    !! plain-cuFFT fallback path) of the shape n_spec. n_spec comes from
    !! cuda_fft_spectral_dims in src/backend/cuda/poisson_fft.f90, the routine
    !! the solver init calls, so the shape per BC class and per number of
    !! ranks is the solver's. Whether a decomposition is valid is checked by
    !! cuda_fft_decomp_reason, not here.
    integer, intent(in) :: n_spec(3)
    integer(i8) :: nbytes8

    nbytes8 = int(n_spec(1), i8)*int(n_spec(2), i8)*int(n_spec(3), i8) &
              *2_i8*int(nbytes, i8)
  end function spectral_slab_bytes

  pure function round_up_bytes(n_bytes, step_bytes) result(nbytes8)
    !! The smallest multiple of step_bytes that is at least n_bytes.
    integer(i8), intent(in) :: n_bytes, step_bytes
    integer(i8) :: nbytes8

    nbytes8 = ((n_bytes + step_bytes - 1_i8)/step_bytes)*step_bytes
  end function round_up_bytes

  pure function predicted_fft_heap_bytes(plan_bytes, worksize_bytes) &
    result(nbytes8)
    !! The NVSHMEM heap of the forward and backward cuFFTMp plans with the
    !! workspace allocated and without the xtdesc buffer: the plan memory plus
    !! the workspace rounded up to 1 GiB. This is the definition of the
    !! measured heap of the normal FFT workspace query.
    integer(i8), intent(in) :: plan_bytes, worksize_bytes
    integer(i8) :: nbytes8

    nbytes8 = plan_bytes + &
              round_up_bytes(worksize_bytes, fft_workspace_step_bytes)
  end function predicted_fft_heap_bytes

  pure function mirror_buffer_bytes_100(n_spec, num_gpus) result(nbytes8)
    !! 100-case multi rank only mirror/exchange buffers, zero at num_gpus<=1 -
    !! matches the nproc>1 guard at
    !! src/backend/cuda/poisson_fft.f90:332-349:
    !!   c_mirror_dev + c_slab_send_dev: (nx_spec, ny_spec, nz_spec) x2
    !!   c_plane_send_dev + c_plane_recv_dev: (nx_spec, nz_spec) x2
    !! all complex(dp). n_spec is the spectral shape of the 100 case.
    integer, intent(in) :: n_spec(3)
    integer, intent(in) :: num_gpus
    integer(i8) :: nbytes8

    if (num_gpus <= 1) then
      nbytes8 = 0_i8
      return
    end if
    nbytes8 = 2_i8*int(n_spec(1), i8)*int(n_spec(2), i8)*int(n_spec(3), i8) &
              *2_i8*int(nbytes, i8) + &
              2_i8*int(n_spec(1), i8)*int(n_spec(3), i8)*2_i8*int(nbytes, i8)
  end function mirror_buffer_bytes_100

  pure function spectral_extra_bytes_110(n_spec, cdims) result(nbytes8)
    !! 110-case bytes on top of the waves_dev slab (never runs at num_gpus>1 in
    !! this solver): a SECOND complex array of the same spectral shape n_spec,
    !! c_dev (src/backend/cuda/poisson_fft.f90:416-418), plus a real,
    !! globally-shaped transposed workspace r_dev_110 (:414), which has the
    !! same total element count as cdims, just permuted.
    integer, intent(in) :: n_spec(3), cdims(3)
    integer(i8) :: nbytes8

    nbytes8 = spectral_slab_bytes(n_spec) &
              + int(cdims(1), i8)*int(cdims(2), i8)*int(cdims(3), i8) &
              *int(nbytes, i8)
  end function spectral_extra_bytes_110

  pure function stretched_y_matrix_bytes(bc_is_010, stretching_y, lowmem_fft, &
                                         n_spec) result(nbytes8)
    !! 010-case, non-uniform y-stretching Poisson coefficient matrices,
    !! src/poisson_fft.f90:176-186 (only 010 sets stretched_y, and only at
    !! num_gpus=1 since 010 error-stops at nproc>1) and :320-323,424-428
    !! (allocation shapes), copied to the device at
    !! src/backend/cuda/poisson_fft.f90:390-411:
    !!   'bottom': a_re, a_im, shape (nx_spec, ny_spec, nz_spec, 5), real(dp)
    !!   otherwise ('centred', 'top-bottom'): a_odd_re/im, a_even_re/im,
    !!     shape (nx_spec, ny_spec/2, nz_spec, 5), real(dp)
    !! doubled by the store_a_*_dev copies when lowmem_fft is .false.
    !! (the default).
    logical, intent(in) :: bc_is_010, lowmem_fft
    character(len=*), intent(in) :: stretching_y
    integer, intent(in) :: n_spec(3)
    integer(i8) :: nbytes8
    integer(i8) :: nx_spec, ny_spec, nz_spec

    if ((.not. bc_is_010) .or. trim(stretching_y) == 'uniform') then
      nbytes8 = 0_i8
      return
    end if

    nx_spec = n_spec(1)
    ny_spec = n_spec(2)
    nz_spec = n_spec(3)

    if (trim(stretching_y) == 'bottom') then
      nbytes8 = 2_i8*nx_spec*ny_spec*nz_spec*5_i8*int(nbytes, i8)
    else
      nbytes8 = 4_i8*nx_spec*(ny_spec/2_i8)*nz_spec*5_i8*int(nbytes, i8)
    end if

    if (.not. lowmem_fft) nbytes8 = 2_i8*nbytes8
  end function stretched_y_matrix_bytes

  pure function halo_bytes(padded_dims, sz, n_halo) result(nbytes8)
    !! CUDA-backend halo-exchange buffer bytes. Mirrors the 24 allocate
    !! statements in src/backend/cuda/backend.f90:128-152: 12 arrays shaped
    !! (sz, n_halo, n_groups) for u/v/w send/recv start/end, plus 12 more
    !! shaped (sz, 1, n_groups) for du/dud/d2u send/recv start/end - all
    !! real(dp). n_groups is the max over X/Y/Z of ny_p*nz_p/sz, nx_p*nz_p/
    !! sz, nx_p*ny_p/sz (src/allocator.f90:81-83, allocator_t%
    !! n_groups_dir/get_n_groups), where *_p are the sz-padded dims (x/y
    !! padded via padded_dim, z untouched, same rule padded_cells uses).
    integer, intent(in) :: padded_dims(3) !! [nx_padded, ny_padded, nz_padded]
    integer, intent(in) :: sz, n_halo
    integer(i8) :: nbytes8
    integer(i8) :: n_groups

    n_groups = max( &
               int(padded_dims(2), i8)*int(padded_dims(3), i8)/int(sz, i8), &
               int(padded_dims(1), i8)*int(padded_dims(3), i8)/int(sz, i8), &
               int(padded_dims(1), i8)*int(padded_dims(2), i8)/int(sz, i8))
    nbytes8 = int(sz, i8)*n_groups*int(nbytes, i8)*12_i8*int(n_halo + 1, i8)
  end function halo_bytes

  pure function padded_halo_bytes(local_dims, sz, n_halo) result(nbytes8)
    !! halo_bytes for a local block of local_dims cells: pads x and y to sz
    !! (z untouched) the way the allocator does, then sizes the halo-exchange
    !! buffers of those padded dims.
    integer, intent(in) :: local_dims(3)
    integer, intent(in) :: sz, n_halo
    integer(i8) :: nbytes8

    nbytes8 = halo_bytes([padded_dim(local_dims(1), sz), &
                          padded_dim(local_dims(2), sz), local_dims(3)], &
                         sz, n_halo)
  end function padded_halo_bytes

  pure subroutine gpu_io_staging_parts(checkpoint_cfg, snapshot_dev, &
                                       checkpoint_dev)
    !! Which writes can take the device path: a snapshot only at unit stride
    !! (striding falls back to the host path in snapshot_manager.f90), a
    !! checkpoint never strides. Shared by gpu_io_staging_bytes and the
    !! x3d2-memcheck report.
    type(checkpoint_config_t), intent(in) :: checkpoint_cfg
    logical, intent(out) :: snapshot_dev, checkpoint_dev

    snapshot_dev = checkpoint_cfg%snapshot_freq > 0 .and. &
                   all(checkpoint_cfg%output_stride == 1)
    checkpoint_dev = checkpoint_cfg%checkpoint_freq > 0
  end subroutine gpu_io_staging_parts

  pure function gpu_io_staging_bytes(local_dims, checkpoint_cfg, &
                                     device_write) result(bytes)
    !! Device memory the GPU-aware ADIOS2 write path (src/io/adios2/io.f90,
    !! write_device_field_adios2) holds on top of the solver's own peak: a
    !! single reused staging buffer (file_handle%device_staging), sized to
    !! one unpadded LOCAL field in whichever precision that write uses.
    !! The buffer is grown, never shrunk (see write_device_field_adios2),
    !! so its steady-state size is the LARGEST single write it is asked to
    !! hold: the build precision whenever a checkpoint is enabled (always
    !! native, and always at least as large as a single-precision
    !! snapshot), otherwise whatever precision the snapshot uses.
    !!
    !! Zero when the device path cannot be taken: device_write false
    !! (build without X3D2_ADIOS2_CUDA, or X3D2_ADIOS2_GPU_WRITE_MODE=
    !! host), or neither a unit-stride snapshot nor a checkpoint is
    !! enabled (striding falls back to the host path in
    !! snapshot_manager.f90; checkpoints never stride).
    !!
    !! The DIR_X->DIR_C reorder borrows a block from the allocator free
    !! list and is not counted here (measurements 2026-09-15 showed
    !! exactly one field of growth).
    integer, intent(in) :: local_dims(3)
    type(checkpoint_config_t), intent(in) :: checkpoint_cfg
    logical, intent(in) :: device_write
    integer(i8) :: bytes

    ! Bytes of one real(sp) value, needed only when the staging buffer is
    ! sized for a single-precision snapshot rather than the build precision.
    integer(i8), parameter :: sp_bytes = int(storage_size(0.0_sp)/8, i8)
    integer(i8) :: n, element_bytes
    logical :: snapshot_dev, checkpoint_dev

    call gpu_io_staging_parts(checkpoint_cfg, snapshot_dev, checkpoint_dev)

    if (.not. device_write .or. .not. (snapshot_dev .or. checkpoint_dev)) then
      bytes = 0_i8
      return
    end if

    n = int(local_dims(1), i8)*int(local_dims(2), i8)*int(local_dims(3), i8)
    element_bytes = int(nbytes, i8)
    if (.not. checkpoint_dev .and. snapshot_dev .and. &
        checkpoint_cfg%snapshot_sp .and. .not. is_sp) element_bytes = sp_bytes
    bytes = n*element_bytes
  end function gpu_io_staging_bytes

  integer function peak_fields_lookup(flow_case_name, time_intg, &
                                      n_species, les_on, ibm_on, &
                                      output_vorticity, &
                                      output_qcriterion, &
                                      lowmem_transeq) result(n)
    !! Static estimate of the allocator's work-field high-water mark
    !! (allocator_t%next_id after one substep - see src/allocator.f90),
    !! evaluated purely from config flags, no build/run required.
    !!
    !! This is grid independent: which get_block calls a run makes depends
    !! only on flow_case_name and which features are switched on, not on
    !! the grid dimensions. Terms below are additive EXCEPT
    !! output_vorticity/output_qcriterion, which are measured jointly - see
    !! that block's comment for why.
    !!
    !! Every number here was measured with the real build+measure probe
    !! (src/memcheck.f90's --build, CUDA backend, regular -not
    !! lowmem_transeq- transeq path) unless flagged TODO, in which case it
    !! is a best-effort value derived by reading the allocation site rather
    !! than measuring it, and should be treated as provisional until
    !! confirmed. --build's self-verification check (check_peak_fields_
    !! table in src/memcheck.f90) compares its own measured next_id against
    !! this function for the same config and warns on any mismatch, so
    !! this table is expected to be corrected/extended over time rather
    !! than treated as final.
    character(*), intent(in) :: flow_case_name
    character(len=3), intent(in) :: time_intg
    integer, intent(in) :: n_species
    logical, intent(in) :: les_on, ibm_on
    logical, intent(in) :: output_vorticity, output_qcriterion
    logical, intent(in) :: lowmem_transeq

    integer :: base, derived_fields_delta, nolds

    ! Base count for the regular (non-lowmem) transeq path, ibm_on=F,
    ! n_species=0, les_on=F, no derived-field output, WITHOUT the history
    ! blocks of the time integrator (3*nolds, added below). Measured
    ! 2026-09-09 with each example's own integrator, NVHPC 25.3, CUDA
    ! backend, one substep, then re-expressed on 2026-10-13 by taking the
    ! history blocks out (nolds = 2 for AB3, 3 for RK3):
    !   tgv      examples/TGV/input.x3d (256^3, num_gpus=1, AB3) -> 20, base 14
    !   generic  examples/generic/input.x3d (AB3)            -> 20, base 14
    !   channel  examples/channel/input.x3d  (128x65x64, RK3) -> 26, base 17
    !            (channel's lazily-allocated BC ghost blocks, per
    !            define_BC_channel, account for the +3 over tgv/generic)
    !   cylinder examples/cylinder/input_ibm_100_LR-64.x3d with
    !            ibm_on=F (257x130x64, AB3)                   -> 23, base 17
    !            confirmed by the self-check (see check_peak_fields_table
    !            in src/memcheck.f90) against BOTH the 100 and 110 BC
    !            variants with ibm_on=T (measured 24 = 23+1, matching the
    !            ibm_on delta below exactly) - a built-in ADIOS2 build
    !            (deps/adios2-install-cuda-v2.12.1) was needed since
    !            cylinder's mask file is ADIOS2 .bp format
    ! Measured 2026-10-13 (PR 326 review, S1a-S1d, S1g) on the tgv input
    ! with the integrator changed: AB3 20, RK3 23, RK4 26, AB2 17, i.e.
    ! 14 + 3*nolds; and on channel with AB3: 23 = 17 + 3*2.
    ! foil is TODO: not yet calibrated, and has no case module anywhere in
    ! this tree (no m_case_foil, no case ('foil') dispatch in xcompact.f90
    ! or in memcheck.f90's build_and_measure) - `x3d2-memcheck
    ! examples/foil/input.x3d --build` cannot measure a real number until
    ! one is added, and validate_build_inputs' flow_case_supported check
    ! skips the real build for it by name rather than attempting it. Falls
    ! back here to channel's base (17, the largest calibrated base) plus the
    ! history blocks of the selected integrator, 17 + 3*nolds: 26 with RK3
    ! and 23 with AB3.  The static (Tier 1) estimate stays on the high side
    ! in the meantime.
    select case (trim(flow_case_name))
    case ('tgv')
      base = 14
    case ('generic')
      base = 14
    case ('channel')
      base = 17
    case ('cylinder')
      base = 17
    case ('abl')
      ! examples/abl/input.x3d, AB3 + smagorinsky: 14 + 6 + 7 = 27, measured 27.
      base = 14
    case default
      ! Covers 'foil' and anything else: see TODO above.
      base = 17
    end select

    ! time_intg: time_intg_t%init allocates nolds history blocks for each of
    ! the 3 velocity components and for each species, held for the whole run.
    nolds = nolds_of(time_intg)

    ! lowmem_transeq: TODO, not yet calibrated - it dispatches to a distinct
    ! code path (transeq_lowmem, src/solver.f90:223-224,410-526) with a
    ! different get_block call pattern than the regular transeq measured
    ! above. base is not adjusted for it (unverified placeholder: assumes
    ! the regular-path base, which may over- or under-count).

    ! output_vorticity/output_qcriterion are NOT simply additive: measured
    ! 2026-09-09 on examples/TGV/input.x3d (256^3, num_gpus=1, base case=20):
    !   vorticity only            -> 21 (+1)
    !   qcriterion only           -> 21 (+1)
    !   vorticity AND qcriterion  -> 24 (+4, not +2)
    ! compute_derived_fields (src/postprocess/postprocess.f90) evidently
    ! holds extra shared gradient temporaries concurrently when computing
    ! both in the same call that it does not need for either alone, so the
    ! combination is looked up directly rather than summed from the two
    ! individual deltas.
    if (output_vorticity .and. output_qcriterion) then
      derived_fields_delta = 4
    else if (output_vorticity .or. output_qcriterion) then
      derived_fields_delta = 1
    else
      derived_fields_delta = 0
    end if

    ! n_species: +(2 + nolds) per species - each species gets exactly one
    ! persistent block, same pattern as u/v/w (src/solver.f90:159-165:
    ! `solver%species(i)%ptr => backend%allocator%get_block(DIR_X)` in a
    ! loop over 1..nspecies, never released), plus its nolds history blocks
    ! (nvars = 3 + n_species in time_intg_t%init), plus one transient
    ! derivative block (base_case.f90 takes one per variable, 1..nvars, for
    ! each substep). Measured 2026-10-13 (PR 326 review, S1e) on the tgv
    ! input with AB3 and n_species = 2: 28 = 14 + 3*2 + 2*(2 + 2).
    !
    ! ibm_on: +1, exact by construction, not measured - ibm_t%init
    ! (src/module/ibm.f90) gets exactly one device block for ep1
    ! (get_block(DIR_X)), regardless of grid size.
    !
    ! les_on: +7, measured 2026-09-09 on
    ! examples/TGV/input_static_smagorinsky.x3d (160^3, model='smagorinsky',
    ! num_gpus=1): base tgv=20 -> 30 with LES on with RK3 (+10 including the
    ! history blocks; 7 once the 3 RK3 history blocks are taken out, and
    ! 2026-10-13 S1f: AB3 -> 27 = 14 + 6 + 7). Caught by
    ! check_peak_fields_table's self-verification warning (src/memcheck.f90
    ! --build): an initial guess of +2
    ! (from reading only les.f90's two persistent blocks, nut and
    ! mixing_length_sq - src/les/les.f90:168-181,222-224) undercounted by
    ! 8, because LES's own transient temporaries (dudx..dwdz etc., 9 blocks
    ! momentarily live in compute_and_set_strain_rate, src/les/les.f90:
    ! 270-285) coincide with the solver's own peak rather than being
    ! computed at a quiet point. Only measured for the tgv/generic base;
    ! not yet confirmed additive on top of channel/cylinder's BC ghost
    ! blocks or combined with ibm_on/species/derived-field output.
    !
    ! keep_pressure: +0, measured (examples/TGV/input.x3d with
    ! output_fields='pressure', snapshot_freq=1 -> still 20) - the pressure
    ! block persisting rather than being released does not raise the
    ! high-water mark for this case, so no term is added for it.
    n = base + 3*nolds + derived_fields_delta + n_species*(2 + nolds) + &
        merge(1, 0, ibm_on) + merge(7, 0, les_on)
  end function peak_fields_lookup

end module m_memory_estimate
