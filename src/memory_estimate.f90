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
  use m_common, only: i8, nbytes
  use m_config, only: checkpoint_config_t, has_output_field
  implicit none

  private
  public :: padded_dim, padded_cells, cell_dims, spectral_slab_bytes, &
            mirror_buffer_bytes_100, output_field_active, &
            peak_fields_lookup, halo_bytes

contains

  pure integer function padded_dim(n, sz) result(n_padded)
    !! Smallest multiple of sz that is >= n. Mirrors the x/y padding in
    !! allocator_init (m_allocator); z is never padded so is not covered here.
    integer, intent(in) :: n, sz

    n_padded = n - 1 + mod(-(n - 1), sz) + sz
  end function padded_dim

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

  pure function spectral_slab_bytes(bc_is_100, bc_is_110, cdims, ng) &
    result(nbytes8)
    !! Bytes of ONE complex spectral array (waves_dev, and c_dev on the
    !! plain-cuFFT fallback path), shape per
    !! src/backend/cuda/poisson_fft.f90:259-291:
    !!   100: (cy/2+1, cx/ng, cz)  - decomposition axis maps to dim2
    !!   110: (cz/2+1, cx, cy)    - never runs at ng>1 in this solver
    !!   000/010: (cx/2+1, cy/ng, cz) - decomposition axis maps to dim2
    !! ng is the number of ranks the z-direction is split over
    !! (mesh%par%nproc_dir(3)); note this differs from the physical
    !! decomposition axis (z), which is why ng divides cx or cy here, not
    !! cz.
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    integer, intent(in) :: ng
    integer(i8) :: nbytes8
    integer(i8) :: nx_spec, ny_spec, nz_spec

    if (bc_is_100) then
      nx_spec = cdims(2)/2 + 1
      ny_spec = cdims(1)/ng
      nz_spec = cdims(3)
    else if (bc_is_110) then
      nx_spec = cdims(3)/2 + 1
      ny_spec = cdims(1)
      nz_spec = cdims(2)
    else
      ! 000, 010
      nx_spec = cdims(1)/2 + 1
      ny_spec = cdims(2)/ng
      nz_spec = cdims(3)
    end if
    nbytes8 = nx_spec*ny_spec*nz_spec*2_i8*int(nbytes, i8)
  end function spectral_slab_bytes

  pure function mirror_buffer_bytes_100(cdims, ng) result(nbytes8)
    !! 100-case multi rank only mirror/exchange buffers, zero at ng<=1 -
    !! matches the nproc>1 guard at
    !! src/backend/cuda/poisson_fft.f90:332-349:
    !!   c_mirror_dev + c_slab_send_dev: (nx_spec, ny_spec, nz_spec) x2
    !!   c_plane_send_dev + c_plane_recv_dev: (nx_spec, nz_spec) x2
    !! all complex(dp).
    integer, intent(in) :: cdims(3)
    integer, intent(in) :: ng
    integer(i8) :: nbytes8
    integer(i8) :: nx_spec, ny_spec, nz_spec

    if (ng <= 1) then
      nbytes8 = 0_i8
      return
    end if
    nx_spec = cdims(2)/2 + 1
    ny_spec = cdims(1)/ng
    nz_spec = cdims(3)
    nbytes8 = 2_i8*nx_spec*ny_spec*nz_spec*2_i8*int(nbytes, i8) + &
              2_i8*nx_spec*nz_spec*2_i8*int(nbytes, i8)
  end function mirror_buffer_bytes_100

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

    n_groups = max(int(padded_dims(2), i8)*int(padded_dims(3), i8)/int(sz, i8), &
                   int(padded_dims(1), i8)*int(padded_dims(3), i8)/int(sz, i8), &
                   int(padded_dims(1), i8)*int(padded_dims(2), i8)/int(sz, i8))
    nbytes8 = int(sz, i8)*n_groups*int(nbytes, i8)*12_i8*int(n_halo + 1, i8)
  end function halo_bytes

  pure logical function output_field_active(checkpoint_cfg, name) result(active)
    !! Whether a snapshot output field is actually driven at runtime - a
    !! field only triggers its associated allocations if BOTH it is listed
    !! in output_fields AND snapshotting is enabled. Mirrors the pattern
    !! used for keep_pressure (src/case/base_case.f90:120-122) and for
    !! vorticity/qcriterion (src/memcheck.f90's build_and_measure),
    !! consolidated here instead of duplicated a third time.
    type(checkpoint_config_t), intent(in) :: checkpoint_cfg
    character(*), intent(in) :: name

    active = has_output_field(checkpoint_cfg, name) .and. &
             checkpoint_cfg%snapshot_freq > 0
  end function output_field_active

  pure integer function peak_fields_lookup(flow_case_name, n_species, &
                                           les_on, ibm_on, &
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
    integer, intent(in) :: n_species
    logical, intent(in) :: les_on, ibm_on
    logical, intent(in) :: output_vorticity, output_qcriterion
    logical, intent(in) :: lowmem_transeq

    integer :: base, derived_fields_delta

    ! Base count for the regular (non-lowmem) transeq path, ibm_on=F,
    ! n_species=0, les_on=F, no derived-field output. Measured 2026-09-09,
    ! NVHPC 25.3, CUDA backend, one substep:
    !   tgv      examples/TGV/input.x3d                (256^3, ng=1) -> 20
    !   generic  examples/generic/input.x3d                          -> 20
    !   channel  examples/channel/input.x3d  (128x65x64)             -> 26
    !            (channel's lazily-allocated BC ghost blocks, per
    !            define_BC_channel, account for the +6 over tgv/generic)
    !   cylinder examples/cylinder/input_ibm_100_LR-64.x3d with
    !            ibm_on=F (257x130x64)                               -> 23
    !            confirmed by the self-check (see check_peak_fields_table
    !            in src/memcheck.f90) against BOTH the 100 and 110 BC
    !            variants with ibm_on=T (measured 24 = 23+1, matching the
    !            ibm_on delta below exactly) - a built-in ADIOS2 build
    !            (deps/adios2-install-cuda-v2.12.1) was needed since
    !            cylinder's mask file is ADIOS2 .bp format
    ! foil is not a supported flow case on this branch (xcompact.f90 has no
    ! 'foil' dispatch), so the default branch below only needs to cover
    ! genuinely unknown/future flow_case_name values. Falls back to
    ! channel's base (the largest measured, i.e. safest over-estimate)
    ! rather than under-estimating with tgv/generic's base.
    select case (trim(flow_case_name))
    case ('tgv')
      base = 20
    case ('generic')
      base = 20
    case ('channel')
      base = 26
    case ('cylinder')
      base = 23
    case default
      ! Covers any unrecognised flow_case_name: see note above.
      base = 26
    end select

    ! lowmem_transeq: TODO, not yet calibrated - it dispatches to a distinct
    ! code path (transeq_lowmem, src/solver.f90:223-224,410-526) with a
    ! different get_block call pattern than the regular transeq measured
    ! above. base is not adjusted for it (unverified placeholder: assumes
    ! the regular-path base, which may over- or under-count).

    ! output_vorticity/output_qcriterion are NOT simply additive: measured
    ! 2026-09-09 on examples/TGV/input.x3d (256^3, ng=1, base case=20):
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

    ! n_species: +1 per species, exact by construction, not measured - each
    ! species gets exactly one persistent block, same pattern as u/v/w
    ! (src/solver.f90:159-165: `solver%species(i)%ptr =>
    ! backend%allocator%get_block(DIR_X)` in a loop over 1..nspecies, never
    ! released).
    !
    ! ibm_on: +1, exact by construction, not measured - ibm_t%init
    ! (src/module/ibm.f90) gets exactly one device block for ep1
    ! (get_block(DIR_X)), regardless of grid size.
    !
    ! les_on: +10, measured 2026-09-09 on
    ! examples/TGV/input_static_smagorinsky.x3d (160^3, model='smagorinsky',
    ! ng=1): base tgv=20 -> 30 with LES on. Caught by
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
    n = base + derived_fields_delta + n_species + merge(1, 0, ibm_on) + &
        merge(10, 0, les_on)
  end function peak_fields_lookup

end module m_memory_estimate
