module m_memcheck_estimate
  !! The per-GPU memory estimate of x3d2-memcheck.
  use m_common, only: dp, i8, nbytes
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memory_estimate, only: padded_cells, padded_halo_bytes, &
                               spectral_slab_bytes, mirror_buffer_bytes_100, &
                               spectral_extra_bytes_110, &
                               stretched_y_matrix_bytes, gpu_io_staging_bytes
  implicit none

  private
  public :: to_gib, classify, num_gpus_unsupported_reason, &
            requested_unsupported_reason, estimate_for_num_gpus, &
            prepare_probe, fields_plus_halo_bytes_n, n_halo, basis_measured, &
            basis_extrapolated, basis_formula, basis_name

  !> <80% of card memory: FITS. 80-95%: BORDERLINE. >95%: DOES_NOT_FIT.
  real(dp), parameter :: FITS_FRACTION = 0.80_dp
  real(dp), parameter :: DOES_NOT_FIT_FRACTION = 0.95_dp
  integer, parameter :: n_halo = 4

  !> How the overhead of an estimate was obtained: from the num_gpus=1 probe or
  !> measurement, that scaled to num_gpus>1, or from the formula alone.
  integer, parameter :: basis_measured = 1, basis_extrapolated = 2, &
                        basis_formula = 3
  character(len=12), parameter :: basis_name(3) = [character(len=12) :: &
                                                   'measured', &
                                                   'extrapolated', &
                                                   'formula-only']

contains

  pure function to_gib(n_bytes) result(gib)
    !! Bytes to GiB, the one conversion every printed figure goes through.
    integer(i8), intent(in) :: n_bytes
    real(dp) :: gib

    gib = real(n_bytes, dp)/1024._dp**3
  end function to_gib

  pure function classify(per_gpu_gib, total_gib) result(verdict)
    !! The FITS, BORDERLINE or DOES_NOT_FIT verdict of every estimate.
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

  function num_gpus_unsupported_reason(ctx, num_gpus) result(reason)
    !! The reason GPU count num_gpus is not a supported decomposition, blank
    !! when it is. The FFT check applies once the FFT path is known.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: num_gpus
    character(len=96) :: reason

    reason = ''
    if (num_gpus < 1) then
      reason = 'nproc_dir(3) must be >= 1'
    else if (num_gpus > 1 .and. .not. ctx%multi_gpu_supported) then
      reason = 'not supported: BC_y is non-periodic - &
               &src/poisson_fft.f90 error-stops at nproc>1'
    else if (num_gpus > 1 .and. (ctx%bc_is_000 .or. ctx%bc_is_100) .and. &
             ctx%device%fft_path_known() .and. &
             (.not. ctx%device%uses_distributed_fft())) then
      ! Plain cuFFT would give a wrong result at num_gpus>1 for 000 and 100.
      reason = 'not supported: the FFT fell back to plain cuFFT, no &
               &multi-GPU transform'
    else if (mod(ctx%gdims(3), num_gpus) /= 0) then
      write (reason, '(a,i0,a)') 'z=', ctx%gdims(3), ' not divisible'
    else if ((ctx%bc_is_000 .or. ctx%bc_is_010) .and. &
             mod(ctx%cdims(2), num_gpus) /= 0) then
      write (reason, '(a,i0,a)') 'y cells=', ctx%cdims(2), ' not divisible &
        &by num_gpus, the solver would truncate the spectral slab'
    else if (ctx%bc_is_100 .and. mod(ctx%cdims(1), num_gpus) /= 0) then
      write (reason, '(a,i0,a)') 'x cells=', ctx%cdims(1), ' not divisible &
        &by num_gpus, src/backend/cuda/poisson_fft.f90 error-stops'
    end if
  end function num_gpus_unsupported_reason

  function requested_unsupported_reason(ctx) result(reason)
    !! The reason the nproc_dir of the input file is not supported, blank
    !! when it is. The CUDA backend only supports nproc_dir = [1,1,num_gpus].
    type(memcheck_ctx_t), intent(in) :: ctx
    character(len=96) :: reason

    reason = num_gpus_unsupported_reason(ctx, ctx%domain_cfg%nproc_dir(3))
    if (ctx%domain_cfg%nproc_dir(1) /= 1 .or. &
        ctx%domain_cfg%nproc_dir(2) /= 1) &
      reason = 'only nproc_dir = [1,1,num_gpus] is supported by the CUDA &
        &backend'
  end function requested_unsupported_reason

  function fields_plus_halo_bytes_n(ctx, num_gpus, npeak) result(nbytes8)
    !! The fields and halo term at num_gpus GPUs for npeak fields.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: num_gpus, npeak
    integer(i8) :: nbytes8
    integer :: local_dims(3)

    local_dims = [ctx%gdims(1), ctx%gdims(2), ctx%gdims(3)/num_gpus]
    nbytes8 = int(npeak, i8)*padded_cells(local_dims, ctx%device%sz()) &
              *int(nbytes, i8) + &
              padded_halo_bytes(local_dims, ctx%device%sz(), n_halo)
  end function fields_plus_halo_bytes_n

  function fields_plus_halo_bytes(ctx, num_gpus) result(nbytes8)
    !! The fields and halo term at num_gpus GPUs for the static peak_fields.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: num_gpus
    integer(i8) :: nbytes8

    nbytes8 = fields_plus_halo_bytes_n(ctx, num_gpus, ctx%peak_fields)
  end function fields_plus_halo_bytes

  function spectral_plus_mirror_bytes(ctx, num_gpus) result(nbytes8)
    !! The spectral slab plus the BC specific extra terms: the mirror buffers
    !! of 100 at num_gpus>1, the second array and the transposed workspace of
    !! 110, and the Poisson coefficient matrices of 010 with y stretching.
    type(memcheck_ctx_t), intent(in) :: ctx
    integer, intent(in) :: num_gpus
    integer(i8) :: nbytes8

    nbytes8 = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, &
                                  num_gpus)
    if (ctx%bc_is_100) then
      nbytes8 = nbytes8 + mirror_buffer_bytes_100(ctx%cdims, num_gpus)
    else if (ctx%bc_is_110) then
      nbytes8 = nbytes8 + spectral_extra_bytes_110(ctx%cdims, num_gpus)
    else if (ctx%bc_is_010) then
      nbytes8 = nbytes8 + stretched_y_matrix_bytes( &
                ctx%bc_is_010, ctx%domain_cfg%stretching(2), &
                ctx%solver_cfg%lowmem_fft, ctx%cdims, num_gpus)
    end if
  end function spectral_plus_mirror_bytes

  subroutine prepare_probe(ctx)
    !! Decides the FFT probe once per run and runs it. The small-grid probe
    !! gives the context and the heap. The full-grid probe follows when the
    !! mode allows and the card has room for the probe's own size. The
    !! decision is in ctx%probe_plan and ctx%probe_reason (see ensure_probe
    !! in the device).
    type(memcheck_ctx_t), intent(inout) :: ctx
    logical :: try_full
    character(len=32) :: status_text
    integer(i8) :: total_b, free_b, probe_bytes

    call ctx%device%ensure_probe(ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, &
                                 .false.)
    ! The full-grid probe builds the FFT plans on the whole global grid, one
    ! GPU, whatever num_gpus is. Skip it when free memory cannot hold that:
    ! seven spectral copies plus the fixed overhead. A cuFFTMp plan that does
    ! not fit aborts inside NVSHMEM and kills this process. On a 40 GB card
    ! the heap took 2.5 to 5 copies, depending on the grid size, and the
    ! xtdesc buffer one. The 984^3 probe aborted with a guess of five copies,
    ! the 1056^3 probe with four and the 1280^3 probe with two, so seven is
    ! the margin. The gate is decided once per run and does not depend on
    ! num_gpus.
    try_full = .false.
    ctx%probe_reason = ''
    if (trim(ctx%run_mode) /= 'STATIC') then
      call ctx%device%mem_info(total_b, free_b)
      probe_bytes = spectral_slab_bytes(ctx%bc_is_100, ctx%bc_is_110, &
                                        ctx%cdims, 1)
      try_full = free_b >= 7_i8*probe_bytes + &
                 ctx%device%overhead_floor_bytes(1, ctx%bc_is_110)
      if (.not. try_full) &
        ctx%probe_reason = 'no full-grid probe: free memory is under seven &
                           &spectral copies plus the overhead'
    end if
    ctx%probe_plan = merge('full ', 'small', try_full)
    call ctx%device%ensure_probe(ctx%bc_is_100, ctx%bc_is_110, ctx%cdims, &
                                 try_full)
    if (ctx%device%full_failed) then
      write (status_text, '(i0)') ctx%device%full_failed_stat
      ctx%probe_plan = 'small'
      ctx%probe_reason = 'full-grid probe returned status '// &
                         trim(status_text)
    end if
  end subroutine prepare_probe

  subroutine estimate_for_num_gpus(ctx, num_gpus, per_gpu_gib, basis, &
                                   verdict, workspace_gib, overhead_gib, &
                                   io_gib, used_gib, npeak)
    !! The per-GPU estimate at num_gpus GPUs. The formula floor stands unless
    !! the full-grid probe gave a result and the floor is not over the limit.
    !! workspace_gib, overhead_gib and io_gib (optional) split the result.
    !! used_gib and npeak (optional) are the num_gpus=1 measurement of the real
    !! build and its peak_fields: the estimate is then that measurement
    !! scaled to num_gpus, and workspace_gib holds the fields and halo alone.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer, intent(in) :: num_gpus
    real(dp), intent(out) :: per_gpu_gib
    integer, intent(out) :: basis
    character(len=12), intent(out) :: verdict
    real(dp), intent(out), optional :: workspace_gib, overhead_gib, io_gib
    real(dp), intent(in), optional :: used_gib
    integer, intent(in), optional :: npeak

    integer(i8) :: base_bytes, spec_bytes, floor_bytes, overhead_bytes, &
                   io_bytes
    real(dp) :: floor_gib, io_gib_local
    integer :: peak

    peak = ctx%peak_fields
    if (present(npeak)) peak = npeak
    base_bytes = fields_plus_halo_bytes_n(ctx, num_gpus, peak)
    spec_bytes = spectral_plus_mirror_bytes(ctx, num_gpus)
    io_bytes = gpu_io_staging_bytes([ctx%gdims(1), ctx%gdims(2), &
                                     ctx%gdims(3)/num_gpus], &
                                    ctx%checkpoint_cfg, &
                                    ctx%gpu_io_device_write)
    io_gib_local = to_gib(io_bytes)
    if (present(io_gib)) io_gib = io_gib_local
    if (present(workspace_gib)) &
      workspace_gib = to_gib(base_bytes + spec_bytes)
    if (present(used_gib)) then
      per_gpu_gib = used_gib + to_gib(base_bytes + spec_bytes + io_bytes - &
                                      fields_plus_halo_bytes_n(ctx, 1, peak) &
                                      - spectral_plus_mirror_bytes(ctx, 1))
      if (ctx%device%fft_path_known() .and. &
          ctx%device%uses_distributed_fft()) &
        per_gpu_gib = per_gpu_gib + &
                      to_gib(ctx%device%overhead_floor_bytes( &
                             num_gpus, ctx%bc_is_110) - &
                             ctx%device%overhead_floor_bytes( &
                             1, ctx%bc_is_110))
      basis = merge(basis_measured, basis_extrapolated, num_gpus == 1)
      if (present(workspace_gib)) workspace_gib = to_gib(base_bytes)
      if (present(overhead_gib)) &
        overhead_gib = per_gpu_gib - to_gib(base_bytes)
      verdict = classify(per_gpu_gib, ctx%card_gib)
      return
    end if
    ! The floor assumes cuFFTMp, which gives the larger lower bound.
    floor_bytes = base_bytes + spec_bytes + io_bytes + &
                  ctx%device%overhead_floor_bytes(num_gpus, ctx%bc_is_110)
    floor_gib = to_gib(floor_bytes)

    if (len_trim(ctx%probe_reason) > 0 .and. .not. ctx%probe_note_shown .and. &
        floor_gib <= DOES_NOT_FIT_FRACTION*ctx%card_gib) then
      print '(a)', trim(ctx%probe_reason)//'; formula-only rows follow.'
      ctx%probe_note_shown = .true.
    end if

    if (floor_gib > DOES_NOT_FIT_FRACTION*ctx%card_gib .or. &
        .not. ctx%device%probe_ran()) then
      per_gpu_gib = floor_gib
      basis = basis_formula
      if (present(overhead_gib)) &
        overhead_gib = to_gib(ctx%device%overhead_floor_bytes( &
                              num_gpus, ctx%bc_is_110)) + io_gib_local
      verdict = classify(per_gpu_gib, ctx%card_gib)
      return
    end if

    overhead_bytes = ctx%device%overhead_query_bytes(num_gpus)
    per_gpu_gib = to_gib(base_bytes + spec_bytes + overhead_bytes + io_bytes)
    basis = merge(basis_measured, basis_extrapolated, num_gpus == 1)
    if (present(overhead_gib)) &
      overhead_gib = to_gib(overhead_bytes) + io_gib_local

    verdict = classify(per_gpu_gib, ctx%card_gib)
  end subroutine estimate_for_num_gpus

end module m_memcheck_estimate
