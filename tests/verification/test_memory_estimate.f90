program test_memory_estimate
  use iso_fortran_env, only: stderr => error_unit
  use mpi
  use m_common, only: dp, i8, DIR_X, VERT, nbytes
  use m_config, only: domain_config_t
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_field, only: flist_t
  use m_solver, only: solver_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice, cudaGetDevice

  implicit none

  !> Fallback work-field count, used only for the pre-flight fit check and for
  !> grids too large to build (where the exact value is moot - they won't fit).
  !> The reported workspace uses the value MEASURED during calibration.
  integer, parameter :: PEAK_FIELDS_GUESS = 20
  integer, parameter :: n_gpu_list(4) = [1, 2, 4, 8]

  logical :: allpass = .true.
  integer :: ierr, irank, nproc, ndevs, devnum
  integer :: gdims(3), calib_lo, calib_hi, peak_fields
  character(len=256) :: input_path
  type(domain_config_t) :: domain_cfg
  logical :: calib_ok
  real(dp) :: kk, alpha, fft_hi
  real(dp) :: gpu_gib, calib_need_gib
  real(dp) :: cal_million_cells(2), cal_workspace(2), cal_dev(2)

  call initialise_mpi()
  call select_device()
  call read_config()
  call calibrate()
  call report()
  call check_result()
  call finalise()

contains

  subroutine initialise_mpi()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    ! Single-GPU calibration (enforced): the cubes keep their full local size,
    ! giving the most representative cuFFT slope. Multi-rank would split them
    ! into smaller local slabs and under-predict.
    if (nproc /= 1) error stop &
      'Run single-rank: mpirun -n 1 (the sweep extrapolates to 1/2/4/8 GPUs).'
  end subroutine initialise_mpi

  subroutine select_device()
    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(irank, ndevs))
    ierr = cudaGetDevice(devnum)
  end subroutine select_device

  subroutine read_config()
    if (command_argument_count() < 1) &
      error stop 'usage: test_memory_estimate <input.x3d>'
    call get_command_argument(1, input_path)

    call domain_cfg%read(nml_file=trim(input_path))
    gdims = domain_cfg%dims_global

    ! Calibration cubes: a quarter and a half of the target side (floored to a
    ! multiple of SZ). The half-side cube is the closest small grid to the
    ! target - best slope. Assumes a roughly cubic target.
    calib_lo = (gdims(1)/4/SZ)*SZ
    calib_hi = (gdims(1)/2/SZ)*SZ
    if (calib_lo < 64) error stop &
      'Global grid too small to calibrate (dims/4 < 64). Use dims >= 256.'
  end subroutine read_config

  subroutine calibrate()
    integer(kind=cuda_count_kind) :: free0_b, total0_b

    ! Will the (dims/2) calibration cube fit on this GPU? Its workspace is a
    ! lower bound and cuFFT needs headroom, so require it under ~70% of the
    ! card. If not, skip calibration (building it would OOM) and report
    ! workspace only.
    ierr = cudaMemGetInfo(free0_b, total0_b)
    gpu_gib = real(total0_b, dp)/1024._dp**3
    calib_need_gib = real(int(PEAK_FIELDS_GUESS, i8)*int(calib_hi, i8)**3 &
                          *int(nbytes, i8), dp)/1024._dp**3
    calib_ok = calib_need_gib < 0.7_dp*gpu_gib

    kk = 0._dp; alpha = 0._dp; peak_fields = PEAK_FIELDS_GUESS
    if (.not. calib_ok) return

    ! Smaller cube first so device(hi) - device(lo) cancels it cleanly. Each
    ! call runs one substep, which both reaches the workspace peak and measures
    ! the peak field count (the same for either cube). Delta fit: K + alpha*x.
    call measure_device(calib_lo, peak_fields, cal_million_cells(1), &
                        cal_workspace(1), cal_dev(1), verbose=.false.)
    call measure_device(calib_hi, peak_fields, cal_million_cells(2), &
                        cal_workspace(2), cal_dev(2), verbose=.true.)
    fft_hi = (cal_dev(2) - cal_dev(1)) - cal_workspace(2)
    alpha = fft_hi/cal_million_cells(2)
    kk = (cal_dev(1) - cal_workspace(1)) - alpha*cal_million_cells(1)

    if (irank == 0) then
      print '(a)', ''
      print '(a)', 'cuFFT/context fit (slope from the hi cube, intercept from &
        &the lo cube):'
      print '(a,f0.3,a,f0.3,a,f0.3)', &
        '  fft_hi [GiB] = (device_hi - device_lo) - workspace_hi = ', &
        cal_dev(2), ' - ', cal_dev(1), ' - ', cal_workspace(2)
      print '(a,f0.3,a)', '              = ', fft_hi, ' GiB'
      print '(a)', '  alpha = fft_hi / local_million_cells_hi   <- SLOPE: GiB &
        &of cuFFT/context per million local cells'
      print '(a,f0.3,a,f0.3,a,f0.6,a)', '        = ', fft_hi, ' / ', &
        cal_million_cells(2), ' = ', alpha, ' GiB/Mcell'
      print '(a)', '  kk = (device_lo - workspace_lo) - alpha*cells_lo   <- &
        &INTERCEPT: fixed cuFFT/context cost'
      print '(a,f0.3,a)', '     = ', kk, ' GiB'
      print '(a)', ''
    end if
  end subroutine calibrate

  subroutine report()
    integer :: k, ng, nx_p, ny_p, local_dims(3)
    integer(i8) :: fb
    real(dp) :: local_million_cells, workspace_gib, oh_gib, total_gib

    if (irank /= 0) return

    print '(a)', '==================================================='
    print '(a,i0,a,i0,a,i0)', 'Per-GPU memory, global grid ', &
      gdims(1), ' x ', gdims(2), ' x ', gdims(3)
    print '(a,f0.1,a)', 'GPU memory: ', gpu_gib, ' GiB total'
    if (calib_ok) then
      print '(a,i0,a)', 'Peak work fields: ', peak_fields, '  (measured)'
      print '(a,i0,a,i0,a)', 'Calibrated cuFFT/context (single-GPU, grids ', &
        calib_lo, '^3 & ', calib_hi, '^3):'
      print '(a,f0.3,a,f0.3,a)', '  overhead[GiB] ~ ', kk, ' + ', alpha, &
        ' * local_million_cells'
    else
      print '(a,i0,a)', 'Peak work fields: ', peak_fields, '  (estimate)'
      print '(a,i0,a)', 'cuFFT/context NOT calibrated: the ', calib_hi, &
        '^3 calibration grid will not fit on this GPU. Workspace only below.'
    end if
    print '(a)', '---------------------------------------------------'

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      print '(a)', ''
      print '(a,i0,a)', '[', ng, ' GPU(s)]'
      if (mod(gdims(3), ng) /= 0) then
        print '(a,i0,a)', '  skipped: global z (', gdims(3), &
          ') not divisible by GPU count'
        cycle
      end if
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]

      ! Padded field size in i8 (the allocator's ngrid is int32 and overflows
      ! for very large grids - e.g. 4096^3 = 2^36 wraps to 0). x and y are
      ! padded up to a multiple of SZ, z is not; matches allocator_init.
      nx_p = local_dims(1) - 1 + mod(-(local_dims(1) - 1), SZ) + SZ
      ny_p = local_dims(2) - 1 + mod(-(local_dims(2) - 1), SZ) + SZ
      fb = int(nx_p, i8)*int(ny_p, i8)*int(local_dims(3), i8)*int(nbytes, i8)
      workspace_gib = real(int(peak_fields, i8)*fb, dp)/1024._dp**3

      ! Unpadded local cell count (millions) - the x variable for the fit, and
      ! the SAME definition measure_device uses (dims**3/1e6), so the two
      ! match. NOT the padded byte-size used for workspace_gib above.
      local_million_cells = real(local_dims(1), dp)*real(local_dims(2), dp) &
                            *real(local_dims(3), dp)/1.0e6_dp

      if (calib_ok) then
        oh_gib = max(kk + alpha*local_million_cells, 0._dp)
      else
        oh_gib = 0._dp
      end if
      total_gib = workspace_gib + oh_gib

      print '(a,i0,a,i0,a,i0)', '  local grid:    ', local_dims(1), ' x ', &
        local_dims(2), ' x ', local_dims(3)
      print '(a,f0.3,a)', '  local cells:   ', local_million_cells, ' M'
      print '(a,f0.2,a)', '  workspace:     ', workspace_gib, ' GiB'
      print '(a,f0.2,a)', '  cuFFT/context: ', oh_gib, ' GiB'
      print '(a,f0.2,a)', '  total:         ', total_gib, ' GiB'
      if (total_gib > gpu_gib) print '(a)', '  ** EXCEEDS GPU memory **'
    end do
    print '(a)', ''

    ! Minimum GPUs of this size for the workspace alone (workspace total is
    ! decomposition-independent; cuFFT/context is extra on top).
    total_gib = real(int(peak_fields, i8)*int(gdims(1), i8)*int(gdims(2), i8) &
                     *int(gdims(3), i8)*int(nbytes, i8), dp)/1024._dp**3
    print '(a,i0,a)', 'Workspace alone needs >= ', ceiling(total_gib/gpu_gib), &
      ' GPU(s) of this size.'
    print '(a)', '==================================================='
  end subroutine report

  subroutine measure_device(dims, npeak, local_million_cells, workspace_gib, &
                            device_gib, verbose)
    !! Build the solver at dims^3 (input's config, single GPU), run one substep
    !! to reach the workspace peak, and return: npeak (peak work-field count),
    !! local Mcells, the per-GPU workspace (GiB) and the ABSOLUTE per-GPU device
    !! usage (GiB). No teardown: the build is left resident for the caller's
    !! delta. dims is small (<= dims/2 that already passed the fit check), so
    !! the allocator's int32 ngrid does not overflow here.
    integer, intent(in) :: dims
    integer, intent(out) :: npeak
    real(dp), intent(out) :: local_million_cells, workspace_gib, device_gib
    logical, intent(in), optional :: verbose

    type(mesh_t), target :: mesh
    type(cuda_allocator_t), target :: cuda_allocator
    type(allocator_t), target :: host_allocator
    class(allocator_t), pointer :: allocator
    type(cuda_backend_t), target :: cuda_backend
    class(base_backend_t), pointer :: backend
    type(solver_t) :: solver
    type(flist_t), allocatable :: curr(:), deriv(:)

    integer :: dims(3), i
    integer :: id0, id1, id2, id3, id4
    integer(kind=cuda_count_kind) :: free_b, total_b
    integer(i8) :: used_b
    logical :: verb

    verb = .false.
    if (present(verbose)) verb = verbose

    mesh = mesh_t([dims, dims, dims], [1, 1, 1], domain_cfg%L_global, &
                  domain_cfg%BC_x, domain_cfg%BC_y, domain_cfg%BC_z, &
                  domain_cfg%stretching, domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)
    cuda_allocator = cuda_allocator_t(dims, SZ)
    allocator => cuda_allocator
    host_allocator = allocator_t(dims, SZ)
    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend

    solver = solver_t(backend, mesh, host_allocator)

    call solver%u%fill(0._dp)
    call solver%v%fill(0._dp)
    call solver%w%fill(0._dp)

    allocate (curr(solver%nvars))
    curr(1)%ptr => solver%u
    curr(2)%ptr => solver%v
    curr(3)%ptr => solver%w
    do i = 1, solver%nspecies
      curr(3 + i)%ptr => solver%species(i)%ptr
    end do
    allocate (deriv(solver%nvars))

    ! One substep: drives the workspace to its high-water mark and exercises the
    ! same allocation path a real time step does. Field values are irrelevant -
    ! the block count is set by control flow, not the data. next_id is the
    ! allocator high-water mark, so the readings below are CUMULATIVE.
    id0 = allocator%next_id
    do i = 1, solver%nvars
      deriv(i)%ptr => allocator%get_block(DIR_X)
    end do
    id1 = allocator%next_id
    call solver%transeq(deriv, curr)
    id2 = allocator%next_id
    call solver%time_integrator%step(curr, deriv, solver%dt)
    id3 = allocator%next_id
    do i = 1, solver%nvars
      call allocator%release_block(deriv(i)%ptr)
    end do
    call solver%pressure_correction(solver%u, solver%v, solver%w)
    id4 = allocator%next_id

    npeak = allocator%next_id

    if (verb .and. irank == 0) then
      print '(a)', ''
      print '(a,i0,a)', 'Peak work-field count (', dims, '^3 cube).'
      print '(a)', 'next_id is the allocator high-water mark: it grows only &
        &when no freed block can'
      print '(a)', 'be reused, so these are CUMULATIVE peaks after each stage &
        &in call order - NOT'
      print '(a)', 'independent per-stage costs (e.g. pressure correction &
        &reuses the just-released'
      print '(a)', 'derivative buffers, so its line need not grow):'
      print '(a,i0)', '  after build (persistent set):       ', id0
      print '(a,i0)', '  after derivative buffers:           ', id1
      print '(a,i0)', '  after transeq:                      ', id2
      print '(a,i0)', '  after time integration:             ', id3
      print '(a,i0)', '  after pressure correction (= peak): ', id4
      print '(a)', ''
    end if

    ierr = cudaMemGetInfo(free_b, total_b)
    used_b = int(total_b - free_b, i8)

    ! Unpadded local cells (millions): the cube is dims^3 and (nproc==1) sits
    ! entirely on one GPU. Same definition as the report loop, so the fit is
    ! consistent.
    local_million_cells = real(dims, dp)**3/1.0e6_dp
    workspace_gib = real(int(npeak, i8)*int(allocator%ngrid, i8) &
                         *int(nbytes, i8), dp)/1024._dp**3
    device_gib = real(used_b, dp)/1024._dp**3
  end subroutine measure_device

  subroutine check_result()
    !! Pass/fail on the robust invariants only. The fit quality is reported as
    !! a warning, not a failure: at the tiny calibration cubes cuFFTMp allocates
    !! in chunks, so the slope can be noisy without anything being wrong.
    if (irank /= 0) return

    if (peak_fields <= 0) then
      allpass = .false.
      write (stderr, '(a)') 'Check: peak work-field count is not positive.'
    end if
    if (gpu_gib <= 0._dp) then
      allpass = .false.
      write (stderr, '(a)') 'Check: GPU memory query returned non-positive &
        &total.'
    end if

    if (calib_ok .and. .not. (alpha > 0._dp .and. kk >= 0._dp)) &
      write (stderr, '(a)') 'Warning: calibration fit looks non-physical &
        &(alpha <= 0 or kk < 0); treat the estimate with caution.'
  end subroutine check_result

  subroutine finalise()
    if (allpass) then
      if (irank == 0) write (stderr, '(a)') 'Memory estimate completed.'
    else
      error stop 'Memory estimate sanity check FAILED.'
    end if

    call MPI_Finalize(ierr)
  end subroutine finalise

end program test_memory_estimate
