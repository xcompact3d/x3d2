program test_memory_estimate
  !! Single-grid device-memory probe (CUDA backend).
  !!
  !! Builds the real flow case at the grid given in the input file on ONE
  !! GPU, runs one substep to reach the allocator's work-field high-water
  !! mark, and reports the device memory in use. This ng=1 measurement does
  !! NOT fit or extrapolate: each run is one clean, independent data point
  !! (fresh process => fresh context). The calibration sweep and the
  !! safest-fit comparison are done by the driver (compare_fits.sh +
  !! fit_memory.py), which runs this on several grids and only ever parses
  !! the MEMINFO/MEMPOINT lines below.
  !!
  !! report_per_gpu, printed after the ng=1 measurement, DOES extrapolate:
  !! it computes an exact per-GPU workspace at each requested GPU count and
  !! applies an exact correction for the two terms that change with GPU
  !! count (the spectral array shrinking, and the 100 case's multi-rank
  !! mirror buffers), on top of the ng=1-measured, ng-invariant overhead.
  !! This table is human-readable commentary, not machine-parsed output.
  !!
  !! Output (machine readable, one MEMINFO + one MEMPOINT per run):
  !!   MEMINFO sz=<SZ> nbytes=<bytes> total_gib=<GPU total>
  !!   MEMPOINT dims=<nx>x<ny>x<nz> mcells=<M> peak_fields=<n> \
  !!            workspace_gib=<exact> used_gib=<measured> status=<OK|SKIP>
  !! A grid whose workspace alone would exceed BUILD_FRACTION of the card is
  !! reported status=SKIP without building, so the driver can still use it as a
  !! target to extrapolate to without risking an OOM here.
  use iso_fortran_env, only: stderr => error_unit
  use mpi
  use cudafor, only: cudaMemGetInfo, cuda_count_kind, &
                     cudaGetDeviceCount, cudaSetDevice, cudaGetDevice
  use m_common, only: dp, i8, VERT, nbytes
  use m_config, only: domain_config_t, solver_config_t, has_output_field
  use m_postprocess, only: compute_derived_fields, compute_pressure_vert
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_case_channel, only: case_channel_t
  use m_case_cylinder, only: case_cylinder_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t
  use m_field, only: flist_t
  use m_test_utils, only: padded_cells
  use m_test_memory_sizing, only: cell_dims, spectral_slab_bytes, &
                                  mirror_buffer_bytes_100
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t

  implicit none

  !> Fallback field count for the pre-build fit check (the report uses the
  !> MEASURED value). Grid-independent, so a fixed guess is fine here.
  integer, parameter :: PEAK_FIELDS_GUESS = 20
  !> Do not build a grid whose workspace alone exceeds this fraction of the
  !> card - leave headroom for context + FFT scratch so the probe never OOMs.
  real(dp), parameter :: BUILD_FRACTION = 0.55_dp
  !> GPU counts to report. The CUDA backend decomposes along z, so on ng GPUs
  !> each GPU holds the local grid [nx, ny, nz/ng].
  integer, parameter :: n_gpu_list(4) = [1, 2, 4, 8]

  logical :: allpass = .true.
  integer :: ierr, irank, nproc, ndevs, devnum, peak_fields
  integer :: gdims(3)
  character(len=256) :: input_path
  type(domain_config_t) :: domain_cfg
  type(solver_config_t) :: solver_cfg
  real(dp) :: total_gib, used_gib, workspace_gib, mcells
  logical :: fits
  !> BC classification (src/backend/cuda/poisson_fft.f90:232-239's four-way
  !> split), derived from the BC strings alone so it is available even for
  !> grids too large to build.
  logical :: bc_is_000, bc_is_010, bc_is_100, bc_is_110
  logical :: periodic_x, periodic_y, periodic_z
  !> Whether the ng=1 build actually used cuFFTMp (only known once a build
  !> has run) - the only real "cuFFTMp available" signal in this codebase.
  logical :: use_cufftmp = .false., use_cufftmp_known = .false.

  call initialise_mpi()
  call select_device()
  call read_config()
  call measure_point()
  call report_per_gpu()
  call check_result()
  call finalise()

contains

  subroutine initialise_mpi()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    ! Single rank: this measures the full input grid on one GPU. The driver
    ! feeds it the (small) calibration grids one at a time.
    if (nproc /= 1) error stop 'Run single-rank: mpirun -n 1.'
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
    call solver_cfg%read(nml_file=trim(input_path))
    gdims = domain_cfg%dims_global

    ! BC classification, mirroring src/backend/cuda/poisson_fft.f90:227-239.
    periodic_x = periodic_dir(domain_cfg%BC_x)
    periodic_y = periodic_dir(domain_cfg%BC_y)
    periodic_z = periodic_dir(domain_cfg%BC_z)
    bc_is_000 = periodic_x .and. periodic_y .and. periodic_z
    bc_is_010 = periodic_x .and. (.not. periodic_y) .and. periodic_z
    bc_is_100 = (.not. periodic_x) .and. periodic_y .and. periodic_z
    bc_is_110 = (.not. periodic_x) .and. (.not. periodic_y) .and. periodic_z
  end subroutine read_config

  logical function periodic_dir(bc_pair) result(is_periodic)
    !! A direction is periodic only if BOTH ends are 'periodic' - mirrors
    !! src/mesh.f90:62-88. Duplicated here (string comparison only, no
    !! mesh_t) so BC classification is available even for grids too large
    !! to build; if src/mesh.f90's BC keyword parsing ever changes, this
    !! must change with it.
    character(len=*), intent(in) :: bc_pair(2)

    is_periodic = trim(bc_pair(1)) == 'periodic' .and. &
                  trim(bc_pair(2)) == 'periodic'
  end function periodic_dir

  subroutine measure_point()
    real(dp) :: ws_guess
    logical :: ibm_missing

    call query_total_gib(total_gib)

    workspace_gib = real(int(PEAK_FIELDS_GUESS, i8)*padded_cells(gdims, SZ) &
                         *int(nbytes, i8), dp)/1024._dp**3
    mcells = real(gdims(1), dp)*real(gdims(2), dp)*real(gdims(3), dp)/1.0e6_dp

    ws_guess = workspace_gib
    fits = ws_guess < BUILD_FRACTION*total_gib

    peak_fields = PEAK_FIELDS_GUESS
    used_gib = 0._dp
    ibm_missing = .false.
    if (fits) then
      call build_grid(gdims, peak_fields, used_gib, ibm_missing)
      if (ibm_missing) then
        ! Cannot estimate this point at all; treat like a SKIP below.
        fits = .false.
      else
        ! Recompute workspace with the MEASURED field count.
        workspace_gib = real(int(peak_fields, i8)*padded_cells(gdims, SZ) &
                             *int(nbytes, i8), dp)/1024._dp**3
      end if
    end if

    if (irank == 0) then
      print '(a,i0,a,i0,a,f0.3)', 'MEMINFO sz=', SZ, ' nbytes=', nbytes, &
        ' total_gib=', total_gib
      if (fits) then
        print '(a,i0,a,i0,a,i0,a,f0.3,a,i0,a,f0.3,a,f0.3,a)', &
          'MEMPOINT dims=', gdims(1), 'x', gdims(2), 'x', gdims(3), &
          ' mcells=', mcells, ' peak_fields=', peak_fields, &
          ' workspace_gib=', workspace_gib, ' used_gib=', used_gib, &
          ' status=OK'
        print '(a,i0,a,i0,a,i0,a,f0.2,a,f0.2,a,f0.2)', &
          'Measured ', gdims(1), 'x', gdims(2), 'x', gdims(3), &
          ':  workspace=', workspace_gib, '  used=', used_gib, &
          '  overhead=', used_gib - workspace_gib
      else if (ibm_missing) then
        print '(a,i0,a,i0,a,i0,a,f0.3,a,i0,a,f0.3,a)', &
          'MEMPOINT dims=', gdims(1), 'x', gdims(2), 'x', gdims(3), &
          ' mcells=', mcells, ' peak_fields=', peak_fields, &
          ' workspace_gib=', workspace_gib, &
          ' used_gib=0.000 status=SKIP'
        print '(a)', &
          'Cannot estimate: this input has ibm_on=T but the matching &
          &ibm_<BC-suffix>.bp mask file was not found in the working &
          &directory; reported as a SKIP target only.'
      else
        print '(a,i0,a,i0,a,i0,a,f0.3,a,i0,a,f0.3,a)', &
          'MEMPOINT dims=', gdims(1), 'x', gdims(2), 'x', gdims(3), &
          ' mcells=', mcells, ' peak_fields=', peak_fields, &
          ' workspace_gib=', workspace_gib, &
          ' used_gib=0.000 status=SKIP'
        print '(a,f0.2,a,f0.1,a)', &
          'Grid too large to build here (workspace ', workspace_gib, &
          ' GiB > build limit); reported as a SKIP target only.'
      end if
    end if
  end subroutine measure_point

  subroutine report_per_gpu()
    !! Print the per-GPU memory for each GPU count. On ng GPUs the CUDA backend
    !! holds the local grid [nx, ny, nz/ng] per device. The field workspace is
    !! exact arithmetic; the overhead (context + the ng-invariant parts of
    !! the FFT scratch: y-stretching arrays, the 110 case's extra real
    !! workspace) is taken from the measured full-grid build and is exact at
    !! every ng, since none of those terms change with ng. Two terms DO
    !! change with ng and are corrected exactly on top of that overhead: the
    !! spectral array shrinking as ng grows (000/100 cases), and the 100
    !! case's multi-rank-only mirror/exchange buffers, which a single-GPU
    !! build never allocates (see m_test_memory_sizing). If the input was
    !! too big to build, only the exact workspace floor is shown and the
    !! overhead is flagged as not measured.
    !!
    !! ng>1 rows are gated on what the solver can actually run: BC_y
    !! non-periodic (010/110) hard error-stops at nproc>1
    !! (src/poisson_fft.f90:178-180,196-198) for every grid size, so no
    !! per-GPU table beyond 1 GPU is printed for those. The 100 case needs
    !! cuFFTMp specifically (src/backend/cuda/poisson_fft.f90:459-464); this
    !! is inferred from the ng=1 build, which is the only place cuFFTMp
    !! availability is ever settled in this codebase (there is no static
    !! capability query). All ng>1 rows also assume nproc_dir=[1,1,ng] (the
    !! Z-only pencil the solver forces everywhere today) - printed
    !! explicitly so this fails loudly, not silently, if 2D decomposition
    !! is ever added.
    integer :: k, ng, local_dims(3), cdims(3)
    real(dp) :: overhead, w_local, per_gpu, spec_delta_gib, mirror_gib
    integer(i8) :: spec_bytes_1
    logical :: multi_gpu_supported
    character(len=5) :: gpu_status

    if (irank /= 0) return

    ! Global CELL dims, and the ng=1 spectral-slab size the ng=1 overhead
    ! already includes - used below to correct for the two terms that
    ! actually change with ng: the spectral array shrinking (000/100), and
    ! the 100 case's multi-rank-only mirror buffers (zero at ng=1, so this
    ! is invisible to a single-GPU build). Every other overhead component
    ! (context, y-stretching, the 110 case's extra real workspace) does not
    ! vary with ng and is already exact in the ng=1 measurement.
    cdims = cell_dims(gdims, [periodic_x, periodic_y, periodic_z])
    spec_bytes_1 = spectral_slab_bytes(bc_is_100, bc_is_110, cdims, 1)

    if (fits) then
      overhead = used_gib - workspace_gib   ! measured context + FFT
    else
      overhead = 0._dp                      ! not measured; workspace floor only
    end if

    print '(a)', '-----------------------------------------------------------'
    print '(a,i0,a,i0,a,i0,a,f0.1,a)', 'Per-GPU memory for ', gdims(1), 'x', &
      gdims(2), 'x', gdims(3), '  (card ', total_gib, ' GiB)'
    if (fits) then
      print '(a,f0.2,a)', 'Overhead (context + FFT) measured at full grid: ', &
        overhead, ' GiB  (exact at 1 GPU, upper bound for more)'
    else
      print '(a)', 'Input too big to build here: showing the EXACT workspace &
        &floor only (overhead NOT included - run a grid that fits to measure it).'
    end if

    multi_gpu_supported = .true.
    if (bc_is_010 .or. bc_is_110) then
      multi_gpu_supported = .false.
      print '(a)', 'BC_y is non-periodic: the solver does not support &
        &nproc>1 for this BC combination (src/poisson_fft.f90:178-180,&
        &196-198). No per-GPU table beyond 1 GPU can be produced, because &
        &the target configuration cannot run there.'
    else if (bc_is_100 .or. bc_is_000) then
      ! Covers the 100 case (needs cuFFTMp to decompose at all) and the
      ! fully-periodic 000 case (the plain-cuFFT fallback performs a purely
      ! local per-rank transform with no cross-rank exchange -
      ! src/backend/cuda/poisson_fft.f90:721-740 - so at nproc>1 it would
      ! run without erroring but compute invalid results).
      if (use_cufftmp_known .and. (.not. use_cufftmp)) then
        multi_gpu_supported = .false.
        print '(a)', 'cuFFTMp unavailable in this environment (the ng=1 &
          &build fell back to plain cuFFT): the plain-cuFFT fallback &
          &either cannot decompose the transform (100 case, &
          &src/backend/cuda/poisson_fft.f90:460-464) or would silently &
          &compute invalid results (000 case) at nproc>1. No per-GPU &
          &table beyond 1 GPU is produced.'
      else if (.not. use_cufftmp_known) then
        print '(a)', 'cuFFTMp availability not verified (grid was not &
          &built): the ng>1 rows below assume it is available.'
      end if
    end if
    if (multi_gpu_supported) then
      print '(a)', 'Assuming nproc_dir=[1,1,ng] (Z-only pencil); no 2D/Y-Z &
        &decomposition exists in this solver today.'
    end if

    print '(a)', '-----------------------------------------------------------'
    print '(a)', ' GPUs   local grid            workspace   +overhead =  per-GPU'

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      if (ng > 1 .and. (.not. multi_gpu_supported)) then
        print '(a,i0,a)', ' ', ng, '     (skipped: not supported for this &
          &BC/environment - see note above)'
        cycle
      end if
      if (mod(gdims(3), ng) /= 0) then
        print '(a,i0,a,i0,a)', ' ', ng, '     (skipped: z=', gdims(3), &
          ' not divisible)'
        cycle
      end if
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]
      w_local = real(int(peak_fields, i8)*padded_cells(local_dims, SZ) &
                     *int(nbytes, i8), dp)/1024._dp**3

      ! Both terms are exactly zero at ng=1, so the ng=1 row is unchanged.
      spec_delta_gib = real(spectral_slab_bytes(bc_is_100, bc_is_110, &
                                                cdims, ng) - spec_bytes_1, &
                            dp)/1024._dp**3
      mirror_gib = 0._dp
      if (bc_is_100) &
        mirror_gib = real(mirror_buffer_bytes_100(cdims, ng), dp) &
                    /1024._dp**3

      per_gpu = w_local + overhead + spec_delta_gib + mirror_gib
      print '(a,i0,a,i0,a,i0,a,i0,a,f8.2,a,f7.2,a,f8.2,a)', ' ', ng, &
        '     ', local_dims(1), 'x', local_dims(2), 'x', local_dims(3), &
        '   ', w_local, '   +', overhead + spec_delta_gib + mirror_gib, &
        '  = ', per_gpu, ' GiB'
      if (per_gpu > total_gib) print '(a)', '        ** EXCEEDS card memory **'

      ! Machine-readable counterpart of the row above, one per printed ng -
      ! absence of a MEMGPU line for a given ng IS the "not runnable there"
      ! signal, so driver scripts do not need a separate flag to parse.
      ! status=OK means overhead was measured (fits); status=FLOOR means the
      ! full grid was too big to build here, so per_gpu_gib is the exact
      ! workspace floor only, a lower bound, not a full estimate.
      gpu_status = 'FLOOR'
      if (fits) gpu_status = 'OK'
      print '(a,i0,a,i0,a,i0,a,i0,a,f0.3,a,a)', &
        'MEMGPU ng=', ng, ' local_dims=', local_dims(1), 'x', &
        local_dims(2), 'x', local_dims(3), ' per_gpu_gib=', per_gpu, &
        ' status=', trim(gpu_status)
    end do
    print '(a)', '-----------------------------------------------------------'
  end subroutine report_per_gpu

  subroutine build_grid(dims_in, npeak, dev_used, ibm_missing)
    !! Build the real flow case at dims_in on this single GPU (the same
    !! dispatch xcompact.f90 uses), run one substep to drive the allocator
    !! to its work-field high-water mark, and return that mark (npeak) plus
    !! the absolute device memory in use (dev_used).
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
    !! case_channel_init and friends re-read get_argument(1) (the process's
    !! own CLI argument) rather than taking domain_cfg as a parameter; this
    !! already resolves to the same input file read_config() parsed above,
    !! but would silently diverge if the probe's CLI contract ever gained a
    !! second positional argument.
    !!
    !! case_init also runs the input's own restart logic (io_mgr%is_restart)
    !! before this subroutine gets control: if a checkpoint file matching
    !! the input's restart path already exists in the CWD, the probe will
    !! measure a restarted state instead of the input's initial conditions.
    !! This mirrors production behaviour, not a probe-specific choice.
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
    logical :: output_vorticity, output_qcriterion

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
    ! inside solver init (src/module/ibm.f90), which MPI_Aborts if that file
    ! is missing. That is correct for production xcompact, but not for an
    ! estimator that should degrade gracefully - check for it here, using
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

    ! Same dispatch as xcompact.f90:111-126.
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
      error stop 'test_memory_estimate: unknown flow_case_name.'
    end select

    ! Solver construction (inside case_init, above) already ran
    ! init_poisson_fft, so the plan's actual cuFFTMp/cuFFT fallback outcome
    ! is settled - capture it for report_per_gpu's safety gating on the 100
    ! case, since there is no static "is cuFFTMp available" query.
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
    if (flow_case%solver%keep_pressure) call compute_pressure_vert(flow_case%solver)
    if (has_output_field(flow_case%io_mgr%snapshot_mgr%config, 'vorticity') &
        .and. flow_case%io_mgr%snapshot_mgr%config%snapshot_freq > 0) &
      output_vorticity = .true.
    if (has_output_field(flow_case%io_mgr%snapshot_mgr%config, 'qcriterion') &
        .and. flow_case%io_mgr%snapshot_mgr%config%snapshot_freq > 0) &
      output_qcriterion = .true.
    if (output_vorticity .or. output_qcriterion) &
      call compute_derived_fields(flow_case%solver, output_vorticity, &
                                  output_qcriterion)

    npeak = allocator%next_id
    call query_used_gib(dev_used)
  end subroutine build_grid

  subroutine query_total_gib(total)
    real(dp), intent(out) :: total
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    total = real(total_b, dp)/1024._dp**3
  end subroutine query_total_gib

  subroutine query_used_gib(used)
    real(dp), intent(out) :: used
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    used = real(total_b - free_b, dp)/1024._dp**3
  end subroutine query_used_gib

  subroutine check_result()
    if (irank /= 0) return

    if (total_gib <= 0._dp) then
      allpass = .false.
      write (stderr, '(a)') 'Check: GPU memory query returned non-positive &
        &total.'
    end if
    if (fits .and. peak_fields <= 0) then
      allpass = .false.
      write (stderr, '(a)') 'Check: peak work-field count is not positive.'
    end if
  end subroutine check_result

  subroutine finalise()
    if (allpass) then
      if (irank == 0) print '(a)', 'Memory estimate: PASSED'
    else
      error stop 'Memory estimate: FAILED (sanity check).'
    end if

    call MPI_Finalize(ierr)
  end subroutine finalise

end program test_memory_estimate
