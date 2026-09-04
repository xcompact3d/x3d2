program test_memory_estimate
  !! Single-grid device-memory probe (CUDA backend).
  !!
  !! Builds the solver at the grid given in the input file on ONE GPU, runs one
  !! substep to reach the allocator's work-field high-water mark, and reports
  !! the device memory in use. It does NOT fit or extrapolate: each run is one
  !! clean, independent data point (fresh process => fresh context). The
  !! calibration sweep and the safest-fit comparison are done by the driver
  !! (compare_fits.sh + fit_memory.py), which runs this on several grids.
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
  use m_common, only: dp, i8, DIR_X, VERT, nbytes
  use m_config, only: domain_config_t
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_field, only: flist_t
  use m_solver, only: solver_t
  use m_test_utils, only: padded_cells
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ

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
  real(dp) :: total_gib, used_gib, workspace_gib, mcells
  logical :: fits

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
    gdims = domain_cfg%dims_global
  end subroutine read_config

  subroutine measure_point()
    real(dp) :: ws_guess

    call query_total_gib(total_gib)

    workspace_gib = real(int(PEAK_FIELDS_GUESS, i8)*padded_cells(gdims, SZ) &
                         *int(nbytes, i8), dp)/1024._dp**3
    mcells = real(gdims(1), dp)*real(gdims(2), dp)*real(gdims(3), dp)/1.0e6_dp

    ws_guess = workspace_gib
    fits = ws_guess < BUILD_FRACTION*total_gib

    peak_fields = PEAK_FIELDS_GUESS
    used_gib = 0._dp
    if (fits) then
      call build_grid(gdims, peak_fields, used_gib)
      ! Recompute workspace with the MEASURED field count.
      workspace_gib = real(int(peak_fields, i8)*padded_cells(gdims, SZ) &
                           *int(nbytes, i8), dp)/1024._dp**3
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
    !! exact arithmetic; the overhead (fixed context + FFT scratch) is taken
    !! from the measured full-grid build, which is exact at ng=1 and a safe
    !! upper bound for ng>1 (a smaller local grid cannot need more FFT scratch
    !! than the full grid). If the input was too big to build, only the exact
    !! workspace floor is shown and the overhead is flagged as not measured.
    integer :: k, ng, local_dims(3)
    real(dp) :: overhead, w_local, per_gpu

    if (irank /= 0) return

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
    print '(a)', '-----------------------------------------------------------'
    print '(a)', ' GPUs   local grid            workspace   +overhead =  per-GPU'

    do k = 1, size(n_gpu_list)
      ng = n_gpu_list(k)
      if (mod(gdims(3), ng) /= 0) then
        print '(a,i0,a,i0,a)', ' ', ng, '     (skipped: z=', gdims(3), &
          ' not divisible)'
        cycle
      end if
      local_dims = [gdims(1), gdims(2), gdims(3)/ng]
      w_local = real(int(peak_fields, i8)*padded_cells(local_dims, SZ) &
                     *int(nbytes, i8), dp)/1024._dp**3
      per_gpu = w_local + overhead
      print '(a,i0,a,i0,a,i0,a,i0,a,f8.2,a,f7.2,a,f8.2,a)', ' ', ng, &
        '     ', local_dims(1), 'x', local_dims(2), 'x', local_dims(3), &
        '   ', w_local, '   +', overhead, '  = ', per_gpu, ' GiB'
      if (per_gpu > total_gib) print '(a)', '        ** EXCEEDS card memory **'
    end do
    print '(a)', '-----------------------------------------------------------'
  end subroutine report_per_gpu

  subroutine build_grid(dims_in, npeak, dev_used)
    !! Build the solver at dims_in on this single GPU, run one substep to drive
    !! the allocator to its work-field high-water mark, and return that mark
    !! (npeak) plus the absolute device memory in use (dev_used).
    integer, intent(in) :: dims_in(3)
    integer, intent(out) :: npeak
    real(dp), intent(out) :: dev_used

    type(mesh_t), target :: mesh
    class(allocator_t), pointer :: allocator
    class(base_backend_t), pointer :: backend
    type(solver_t) :: solver
    type(flist_t), allocatable :: curr(:), deriv(:)
    integer :: dims(3), i
    type(cuda_allocator_t), target :: cuda_allocator
    type(cuda_backend_t), target :: cuda_backend
    type(allocator_t), target :: host_allocator

    mesh = mesh_t(dims_in, [1, 1, 1], domain_cfg%L_global, &
                  domain_cfg%BC_x, domain_cfg%BC_y, domain_cfg%BC_z, &
                  domain_cfg%stretching, domain_cfg%beta, use_2decomp=.false.)
    dims = mesh%get_dims(VERT)

    cuda_allocator = cuda_allocator_t(dims, SZ)
    allocator => cuda_allocator
    host_allocator = allocator_t(dims, SZ)
    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend
    solver = solver_t(backend, mesh, host_allocator)

    ! get_block leaves data_loc unset; every case's initial_conditions() sets
    ! it to VERT before the first transeq (e.g. case/tgv.f90).
    call solver%u%set_data_loc(VERT)
    call solver%v%set_data_loc(VERT)
    call solver%w%set_data_loc(VERT)
    call solver%u%fill(0._dp)
    call solver%v%fill(0._dp)
    call solver%w%fill(0._dp)

    allocate (curr(solver%nvars))
    curr(1)%ptr => solver%u
    curr(2)%ptr => solver%v
    curr(3)%ptr => solver%w
    do i = 1, solver%nspecies
      curr(3 + i)%ptr => solver%species(i)%ptr
      call solver%species(i)%ptr%set_data_loc(VERT)
      call solver%species(i)%ptr%fill(0._dp)
    end do
    allocate (deriv(solver%nvars))

    ! One substep drives the workspace to its high-water mark via the same
    ! allocation path a real time step takes; field values are irrelevant.
    do i = 1, solver%nvars
      deriv(i)%ptr => allocator%get_block(DIR_X)
    end do
    call solver%transeq(deriv, curr)
    call solver%time_integrator%step(curr, deriv, solver%dt)
    do i = 1, solver%nvars
      call allocator%release_block(deriv(i)%ptr)
    end do
    call solver%pressure_correction(solver%u, solver%v, solver%w)

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
