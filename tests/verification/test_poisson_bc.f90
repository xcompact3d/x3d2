program test_poisson
  !! Poisson Solver Validation Test (self-contained, no input files)
  !!
  !! Validates the Poisson solver across the boundary condition
  !! configurations listed in the table cases (4 at the moment):
  !!   ppp : all periodic          (128 x 64 x 32)
  !!   pdp : y-dirichlet           (128 x 65 x 128)
  !!   dpp : x-dirichlet           (129 x 64 x 128)
  !!   ddp : x,y-dirichlet         (129 x 257 x 64)
  !! The code is one letter per direction in x,y,z order: d = dirichlet,
  !! n = neumann, p = periodic.
  !!
  !! For each configuration, runs the cosine tests listed in the table tests
  !! (14 at the moment, n=2,3):
  !!   COS_X, COS_Y, COS_Z, COS_XY, COS_XZ, COS_YZ, COS_XYZ
  !!
  !! Each test performs two checks:
  !!   Check 1: Poisson solution vs analytical (relative L2 error)
  !!   Check 2: div(grad(p)) recovers original RHS f (round-trip relative L2
  !!            error)
  !!
  !! Each case then performs one more check, with a smooth p that is not a
  !! cosine mode:
  !!   Check 3: discrete round trip. The RHS is div(grad(p)) with the discrete
  !!            operators, and the solution must give p back up to a constant
  !!            (relative L2 norm). This holds on a stretched mesh too.
  !!
  !! Layout: setup_case builds one poisson_ctx_t for each configuration
  !! (mesh, backend, allocators, tdsops, vector_calculus, poisson_fft) and
  !! every test of that configuration reuses it. run_single_test(ctx, test)
  !! runs one row of the table tests and returns a test_result_t. It uploads
  !! the right hand side with upload_cosine, then Check 1 is solution_error
  !! and Check 2 is divgrad_error. A new check is one function plus one line
  !! in run_single_test.
  !!
  !! Self-contained: bypasses solver_t entirely; sets up mesh, backend,
  !! allocator, tdsops, vector_calculus and poisson_fft directly.
  !! No input files required.
  !!
  !! Command-line arguments, both optional:
  !!   --bc <code>                  run only the case with that BC code
  !!                                (default: every case)
  !!   --nproc_dir <px>,<py>,<pz>   decomposition, as nproc_dir in the
  !!                                domain_settings namelist
  !!                                (default: 1,1,nproc)
  !!   --stretch <mode>             run only the cases with that stretching in
  !!                                y: uniform, centred, top-bottom or bottom
  !!                                (default: uniform)
  !!   --beta <real>                stretching parameter beta of the stretched
  !!                                cases (default: the beta of the case)
  !!   --dump-matrix                rank 0 writes the stretching matrices of
  !!                                each stretched case to x3d2_<stretch>.dump
  !!                                (the format of the Incompact3d dumps)
  !! The stretched cases need the CUDA backend and 1 rank for BC code pdp. The
  !! other BC codes with stretching in y, and every stretched case on the
  !! OpenMP backend, end in the abort of the solver. On a stretched case Check
  !! 1 is printed and not judged, because cos(n*pi*y) is not an eigenfunction
  !! of the stretched discrete operator.
  !! A run with one case per process is what a multi rank ctest entry needs,
  !! because with the OpenMP backend 2decomp can be initialised only once per
  !! process and some cases still stop on more than one rank. What the solver
  !! accepts is left to its own checks. Grids stay per case, so keep at least
  !! 32 cells per rank in a split direction (see the warning in
  !! src/tdsops.f90); the dpp case (129 x 64 x 128) carries a y split on 2
  !! ranks, not on 4.
  !!
  !! Adding a BC scenario is one row in cases, and a new run is one
  !! define_verification_test line in tests/CMakeLists.txt.
  !!
  !! NOTE: Dirichlet directions require odd dims_global (e.g. 65)

  use iso_fortran_env, only: stderr => error_unit
  use ieee_arithmetic, only: ieee_is_nan

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_C, DIR_X, DIR_Y, DIR_Z, &
                      CELL, RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_Bcast
  use m_tdsops_setup, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_vector_calculus, only: vector_calculus_t
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error

  implicit none

  ! One BC configuration: the grid and the BC of each direction (the same at
  ! both ends of a direction)
  type :: bc_case_t
    character(len=30) :: title
    integer :: dims(3)
    character(len=9) :: bc(3)
    character(len=10) :: stretching ! in y: uniform|centred|top-bottom|bottom
    real(dp) :: beta ! stretching parameter, not used when uniform
  end type bc_case_t

  ! beta of examples/channel/input.x3d
  real(dp), parameter :: stretch_beta = 0.259065151_dp

  type(bc_case_t), parameter :: cases(*) = [ &
    bc_case_t('all periodic', [128, 64, 32], &
              [character(len=9) :: 'periodic', 'periodic', 'periodic'], &
              'uniform', 0.0_dp), &
    ! z carries the decomposition, so it needs enough cells per subdomain
    ! for the distributed compact operators. At 32 cells over 2 ranks (16
    ! per rank) the relative div(grad(p)) error is too large, so 128 keeps
    ! 64 cells per rank at 2 ranks, as for the x-dirichlet case below.
    bc_case_t('y-dirichlet', [128, 65, 128], &
              [character(len=9) :: 'periodic', 'dirichlet', 'periodic'], &
              'uniform', 0.0_dp), &
    ! z carries the decomposition, so it needs enough cells per subdomain
    ! for the distributed compact operators. At 32 cells over 2 ranks the
    ! discarded coupling in "interpolate" is 6.2e-08, which puts a floor
    ! under the relative div(grad(p)) error. 128 keeps 64 cells per rank at
    ! 2 ranks, matching the >=64 rule of thumb used elsewhere.
    bc_case_t('x-dirichlet', [129, 64, 128], &
              [character(len=9) :: 'dirichlet', 'periodic', 'periodic'], &
              'uniform', 0.0_dp), &
    bc_case_t('x,y-dirichlet', [129, 257, 64], &
              [character(len=9) :: 'dirichlet', 'dirichlet', 'periodic'], &
              'uniform', 0.0_dp), &
    ! Stretched in y. Only pdp is solved with the stretched operator, the
    ! other BC codes end in the abort of the solver. 64 cells in x and z.
    bc_case_t('y-dirichlet centred', [64, 65, 64], &
              [character(len=9) :: 'periodic', 'dirichlet', 'periodic'], &
              'centred', stretch_beta), &
    bc_case_t('y-dirichlet top-bottom', [64, 65, 64], &
              [character(len=9) :: 'periodic', 'dirichlet', 'periodic'], &
              'top-bottom', stretch_beta), &
    bc_case_t('y-dirichlet bottom', [64, 65, 64], &
              [character(len=9) :: 'periodic', 'dirichlet', 'periodic'], &
              'bottom', stretch_beta), &
    bc_case_t('all periodic centred', [64, 64, 64], &
              [character(len=9) :: 'periodic', 'periodic', 'periodic'], &
              'centred', stretch_beta), &
    bc_case_t('x-dirichlet centred', [65, 64, 64], &
              [character(len=9) :: 'dirichlet', 'periodic', 'periodic'], &
              'centred', stretch_beta), &
    bc_case_t('x,y-dirichlet centred', [65, 65, 64], &
              [character(len=9) :: 'dirichlet', 'dirichlet', 'periodic'], &
              'centred', stretch_beta)]

  ! One cosine test: the directions that carry cos(n*pi*x_d) and the wavenumber
  type :: cosine_test_t
    character(len=10) :: name
    logical :: uses(3)
    integer :: n
  end type cosine_test_t

  ! The n=2 block first, then the n=3 block, as the summary lists them
  type(cosine_test_t), parameter :: tests(*) = [ &
    cosine_test_t('COS_X', [.true., .false., .false.], 2), &
    cosine_test_t('COS_Y', [.false., .true., .false.], 2), &
    cosine_test_t('COS_Z', [.false., .false., .true.], 2), &
    cosine_test_t('COS_XY', [.true., .true., .false.], 2), &
    cosine_test_t('COS_XZ', [.true., .false., .true.], 2), &
    cosine_test_t('COS_YZ', [.false., .true., .true.], 2), &
    cosine_test_t('COS_XYZ', [.true., .true., .true.], 2), &
    cosine_test_t('COS_X', [.true., .false., .false.], 3), &
    cosine_test_t('COS_Y', [.false., .true., .false.], 3), &
    cosine_test_t('COS_Z', [.false., .false., .true.], 3), &
    cosine_test_t('COS_XY', [.true., .true., .false.], 3), &
    cosine_test_t('COS_XZ', [.true., .false., .true.], 3), &
    cosine_test_t('COS_YZ', [.false., .true., .true.], 3), &
    cosine_test_t('COS_XYZ', [.true., .true., .true.], 3)]

  ! The outcome of the checks that run once for each BC configuration. The
  ! default is what a skipped configuration keeps
  type :: case_result_t
    logical :: passed = .true.
    real(dp) :: roundtrip_err = 0.0_dp
  end type case_result_t

  ! The outcome of one cosine test on one BC configuration. The defaults are
  ! what a skipped configuration keeps: it must not register as a failure in
  ! the verdict
  type :: test_result_t
    logical :: passed = .true.
    logical :: xfail = .false.
    real(dp) :: poisson_err = 0.0_dp, divgrad_err = 0.0_dp
  end type test_result_t

  ! What one BC case sets up once and every test of that case reuses. The
  ! runtime holds pointers to its own members and the backend points at the
  ! mesh, so a context is declared target and is never copied
  type :: poisson_ctx_t
    type(backend_runtime_t) :: runtime
    class(base_backend_t), pointer :: backend => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(mesh_t) :: mesh
    type(dirps_t), pointer :: xdirps => null(), ydirps => null(), &
                              zdirps => null()
    type(vector_calculus_t) :: vector_calculus
    logical :: periodic(3)
    logical :: stretched ! stretching in y, Check 1 is not judged
  end type poisson_ctx_t

  ! ERROR_TOLERANCE bounds the relative L2 error ||p - p_ref|| / ||p_ref|| of
  ! the Poisson solution (Check 1). It must sit above the error of the
  ! passing cases and below the n=3 periodic aliasing error that the XFAIL
  ! logic relies on detecting; a looser tolerance turns XFAILs into
  ! unexpected passes. Measured over the passing cases of every BC case, at 1
  ! and 2 ranks:
  !   double precision, OpenMP backend: 4.4699e-08
  !   double precision, CUDA backend:   4.4699e-08
  !   single precision, OpenMP backend, nvfortran build: 1.9487e-06
  !   single precision, OpenMP backend, gfortran build:  1.5340e-06
  !   single precision, CUDA backend:                    2.7265e-07
  ! The XFAIL rows measure at least 8.5654e-01 in double precision and
  ! 8.5640e-01 in single precision. The tolerance is the {1,2,5}x10^k value
  ! just above 10x the largest passing measurement: a margin of about 11 in
  ! double precision and 10 in single precision. The XFAIL rows stay rejected.
#ifdef SINGLE_PREC
  real(dp), parameter :: ERROR_TOLERANCE = 2.0e-5_dp
#else
  real(dp), parameter :: ERROR_TOLERANCE = 5.0e-7_dp
#endif

  ! DIVGRAD_TOLERANCE bounds the relative L2 error of the div(grad(p)) round
  ! trip (Check 2) against the original right hand side. It goes through the
  ! staggered gradient and divergence operators, which carry their own
  ! roundoff on top of the Poisson solve. Largest value measured over the
  ! passing cases:
  !   double precision, OpenMP backend: 9.7767e-11
  !   double precision, CUDA backend:   9.3749e-13
  !   single precision, OpenMP backend, nvfortran build: 8.4237e-03
  !   single precision, OpenMP backend, gfortran build:  4.9875e-03
  !   single precision, CUDA backend:                    5.9888e-04
  ! The tolerance is the smallest {1,2,5}x10^k value at least 10x the largest
  ! of these, a margin of about 10 in double precision and 12 in single
  ! precision. It is independent of ERROR_TOLERANCE, on which the n=3
  ! aliasing (XFAIL) detection depends.
  !
  ! The OpenMP backend sits above the CUDA backend, by about 8x (gfortran
  ! build) to 14x (nvfortran build) in single precision and about 100x in
  ! double precision. The gfortran OpenMP builds were configured with
  ! 2decomp's FFT_Choice=generic. No other FFT engine was measured, so the
  ! gap is not attributed to one.
#ifdef SINGLE_PREC
  real(dp), parameter :: DIVGRAD_TOLERANCE = 1.0e-1_dp
#else
  real(dp), parameter :: DIVGRAD_TOLERANCE = 1.0e-9_dp
#endif

  ! Check 3, relative L2 norm of the solution against p
#ifdef SINGLE_PREC
  real(dp), parameter :: ROUNDTRIP_TOLERANCE = 1.0e-5_dp
#else
  real(dp), parameter :: ROUNDTRIP_TOLERANCE = 1.0e-10_dp
#endif

  integer :: nrank, nproc
  integer :: ic, idx, iarg, ios
  integer :: nproc_dir(3)
  character(len=32) :: arg
  character(len=len(arg)) :: only_bc
  character(len=len(cases(1)%stretching)) :: only_stretch
  real(dp) :: beta_override ! <= 0: use the beta of the case
  character(len=4*size(cases)) :: codes
  logical :: config_run(size(cases))
  logical :: allpass
  logical :: dump_matrix

  ! Per-config results for final summary
  type(test_result_t) :: results(size(tests), size(cases))
  type(case_result_t) :: case_results(size(cases))

  ! Initialise MPI
  call initialise_mpi(nrank, nproc)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  ! Optional --bc <code> restricts the run to the one configuration with that
  ! code: one letter per direction in x,y,z order, d = dirichlet, n = neumann,
  ! p = periodic (ppp, pdp, dpp, ddp). This is what lets a multi rank run
  ! exercise one configuration at a time. With no argument every
  ! configuration runs, as before.
  only_bc = 'all'
  only_stretch = 'uniform'
  beta_override = -1.0_dp
  dump_matrix = .false.
  nproc_dir = [1, 1, nproc]
  do iarg = 1, command_argument_count()
    call get_command_argument(iarg, arg)
    if (trim(arg) == '--bc') then
      call get_command_argument(iarg + 1, arg)
      only_bc = trim(arg)
    else if (trim(arg) == '--stretch') then
      call get_command_argument(iarg + 1, arg)
      only_stretch = trim(arg)
    else if (trim(arg) == '--beta') then
      call get_command_argument(iarg + 1, arg)
      read (arg, *, iostat=ios) beta_override
      if (ios /= 0) then
        if (nrank == 0) write (stderr, '(A)') &
          'test_poisson_bc: --beta needs a real number'
        error stop 1
      end if
    else if (trim(arg) == '--dump-matrix') then
      dump_matrix = .true.
    else if (trim(arg) == '--nproc_dir') then
      ! --nproc_dir <px>,<py>,<pz> sets the decomposition, as nproc_dir does
      ! in the domain_settings namelist; the default is all ranks along z
      call get_command_argument(iarg + 1, arg)
      read (arg, *, iostat=ios) nproc_dir
      if (ios /= 0) then
        if (nrank == 0) write (stderr, '(A)') &
          'test_poisson_bc: --nproc_dir needs three integers, <px>,<py>,<pz>'
        error stop 1
      end if
    end if
  end do

  if (product(nproc_dir) /= nproc) then
    if (nrank == 0) write (stderr, '(A,I0,A,I0,A)') &
      'test_poisson_bc: --nproc_dir asks for ', product(nproc_dir), &
      ' ranks but the run has ', nproc, ' ranks'
    error stop 1
  end if

  do ic = 1, size(cases)
    config_run(ic) = (only_bc == 'all' &
                      .or. only_bc == bc_code(cases(ic))) &
                     .and. only_stretch == cases(ic)%stretching
  end do

  ! An unknown code or stretching would otherwise run nothing and pass
  if (.not. any(config_run)) then
    codes = ''
    do ic = 1, size(cases)
      if (only_stretch == cases(ic)%stretching) &
        codes = trim(codes)//' '//bc_code(cases(ic))
    end do
    if (nrank == 0) write (stderr, '(A)') &
      'test_poisson_bc: no case for --bc '''//trim(only_bc)// &
      ''' --stretch '''//trim(only_stretch)//''', the codes are:'//trim(codes)
    error stop 1
  end if

  do ic = 1, size(cases)
    if (config_run(ic)) call run_config(ic, cases(ic), nproc_dir)
  end do

  ! ---- Grand summary ----
  if (nrank == 0) then
    write (stderr, '(A)') ''
    write (stderr, '(A)') &
      '  ======================================================='
    write (stderr, '(A)') &
      '                     GRAND SUMMARY'
    write (stderr, '(A)') &
      '  ======================================================='
    write (stderr, '(A)') ''
    write (stderr, '(2X,A22,2X,A10,A4,A16,A16,2X,A6,2X,A8,2X,A)') &
      'BC', 'Type      ', '  n ', 'Poisson rel L2', &
      'DivGrad rel L2', 'Result', 'Expected', ''
    write (stderr, '(A)') ''
    do ic = 1, size(cases)
      if (.not. config_run(ic)) cycle
      call print_config_results(ic)
      if (ic < size(cases)) write (stderr, '(A)') ''
    end do
    write (stderr, '(A)') ''
  end if

  ! Final verdict: allpass ignores expected failures (XFAIL)
  ! A test is a real failure if:
  !   - it failed AND was not expected to fail (unexpected failure)
  !   - it passed AND was expected to fail (unexpected pass / XPASS)
  allpass = .true.
  do ic = 1, size(cases)
    do idx = 1, size(tests)
      if (results(idx, ic)%passed .neqv. (.not. results(idx, ic)%xfail)) then
        allpass = .false.
      end if
    end do
  end do
  if (.not. all(case_results%passed)) allpass = .false.

  call finalise_test(allpass, nrank)

contains

  ! ================================================================
  ! Run all the cosine tests for one BC configuration
  ! ================================================================
  subroutine run_config(config_id, cfg, nproc_dir)
    integer, intent(in) :: config_id
    type(bc_case_t), intent(in) :: cfg
    integer, intent(in) :: nproc_dir(3)

    type(poisson_ctx_t), target :: ctx

    integer :: idx

    call setup_case(ctx, cfg, nproc_dir)

    if (dump_matrix .and. ctx%stretched .and. nrank == 0) &
      call dump_stretching_matrix(ctx, cfg)

    ! Run all the cosine tests
    do idx = 1, size(tests)
      results(idx, config_id) = run_single_test(ctx, tests(idx))
    end do

    case_results(config_id) = run_case_checks(ctx)

  end subroutine run_config

  ! ================================================================
  ! Write the stretching matrices of one case, rank 0 only. The indices are
  ! local, which are global on rank 0 (the supported run has one rank).
  ! Entries the solver never writes are written as they are.
  ! ================================================================
  subroutine dump_stretching_matrix(ctx, cfg)
    type(poisson_ctx_t), intent(in) :: ctx
    type(bc_case_t), intent(in) :: cfg

    integer :: unit

    associate (fft => ctx%backend%poisson_fft)
      open (newunit=unit, file='x3d2_'//trim(cfg%stretching)//'.dump', &
            status='replace', action='write')
      call write_dump_header(unit, cfg, fft%nx_spec, fft%nz_spec, &
                             ctx%mesh%geo%beta(2), ctx%mesh%geo%alpha(2), &
                             fft%ny_glob)
      if (fft%stretched_y_sym) then
        call write_dump_block(unit, 'odd', fft%a_odd_re, fft%a_odd_im)
        call write_dump_block(unit, 'even', fft%a_even_re, fft%a_even_im)
      else
        call write_dump_block(unit, 'single', fft%a_re, fft%a_im)
      end if
      close (unit)
    end associate
  end subroutine dump_stretching_matrix

  subroutine write_dump_header(unit, cfg, nx_spec, nz_spec, beta, alpha, ny)
    integer, intent(in) :: unit, nx_spec, nz_spec, ny
    type(bc_case_t), intent(in) :: cfg
    real(dp), intent(in) :: beta, alpha

    ! mode is the istret of Incompact3d: centred 1, top-bottom 2, bottom 3
    write (unit, '(a,1x,i0)') 'mode', &
      findloc([character(len=10) :: 'centred', 'top-bottom', 'bottom'], &
              cfg%stretching, 1)
    write (unit, '(a,1x,i0)') 'ny', ny
    write (unit, '(a,1x,es24.16)') 'beta', beta
    write (unit, '(a,1x,es24.16)') 'alpha', alpha
    write (unit, '(a,2(1x,i0))') 'spec_x_range', 1, nx_spec
    write (unit, '(a,2(1x,i0))') 'spec_z_range', 1, nz_spec
  end subroutine write_dump_header

  subroutine write_dump_block(unit, name, re, im)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: re(:, :, :, :), im(:, :, :, :)

    integer :: i, j, k, d

    write (unit, '(a,1x,a)') 'array', name
    write (unit, '(a,4(1x,i0))') 'lbound', lbound(re)
    write (unit, '(a,4(1x,i0))') 'shape', shape(re)
    do d = 1, size(re, 4)
      do k = 1, size(re, 3)
        do j = 1, size(re, 2)
          do i = 1, size(re, 1)
            write (unit, '(4(i0,1x),2es24.16)') i, j, k, d, &
              re(i, j, k, d), im(i, j, k, d)
          end do
        end do
      end do
    end do
  end subroutine write_dump_block

  ! ================================================================
  ! Set up the context of one BC configuration, once for all its tests
  ! ================================================================
  subroutine setup_case(ctx, cfg, nproc_dir)
    type(poisson_ctx_t), target, intent(inout) :: ctx
    type(bc_case_t), intent(in) :: cfg
    integer, intent(in) :: nproc_dir(3)

    integer :: dims_global(3)
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
    real(dp) :: L_global(3), beta(3)
    logical :: use_2decomp

    dims_global = cfg%dims
    beta = [0.0_dp, cfg%beta, 0.0_dp]
    if (beta_override > 0.0_dp .and. cfg%stretching /= 'uniform') &
      beta(2) = beta_override
    BC_x = [cfg%bc(1), cfg%bc(1)]
    BC_y = [cfg%bc(2), cfg%bc(2)]
    BC_z = [cfg%bc(3), cfg%bc(3)]

    if (nrank == 0) then
      write (stderr, '(A)') ''
      write (stderr, '(A,A)') '  === BC ', &
        bc_code(cfg)//' ('//trim(cfg%title)//')'
      write (stderr, '(4X,A,I0,A,I0,A,I0)') &
        'Grid: ', dims_global(1), ' x ', dims_global(2), &
        ' x ', dims_global(3)
      write (stderr, '(4X,A,A,A,F0.9)') &
        'Stretching in y: ', trim(cfg%stretching), ', beta = ', beta(2)
      write (stderr, '(4X,A,I0,A,I0,A,I0)') &
        'nproc_dir: ', nproc_dir(1), ', ', nproc_dir(2), ', ', nproc_dir(3)
      write (stderr, '(4X,A,A,A,A,A,A,A,A,A,A)') &
        'BC: x=[', trim(BC_x(1)), ',', trim(BC_x(2)), &
        '] y=[', trim(BC_y(1)), ',', trim(BC_y(2)), &
        '] z=[', trim(BC_z(1)), ','//trim(BC_z(2))//']'
    end if

    ! Setup domain decomposition
    L_global = [1.0_dp, 1.0_dp, 1.0_dp]

    ! Decide whether 2decomp is used
    use_2decomp = .not. backend_is_cuda

    ctx%mesh = mesh_t(dims_global, nproc_dir, L_global, &
                      BC_x, BC_y, BC_z, &
                      stretching=[character(len=10) :: 'uniform', &
                                  cfg%stretching, 'uniform'], &
                      beta=beta, use_2decomp=use_2decomp)

    call ctx%runtime%init(ctx%mesh)
    ctx%backend => ctx%runtime%backend
    ctx%host_allocator => ctx%runtime%host_allocator

    ! Setup tdsops directly (like test_fft.f90)
    allocate (ctx%xdirps, ctx%ydirps, ctx%zdirps)
    ctx%xdirps%dir = DIR_X
    ctx%ydirps%dir = DIR_Y
    ctx%zdirps%dir = DIR_Z
    call allocate_tdsops(ctx%xdirps, ctx%backend, ctx%mesh, 'compact6', &
                         'compact6', 'classic', 'compact6')
    call allocate_tdsops(ctx%ydirps, ctx%backend, ctx%mesh, 'compact6', &
                         'compact6', 'classic', 'compact6')
    call allocate_tdsops(ctx%zdirps, ctx%backend, ctx%mesh, 'compact6', &
                         'compact6', 'classic', 'compact6')

    ! Setup vector calculus and Poisson FFT directly
    ctx%vector_calculus = vector_calculus_t(ctx%backend)
    call ctx%backend%init_poisson_fft(ctx%mesh, ctx%xdirps, ctx%ydirps, &
                                      ctx%zdirps)

    ! Determine which directions are periodic
    ctx%periodic = (cfg%bc == 'periodic')
    ctx%stretched = (cfg%stretching /= 'uniform')
  end subroutine setup_case

  ! ================================================================
  ! The --bc code of a case: the first letter of the BC of x, y and z
  ! ================================================================
  pure function bc_code(c) result(code)
    type(bc_case_t), intent(in) :: c
    character(len=3) :: code

    code = c%bc(1)(1:1)//c%bc(2)(1:1)//c%bc(3)(1:1)
  end function bc_code

  ! ================================================================
  ! Print results for one config in the grand summary
  ! ================================================================
  subroutine print_config_results(config_id)
    integer, intent(in) :: config_id
    integer :: idx
    character(len=4) :: result_str, expected_str
    character(len=2) :: verdict_str
    logical :: passed, xfail

    do idx = 1, size(tests)
      passed = results(idx, config_id)%passed
      xfail = results(idx, config_id)%xfail

      result_str = merge('PASS', 'FAIL', passed)
      expected_str = merge('FAIL', 'PASS', xfail)

      if (passed .neqv. xfail) then
        verdict_str = 'OK'
      else
        verdict_str = '!!'
      end if

      write (stderr, '(2X,A22,2X,A10,I4,ES16.4,ES16.4,2X,A4,4X,A4,4X,A)') &
        bc_code(cases(config_id))//' '//cases(config_id)%stretching, &
        tests(idx)%name, tests(idx)%n, &
        results(idx, config_id)%poisson_err, &
        results(idx, config_id)%divgrad_err, &
        result_str, expected_str, trim(verdict_str)
    end do
    write (stderr, '(4X,A,ES12.4,2X,A)') &
      'Check 3 round trip, relative L2: ', &
      case_results(config_id)%roundtrip_err, &
      merge('PASS', 'FAIL', case_results(config_id)%passed)
  end subroutine print_config_results

  pure function is_expected_fail(test, periodic) result(xfail)
    !! Determine if a test is expected to fail.
    !!
    !! n=3 on even-sized periodic grids (64 cells) does not resolve
    !! cos(3*pi*x) cleanly due to aliasing. A test is XFAIL when n=3
    !! AND any direction involved in the test function uses periodic BCs.
    type(cosine_test_t), intent(in) :: test
    logical, intent(in) :: periodic(3)
    logical :: xfail

    xfail = test%n == 3 .and. any(test%uses .and. periodic)
  end function is_expected_fail

  ! ================================================================
  ! Fill a field with the cosine test function divided by a constant
  ! ================================================================
  subroutine fill_cosine_field(mesh, host_field, test, divisor)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(inout) :: host_field
    type(cosine_test_t), intent(in) :: test
    real(dp), intent(in) :: divisor

    integer :: i, j, k, d, dims(3)
    real(dp) :: coords(3), n_pi, val

    dims = mesh%get_dims(CELL)
    n_pi = real(test%n, dp)*pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = mesh%get_coordinates(i, j, k, CELL)
          val = 1.0_dp
          do d = 1, 3
            if (test%uses(d)) val = val*cos(n_pi*coords(d))
          end do
          host_field%data(i, j, k) = val/divisor
        end do
      end do
    end do
  end subroutine fill_cosine_field

  ! ================================================================
  ! Value of the global (1,1,1) point, known on every rank
  ! ================================================================
  function global_first_value(mesh, host_field) result(first_value)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(in) :: host_field
    real(dp) :: first_value
    integer :: ierr

    ! Rank 0 owns the global (1,1,1) point and broadcasts its value
    first_value = 0.0_dp
    if (mesh%par%is_root()) first_value = host_field%data(1, 1, 1)
    call MPI_Bcast(first_value, 1, MPI_X3D2_DP, 0, MPI_COMM_WORLD, ierr)
  end function global_first_value

  ! ================================================================
  ! Create the cosine test function on host and transfer it to a device field
  ! ================================================================
  subroutine upload_cosine(ctx, test, f)
    type(poisson_ctx_t), target, intent(in) :: ctx
    type(cosine_test_t), intent(in) :: test
    class(field_t), intent(inout) :: f

    class(field_t), pointer :: host_field

    host_field => ctx%host_allocator%get_block(DIR_C)
    call fill_cosine_field(ctx%mesh, host_field, test, 1.0_dp)
    call ctx%backend%set_field_data(f, host_field%data, DIR_C)
    call f%set_data_loc(CELL)
    call ctx%host_allocator%release_block(host_field)
  end subroutine upload_cosine

  ! ================================================================
  ! Compare the Poisson solution on a device field with the analytical one
  ! ================================================================
  function solution_error(ctx, test, f_device) result(err)
    type(poisson_ctx_t), target, intent(in) :: ctx
    type(cosine_test_t), intent(in) :: test
    class(field_t), intent(in) :: f_device
    real(dp) :: err

    class(field_t), pointer :: host_field, host_analytical
    integer :: dims(3)
    real(dp) :: n_pi_sq, first_value

    dims = ctx%mesh%get_dims(CELL)
    n_pi_sq = (real(test%n, dp)*pi)**2

    host_field => ctx%host_allocator%get_block(DIR_C)
    call ctx%backend%get_field_data(host_field%data, f_device)

    ! Remove arbitrary constant (Poisson solution unique up to a constant)
    first_value = global_first_value(ctx%mesh, host_field)
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) - first_value

    host_analytical => ctx%host_allocator%get_block(DIR_C)
    call fill_cosine_field(ctx%mesh, host_analytical, test, &
                           -(real(count(test%uses), dp)*n_pi_sq))

    ! Remove same constant from analytical
    first_value = global_first_value(ctx%mesh, host_analytical)
    host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) - first_value

    ! Compute pointwise difference
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

    err = relative_l2_error( &
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
          host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)))

    call ctx%host_allocator%release_block(host_analytical)
    call ctx%host_allocator%release_block(host_field)
  end function solution_error

  ! ================================================================
  ! Gradient of a cell centred field at the vertices, with the operators of
  ! the context. The input is in DIR_Z layout, the outputs in DIR_X layout
  ! ================================================================
  subroutine ctx_gradient(ctx, dpdx, dpdy, dpdz, p_z)
    type(poisson_ctx_t), target, intent(in) :: ctx
    class(field_t), intent(inout) :: dpdx, dpdy, dpdz
    class(field_t), intent(in) :: p_z

    ! gradient_p2v: pressure (cell) -> velocity (vert) gradient
    call ctx%vector_calculus%gradient_c2v( &
      dpdx, dpdy, dpdz, p_z, &
      ctx%xdirps%stagder_p2v, ctx%xdirps%interpl_p2v, &
      ctx%ydirps%stagder_p2v, ctx%ydirps%interpl_p2v, &
      ctx%zdirps%stagder_p2v, ctx%zdirps%interpl_p2v &
      )
  end subroutine ctx_gradient

  ! ================================================================
  ! Divergence of a vertex field at the cell centres, with the operators of
  ! the context. The inputs are in DIR_X layout, the output in DIR_Z layout
  ! ================================================================
  subroutine ctx_divergence(ctx, div_z, u, v, w)
    type(poisson_ctx_t), target, intent(in) :: ctx
    class(field_t), intent(inout) :: div_z
    class(field_t), intent(in) :: u, v, w

    ! divergence_v2p: velocity (vert) -> cell divergence
    call ctx%vector_calculus%divergence_v2c( &
      div_z, u, v, w, &
      ctx%xdirps%stagder_v2p, ctx%xdirps%interpl_v2p, &
      ctx%ydirps%stagder_v2p, ctx%ydirps%interpl_v2p, &
      ctx%zdirps%stagder_v2p, ctx%zdirps%interpl_v2p &
      )
  end subroutine ctx_divergence

  ! ================================================================
  ! div(grad(p)) of a cell centred field in DIR_C layout. The result is a new
  ! block in DIR_X layout, which the caller releases
  ! ================================================================
  subroutine divgrad_apply(ctx, f_device, f_out)
    type(poisson_ctx_t), target, intent(in) :: ctx
    class(field_t), intent(in) :: f_device
    class(field_t), pointer, intent(out) :: f_out

    class(field_t), pointer :: f_result
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input

    gradient_input => ctx%backend%allocator%get_block(DIR_Z)
    call ctx%backend%reorder(gradient_input, f_device, RDR_C2Z)

    dpdx => ctx%backend%allocator%get_block(DIR_X)
    dpdy => ctx%backend%allocator%get_block(DIR_X)
    dpdz => ctx%backend%allocator%get_block(DIR_X)

    call ctx_gradient(ctx, dpdx, dpdy, dpdz, gradient_input)
    call ctx%backend%allocator%release_block(gradient_input)

    f_result => ctx%backend%allocator%get_block(DIR_Z)
    call ctx_divergence(ctx, f_result, dpdx, dpdy, dpdz)

    call ctx%backend%allocator%release_block(dpdx)
    call ctx%backend%allocator%release_block(dpdy)
    call ctx%backend%allocator%release_block(dpdz)

    f_out => ctx%backend%allocator%get_block(DIR_X)
    call ctx%backend%reorder(f_out, f_result, RDR_Z2X)
    call ctx%backend%allocator%release_block(f_result)
  end subroutine divgrad_apply

  ! ================================================================
  ! Recover the RHS from div(grad(p)) and return its L2 error. Only reads the
  ! fields it is given; the caller keeps them and releases them
  ! ================================================================
  function divgrad_error(ctx, f_device, f_reference) result(err)
    type(poisson_ctx_t), target, intent(in) :: ctx
    class(field_t), intent(in) :: f_device, f_reference
    real(dp) :: err

    integer :: dims(3)
    class(field_t), pointer :: f_diff, host_field, host_ref

    call divgrad_apply(ctx, f_device, f_diff)

    ! Compute error: div(grad(p)) - f_original
    call ctx%backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_diff)

    host_field => ctx%host_allocator%get_block(DIR_C)
    call ctx%backend%get_field_data(host_field%data, f_diff)

    host_ref => ctx%host_allocator%get_block(DIR_C)
    call ctx%backend%get_field_data(host_ref%data, f_reference)

    dims = ctx%mesh%get_dims(CELL)
    err = relative_l2_error( &
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
          host_ref%data(1:dims(1), 1:dims(2), 1:dims(3)))

    ! Cleanup
    call ctx%backend%allocator%release_block(f_diff)
    call ctx%host_allocator%release_block(host_ref)
    call ctx%host_allocator%release_block(host_field)
  end function divgrad_error

  ! ================================================================
  ! The smooth p of Checks 3 and 5 sampled at the cell centres:
  ! exp(cos(ky*y)) * cos(kx*x) * cos(kz*z), with k = 2 pi along a periodic
  ! direction and k = pi along a wall, so that dp/dn = 0 at the walls
  ! ================================================================
  subroutine fill_smooth_p(ctx, host_field)
    type(poisson_ctx_t), target, intent(in) :: ctx
    class(field_t), intent(inout) :: host_field

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), kw(3)

    dims = ctx%mesh%get_dims(CELL)
    kw = pi
    where (ctx%periodic) kw = 2.0_dp*pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = ctx%mesh%get_coordinates(i, j, k, CELL)
          host_field%data(i, j, k) = cos(kw(1)*coords(1)) &
                                     *exp(cos(kw(2)*coords(2))) &
                                     *cos(kw(3)*coords(3))
        end do
      end do
    end do
  end subroutine fill_smooth_p

  ! ================================================================
  ! Relative L2 norm of (solution - reference), both up to a constant: the
  ! global (1,1,1) value is removed from each. Overwrites host_sol
  ! ================================================================
  function roundtrip_rel_l2(mesh, host_sol, host_ref) result(err)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(inout) :: host_sol
    class(field_t), intent(in) :: host_ref
    real(dp) :: err

    integer :: dims(3)
    real(dp) :: first_sol, first_ref

    dims = mesh%get_dims(CELL)
    first_sol = global_first_value(mesh, host_sol)
    first_ref = global_first_value(mesh, host_ref)
    host_sol%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      (host_sol%data(1:dims(1), 1:dims(2), 1:dims(3)) - first_sol) &
      - (host_ref%data(1:dims(1), 1:dims(2), 1:dims(3)) - first_ref)

    err = relative_l2_error( &
          host_sol%data(1:dims(1), 1:dims(2), 1:dims(3)), &
          host_ref%data(1:dims(1), 1:dims(2), 1:dims(3)))
  end function roundtrip_rel_l2

  ! ================================================================
  ! Check 3, the discrete round trip. The RHS is div(grad(p)) with the
  ! operators of ctx and the solve uses the solver of solve_ctx, which is ctx
  ! itself unless a caller wants another solver on the same RHS
  ! ================================================================
  function roundtrip_error(ctx, solve_ctx) result(err)
    type(poisson_ctx_t), target, intent(in) :: ctx, solve_ctx
    real(dp) :: err

    class(field_t), pointer :: p_host, p_device, rhs_x, rhs_host
    class(field_t), pointer :: f_solve, temp, sol_host

    p_host => ctx%host_allocator%get_block(DIR_C)
    call fill_smooth_p(ctx, p_host)

    p_device => ctx%backend%allocator%get_block(DIR_C, CELL)
    call ctx%backend%set_field_data(p_device, p_host%data, DIR_C)
    call p_device%set_data_loc(CELL)
    call divgrad_apply(ctx, p_device, rhs_x)
    call ctx%backend%allocator%release_block(p_device)

    rhs_host => ctx%host_allocator%get_block(DIR_C)
    call ctx%backend%get_field_data(rhs_host%data, rhs_x)
    call ctx%backend%allocator%release_block(rhs_x)

    f_solve => solve_ctx%backend%allocator%get_block(DIR_C, CELL)
    call solve_ctx%backend%set_field_data(f_solve, rhs_host%data, DIR_C)
    call f_solve%set_data_loc(CELL)
    call ctx%host_allocator%release_block(rhs_host)

    temp => solve_ctx%backend%allocator%get_block(DIR_C)
    call solve_ctx%backend%poisson_fft%solve_poisson(f_solve, temp)
    call solve_ctx%backend%allocator%release_block(temp)

    sol_host => solve_ctx%host_allocator%get_block(DIR_C)
    call solve_ctx%backend%get_field_data(sol_host%data, f_solve)
    call solve_ctx%backend%allocator%release_block(f_solve)

    err = roundtrip_rel_l2(ctx%mesh, sol_host, p_host)

    call solve_ctx%host_allocator%release_block(sol_host)
    call ctx%host_allocator%release_block(p_host)
  end function roundtrip_error

  ! ================================================================
  ! The checks that run once for each BC configuration
  ! ================================================================
  function run_case_checks(ctx) result(res)
    type(poisson_ctx_t), target, intent(in) :: ctx
    type(case_result_t) :: res

    res%roundtrip_err = roundtrip_error(ctx, ctx)
    res%passed = (res%roundtrip_err <= ROUNDTRIP_TOLERANCE)

    if (ctx%mesh%par%is_root()) then
      write (stderr, '(4X,A,ES12.4,2X,A)') &
        'Check 3 round trip, relative L2: ', res%roundtrip_err, &
        merge('PASS', 'FAIL', res%passed)
    end if
  end function run_case_checks

  ! ================================================================
  ! Run a single Poisson test (2 checks)
  ! ================================================================
  function run_single_test(ctx, test) result(res)
    type(poisson_ctx_t), target, intent(in) :: ctx
    type(cosine_test_t), intent(in) :: test
    type(test_result_t) :: res

    class(field_t), pointer :: f_device, f_reference
    class(field_t), pointer :: temp
    real(dp) :: poisson_error_norm, div_grad_error_norm
    logical :: poisson_passed, div_grad_passed, aliased

    if (ctx%mesh%par%is_root()) then
      write (stderr, '(4X,A,A,A,I1)') &
        'Running: ', trim(test%name), '  n = ', test%n
    end if

    ! Allocate fields
    f_device => ctx%backend%allocator%get_block(DIR_C, CELL)
    f_reference => ctx%backend%allocator%get_block(DIR_X)

    ! Create test function on host and transfer to device
    call upload_cosine(ctx, test, f_device)

    ! Store reference copy (in DIR_X layout) for div-grad check later
    call ctx%backend%reorder(f_reference, f_device, RDR_C2X)

    ! ---- Solve Poisson equation ----
    temp => ctx%backend%allocator%get_block(DIR_C)
    call ctx%backend%poisson_fft%solve_poisson(f_device, temp)
    call ctx%backend%allocator%release_block(temp)

    ! ---- Check 1: Poisson solution vs analytical ----
    poisson_error_norm = solution_error(ctx, test, f_device)

    ! Printed and not judged on a stretched mesh
    poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE) &
                     .or. (ctx%stretched &
                           .and. .not. ieee_is_nan(poisson_error_norm))

    ! ---- Check 2: div(grad(p)) vs original RHS ----
    div_grad_error_norm = divgrad_error(ctx, f_device, f_reference)

    ! The n=3 periodic tests are XFAIL for the aliasing of Check 1. On a
    ! stretched mesh Check 1 is not judged, so these tests are not judged
    ! at all
    aliased = ctx%stretched .and. is_expected_fail(test, ctx%periodic)
    div_grad_passed = (div_grad_error_norm <= DIVGRAD_TOLERANCE) .or. aliased

    call ctx%backend%allocator%release_block(f_device)
    call ctx%backend%allocator%release_block(f_reference)

    ! Report per-test result
    if (ctx%mesh%par%is_root()) then
      write (stderr, '(6X,A,ES12.4,A,A)') &
        'Poisson rel L2: ', poisson_error_norm, '  ', &
        merge('PASS', 'FAIL', poisson_passed)//merge(' (not judged)', &
                                                     '             ', &
                                                     ctx%stretched)
      write (stderr, '(6X,A,ES12.4,A,A)') &
        'DivGrad rel L2: ', div_grad_error_norm, '  ', &
        merge('PASS', 'FAIL', div_grad_passed)//merge(' (not judged)', &
                                                      '             ', aliased)
    end if

    res%passed = poisson_passed .and. div_grad_passed
    res%poisson_err = poisson_error_norm
    res%divgrad_err = div_grad_error_norm
    res%xfail = is_expected_fail(test, ctx%periodic) .and. .not. ctx%stretched

  end function run_single_test

end program test_poisson
