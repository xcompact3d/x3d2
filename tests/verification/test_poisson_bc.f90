program test_poisson
  !! Poisson Solver Validation Test (self-contained, no input files)
  !!
  !! Validates the Poisson solver across 4 boundary condition configurations:
  !!   Config 000 : all periodic          (128 x 64 x 32)
  !!   Config 010 : y-dirichlet           (128 x 65 x 32)
  !!   Config 100 : x-dirichlet           (129 x 64 x 128)
  !!   Config 110 : x,y-dirichlet         (129 x 257 x 64)
  !!
  !! For each configuration, runs 8 cosine test cases (n=2,3):
  !!   COS_X, COS_Y, COS_XY, COS_XYZ
  !!
  !! Each test performs two checks:
  !!   Check 1: Poisson solution vs analytical (L2 norm)
  !!   Check 2: div(grad(p)) recovers original RHS f (round-trip L2 norm)
  !!
  !! Total: 4 configs x 8 cases = 32 tests
  !!
  !! Self-contained: bypasses solver_t entirely; sets up mesh, backend,
  !! allocator, tdsops, vector_calculus and poisson_fft directly.
  !! No input files or command-line arguments required.
  !!
  !! NOTE: Dirichlet directions require odd dims_global (e.g. 65)

  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda
  use m_common, only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, CELL, &
                      RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_vector_calculus, only: vector_calculus_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  ! One BC configuration: the grid and the BC of each direction (the same at
  ! both ends of a direction)
  type :: bc_case_t
    character(len=20) :: title
    integer :: dims(3)
    character(len=9) :: bc(3)
  end type bc_case_t

  type(bc_case_t), parameter :: cases(*) = [ &
    bc_case_t('all periodic', [128, 64, 32], &
              [character(len=9) :: 'periodic', 'periodic', 'periodic']), &
    bc_case_t('y-dirichlet', [128, 65, 32], &
              [character(len=9) :: 'periodic', 'dirichlet', 'periodic']), &
    ! z carries the decomposition, so it needs enough cells per subdomain
    ! for the distributed compact operators. At 32 cells over 2 ranks the
    ! discarded coupling in "interpolate" is 6.2e-08, which puts a 8.7e-11
    ! floor under the div(grad(p)) check against a 1e-11 tolerance. 128
    ! keeps 64 cells per rank at 2 ranks, matching the >=64 rule of thumb
    ! used elsewhere.
    bc_case_t('x-dirichlet', [129, 64, 128], &
              [character(len=9) :: 'dirichlet', 'periodic', 'periodic']), &
    bc_case_t('x,y-dirichlet', [129, 257, 64], &
              [character(len=9) :: 'dirichlet', 'dirichlet', 'periodic'])]

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
    cosine_test_t('COS_XY', [.true., .true., .false.], 2), &
    cosine_test_t('COS_XYZ', [.true., .true., .true.], 2), &
    cosine_test_t('COS_X', [.true., .false., .false.], 3), &
    cosine_test_t('COS_Y', [.false., .true., .false.], 3), &
    cosine_test_t('COS_XY', [.true., .true., .false.], 3), &
    cosine_test_t('COS_XYZ', [.true., .true., .true.], 3)]

  ! The outcome of one cosine test on one BC configuration. The defaults are
  ! what a skipped configuration keeps: it must not register as a failure in
  ! the verdict
  type :: test_result_t
    logical :: passed = .true.
    logical :: xfail = .false.
    real(dp) :: poisson_err = 0.0_dp, divgrad_err = 0.0_dp
  end type test_result_t

  ! The single precision tolerance must sit between the roundoff floor of
  ! the passing cases (~3e-7 in this normalised norm, norm2/N) and the
  ! n=3 periodic aliasing error (~2.4e-6) that the XFAIL logic relies on
  ! detecting; a looser tolerance turns XFAILs into unexpected passes.
#ifdef SINGLE_PREC
  real(dp), parameter :: ERROR_TOLERANCE = 1.0e-6_dp
#else
  real(dp), parameter :: ERROR_TOLERANCE = 1.0e-11_dp
#endif

  ! The div(grad(p)) round trip (Check 2) goes through the staggered
  ! gradient and divergence operators, which carry their own single
  ! precision roundoff on top of the Poisson solve; on the OpenMP backend
  ! that reaches 5.0e-6 on config 110 (COS_Y n=2) and 1.2e-6 on config 100
  ! (COS_X n=2), while the CUDA backend stays below 3.5e-7 on the same
  ! grids. A looser tolerance here, about 2x the observed OpenMP maximum,
  ! keeps Check 2 meaningful without disturbing the Poisson tolerance
  ! above, which the n=3 aliasing (XFAIL) detection depends on.
  !
  ! The higher OpenMP roundoff traces to the FFT engine: this build uses
  ! 2decomp's generic FFT engine (the build does not forward FFT_Choice to
  ! the 2decomp sub-build), whose single precision roundoff on config 110
  ! is about 10x FFTW's (Check 2 measured at 2.4e-6 with the generic
  ! engine vs 2.2e-7 with FFTW, both below the tolerance above).
#ifdef SINGLE_PREC
  real(dp), parameter :: DIVGRAD_TOLERANCE = 1.0e-5_dp
#else
  real(dp), parameter :: DIVGRAD_TOLERANCE = 1.0e-11_dp
#endif

  integer :: nrank, nproc
  integer :: ic, idx, iarg
  character(len=32) :: arg
  character(len=3) :: only_bc
  logical :: config_run(size(cases))
  logical :: allpass

  ! Per-config results for final summary
  type(test_result_t) :: results(size(tests), size(cases))

  ! Initialise MPI
  call initialise_mpi(nrank, nproc)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  ! Optional --bc <code> restricts the run to the one configuration with that
  ! code: one letter per direction in x,y,z order, d = dirichlet, n = neumann,
  ! p = periodic (ppp, pdp, dpp, ddp). This is what lets a multi rank run
  ! exercise dpp without tripping the single rank stops still in place for pdp
  ! and ddp in src/poisson_fft.f90. With no argument every configuration runs,
  ! as before.
  only_bc = 'all'
  do iarg = 1, command_argument_count() - 1
    call get_command_argument(iarg, arg)
    if (trim(arg) == '--bc') then
      call get_command_argument(iarg + 1, arg)
      only_bc = trim(arg)
    end if
  end do

  do ic = 1, size(cases)
    config_run(ic) = (only_bc == 'all' &
                      .or. only_bc == bc_code(cases(ic)))
  end do

  do ic = 1, size(cases)
    if (config_run(ic)) call run_config(ic, cases(ic))
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
    write (stderr, '(2X,A5,2X,A10,A4,A14,A14,2X,A6,2X,A8,2X,A)') &
      'BC', 'Type      ', '  n ', '  Poisson L2  ', &
      '  DivGrad L2  ', 'Result', 'Expected', ''
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

  call finalise_test(allpass, nrank)

contains

  ! ================================================================
  ! Run all 8 cosine tests for one BC configuration
  ! ================================================================
  subroutine run_config(config_id, cfg)
    integer, intent(in) :: config_id
    type(bc_case_t), intent(in) :: cfg

    type(backend_runtime_t), target :: runtime
    class(base_backend_t), pointer :: backend
    type(allocator_t), pointer :: host_allocator
    type(mesh_t), target :: mesh
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    type(vector_calculus_t) :: vector_calculus

    integer :: dims_global(3)
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
    integer :: nproc_dir(3)
    real(dp) :: L_global(3)
    logical :: use_2decomp

    integer :: idx
    logical :: passed
    logical :: periodic(3)

    dims_global = cfg%dims
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
      write (stderr, '(4X,A,A,A,A,A,A,A,A,A,A)') &
        'BC: x=[', trim(BC_x(1)), ',', trim(BC_x(2)), &
        '] y=[', trim(BC_y(1)), ',', trim(BC_y(2)), &
        '] z=[', trim(BC_z(1)), ','//trim(BC_z(2))//']'
    end if

    ! Setup domain decomposition
    nproc_dir = [1, 1, nproc]
    L_global = [1.0_dp, 1.0_dp, 1.0_dp]

    ! Decide whether 2decomp is used
    use_2decomp = .not. backend_is_cuda

    mesh = mesh_t(dims_global, nproc_dir, L_global, &
                  BC_x, BC_y, BC_z, &
                  use_2decomp=use_2decomp)

    call runtime%init(mesh)
    backend => runtime%backend
    host_allocator => runtime%host_allocator

    ! Setup tdsops directly (like test_fft.f90)
    allocate (xdirps, ydirps, zdirps)
    xdirps%dir = DIR_X
    ydirps%dir = DIR_Y
    zdirps%dir = DIR_Z
    call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')
    call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')
    call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')

    ! Setup vector calculus and Poisson FFT directly
    vector_calculus = vector_calculus_t(backend)
    call backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)

    ! Determine which directions are periodic
    periodic = (cfg%bc == 'periodic')

    ! Run all 8 cosine tests
    do idx = 1, size(tests)
      results(idx, config_id)%xfail = is_expected_fail(tests(idx), periodic)

      call run_single_test(backend, host_allocator, mesh, &
                           xdirps, ydirps, zdirps, vector_calculus, &
                           tests(idx), passed, &
                           results(idx, config_id)%poisson_err, &
                           results(idx, config_id)%divgrad_err)
      results(idx, config_id)%passed = passed
    end do

  end subroutine run_config

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

      write (stderr, '(2X,A5,2X,A10,I4,ES14.4,ES14.4,2X,A4,4X,A4,4X,A)') &
        bc_code(cases(config_id)), &
        tests(idx)%name, tests(idx)%n, &
        results(idx, config_id)%poisson_err, &
        results(idx, config_id)%divgrad_err, &
        result_str, expected_str, trim(verdict_str)
    end do
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
  ! Compute normalized L2 error norm
  ! ================================================================
  function compute_error_norm(mesh, host_field) result(error_norm)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(in) :: host_field
    real(dp) :: error_norm
    integer :: dims(3)

    dims = mesh%get_dims(CELL)
    error_norm = norm2(host_field%data(1:dims(1), 1:dims(2), 1:dims(3))) &
                 /product(dims)
  end function compute_error_norm

  ! ================================================================
  ! Run a single Poisson test (2 checks)
  ! ================================================================
  subroutine run_single_test(backend, host_allocator, mesh, &
                             xdirps, ydirps, zdirps, vector_calculus, &
                             test, test_passed, &
                             poisson_err_out, divgrad_err_out)
    class(base_backend_t), pointer, intent(in) :: backend
    type(allocator_t), pointer, intent(in) :: host_allocator
    type(mesh_t), intent(in) :: mesh
    type(dirps_t), pointer, intent(in) :: xdirps, ydirps, zdirps
    type(vector_calculus_t), intent(in) :: vector_calculus
    type(cosine_test_t), intent(in) :: test
    logical, intent(out) :: test_passed
    real(dp), intent(out) :: poisson_err_out, divgrad_err_out

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, host_analytical, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3)
    real(dp) :: poisson_error_norm, div_grad_error_norm, n_pi_sq
    logical :: poisson_passed, div_grad_passed

    dims = mesh%get_dims(CELL)
    n_pi_sq = (real(test%n, dp)*pi)**2

    if (mesh%par%is_root()) then
      write (stderr, '(4X,A,A,A,I1)') &
        'Running: ', trim(test%name), '  n = ', test%n
    end if

    ! Allocate fields
    f_device => backend%allocator%get_block(DIR_C, CELL)
    f_reference => backend%allocator%get_block(DIR_X)
    host_field => host_allocator%get_block(DIR_C)

    ! Create test function on host and transfer to device
    call fill_cosine_field(mesh, host_field, test, 1.0_dp)
    call backend%set_field_data(f_device, host_field%data, DIR_C)
    call f_device%set_data_loc(CELL)
    call host_allocator%release_block(host_field)

    ! Store reference copy (in DIR_X layout) for div-grad check later
    call backend%reorder(f_reference, f_device, RDR_C2X)

    ! ---- Solve Poisson equation ----
    temp => backend%allocator%get_block(DIR_C)
    call backend%poisson_fft%solve_poisson(f_device, temp)
    call backend%allocator%release_block(temp)

    ! ---- Check 1: Poisson solution vs analytical ----
    host_field => host_allocator%get_block(DIR_C)
    call backend%get_field_data(host_field%data, f_device)

    ! Remove arbitrary constant (Poisson solution unique up to a constant)
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_field%data(1, 1, 1)

    host_analytical => host_allocator%get_block(DIR_C)
    call fill_cosine_field(mesh, host_analytical, test, &
                           -(real(count(test%uses), dp)*n_pi_sq))

    ! Remove same constant from analytical
    host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1, 1, 1)

    ! Compute pointwise difference
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

    poisson_error_norm = compute_error_norm(mesh, host_field)

    call host_allocator%release_block(host_analytical)
    call host_allocator%release_block(host_field)

    poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE)

    ! ---- Check 2: div(grad(p)) vs original RHS ----
    gradient_input => backend%allocator%get_block(DIR_Z)
    call backend%reorder(gradient_input, f_device, RDR_C2Z)
    call backend%allocator%release_block(f_device)

    dpdx => backend%allocator%get_block(DIR_X)
    dpdy => backend%allocator%get_block(DIR_X)
    dpdz => backend%allocator%get_block(DIR_X)

    ! gradient_p2v: pressure (cell) -> velocity (vert) gradient
    call vector_calculus%gradient_c2v( &
      dpdx, dpdy, dpdz, gradient_input, &
      xdirps%stagder_p2v, xdirps%interpl_p2v, &
      ydirps%stagder_p2v, ydirps%interpl_p2v, &
      zdirps%stagder_p2v, zdirps%interpl_p2v &
      )
    call backend%allocator%release_block(gradient_input)

    f_result => backend%allocator%get_block(DIR_Z)

    ! divergence_v2p: velocity (vert) -> cell divergence
    call vector_calculus%divergence_v2c( &
      f_result, dpdx, dpdy, dpdz, &
      xdirps%stagder_v2p, xdirps%interpl_v2p, &
      ydirps%stagder_v2p, ydirps%interpl_v2p, &
      zdirps%stagder_v2p, zdirps%interpl_v2p &
      )

    call backend%allocator%release_block(dpdx)
    call backend%allocator%release_block(dpdy)
    call backend%allocator%release_block(dpdz)

    f_device => backend%allocator%get_block(DIR_X)
    call backend%reorder(f_device, f_result, RDR_Z2X)
    call backend%allocator%release_block(f_result)

    ! Compute error: div(grad(p)) - f_original
    call backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

    host_field => host_allocator%get_block(DIR_C)
    call backend%get_field_data(host_field%data, f_device)
    div_grad_error_norm = compute_error_norm(mesh, host_field)

    ! Cleanup
    call backend%allocator%release_block(f_device)
    call backend%allocator%release_block(f_reference)
    call host_allocator%release_block(host_field)

    div_grad_passed = (div_grad_error_norm <= DIVGRAD_TOLERANCE)

    ! Report per-test result
    if (mesh%par%is_root()) then
      write (stderr, '(6X,A,ES12.4,A,A)') &
        'Poisson L2: ', poisson_error_norm, '  ', &
        merge('PASS', 'FAIL', poisson_passed)
      write (stderr, '(6X,A,ES12.4,A,A)') &
        'DivGrad L2: ', div_grad_error_norm, '  ', &
        merge('PASS', 'FAIL', div_grad_passed)
    end if

    test_passed = poisson_passed .and. div_grad_passed
    poisson_err_out = poisson_error_norm
    divgrad_err_out = div_grad_error_norm

  end subroutine run_single_test

end program test_poisson
