program test_poisson
  !! Poisson Solver Validation Test (self-contained, no input files)
  !!
  !! Validates the Poisson solver across 4 boundary condition configurations:
  !!   Config 000 : all periodic          (64 x 64 x 64)
  !!   Config 010 : y-dirichlet           (64 x 65 x 64)
  !!   Config 100 : x-dirichlet           (65 x 64 x 64)
  !!   Config 110 : x,y-dirichlet         (65 x 65 x 64)
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
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, VERT, CELL, &
                      RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_vector_calculus, only: vector_calculus_t

#ifdef CUDA
  use cudafor

  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  ! Test type identifiers
  integer, parameter :: TEST_COS_X = 1
  integer, parameter :: TEST_COS_Y = 2
  integer, parameter :: TEST_COS_XY = 3
  integer, parameter :: TEST_COS_XYZ = 4

  integer, parameter :: NUM_TYPES = 4
  integer, parameter :: NUM_NS = 2
  integer, parameter :: NUM_TESTS = NUM_TYPES*NUM_NS
  integer, parameter :: NUM_CONFIGS = 4
  integer, parameter :: TOTAL_TESTS = NUM_CONFIGS*NUM_TESTS

  real(dp), parameter :: ERROR_TOLERANCE = 1.0e-11_dp

  integer :: nrank, nproc, ierr
  integer :: ic, idx
  logical :: allpass
  character(32) :: backend_name

  ! Per-config results for final summary
  logical :: all_results(NUM_TESTS, NUM_CONFIGS)
  logical :: all_xfail(NUM_TESTS, NUM_CONFIGS)
  real(dp) :: all_poisson_errs(NUM_TESTS, NUM_CONFIGS)
  real(dp) :: all_divgrad_errs(NUM_TESTS, NUM_CONFIGS)
  character(len=3) :: config_labels(NUM_CONFIGS)

  ! BC configuration arrays
  integer :: dims_global(3)
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)

  ! Initialise MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

#ifdef CUDA
  block
    integer :: ndevs, devnum
    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(devnum)
  end block
  backend_name = "CUDA"
#else
  backend_name = "OMP"
#endif

  config_labels = ['000', '010', '100', '110']

  ! ---- Config 000: all periodic (64 x 64 x 64) ----
  dims_global = [64, 64, 64]
  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']
  call run_config(1, '000 (all periodic)', dims_global, &
                  BC_x, BC_y, BC_z)

  ! ---- Config 010: y-dirichlet (64 x 65 x 64) ----
  dims_global = [64, 65, 64]
  BC_x = ['periodic ', 'periodic ']
  BC_y = ['dirichlet', 'dirichlet']
  BC_z = ['periodic ', 'periodic ']
  call run_config(2, '010 (y-dirichlet)', dims_global, &
                  BC_x, BC_y, BC_z)

  ! ---- Config 100: x-dirichlet (65 x 64 x 64) ----
  dims_global = [65, 64, 64]
  BC_x = ['dirichlet', 'dirichlet']
  BC_y = ['periodic ', 'periodic ']
  BC_z = ['periodic ', 'periodic ']
  call run_config(3, '100 (x-dirichlet)', dims_global, &
                  BC_x, BC_y, BC_z)

  ! ---- Config 110: x,y-dirichlet (65 x 65 x 64) ----
  dims_global = [65, 65, 64]
  BC_x = ['dirichlet', 'dirichlet']
  BC_y = ['dirichlet', 'dirichlet']
  BC_z = ['periodic ', 'periodic ']
  call run_config(4, '110 (x,y-dirichlet)', dims_global, &
                  BC_x, BC_y, BC_z)

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
      'Conf', 'Type      ', '  n ', '  Poisson L2  ', &
      '  DivGrad L2  ', 'Result', 'Expected', ''
    write (stderr, '(A)') ''
    do ic = 1, NUM_CONFIGS
      call print_config_results(ic)
      if (ic < NUM_CONFIGS) write (stderr, '(A)') ''
    end do
    write (stderr, '(A)') ''
  end if

  ! Final verdict: allpass ignores expected failures (XFAIL)
  ! A test is a real failure if:
  !   - it failed AND was not expected to fail (unexpected failure)
  !   - it passed AND was expected to fail (unexpected pass / XPASS)
  allpass = .true.
  do ic = 1, NUM_CONFIGS
    do idx = 1, NUM_TESTS
      if (all_results(idx, ic) .neqv. (.not. all_xfail(idx, ic))) then
        allpass = .false.
      end if
    end do
  end do

  if (allpass) then
    if (nrank == 0) then
      write (stderr, '(A)') 'ALL TESTS PASSED SUCCESSFULLY.'
    end if
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

contains

  ! ================================================================
  ! Run all 8 cosine tests for one BC configuration
  ! ================================================================
  subroutine run_config(config_id, config_name, dims_global, &
                        BC_x, BC_y, BC_z)
    integer, intent(in) :: config_id
    character(len=*), intent(in) :: config_name
    integer, intent(in) :: dims_global(3)
    character(len=*), intent(in) :: BC_x(2), BC_y(2), BC_z(2)

    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    type(allocator_t), pointer :: host_allocator
    type(mesh_t), target :: mesh
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    type(vector_calculus_t) :: vector_calculus

#ifdef CUDA
    type(cuda_backend_t), target :: cuda_backend
    type(cuda_allocator_t), target :: cuda_allocator
#else
    type(omp_backend_t), target :: omp_backend
#endif
    type(allocator_t), target :: omp_allocator

    integer :: nproc_dir(3)
    real(dp) :: L_global(3)
    logical :: use_2decomp

    integer :: n, t, idx
    logical :: passed
    logical :: x_periodic, y_periodic, z_periodic
    integer :: test_types(NUM_TYPES), test_ns(NUM_NS)

    if (nrank == 0) then
      write (stderr, '(A)') ''
      write (stderr, '(A,A)') '  === Config ', config_name
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
#ifdef CUDA
    use_2decomp = .false.
#else
    use_2decomp = .true.
#endif

    mesh = mesh_t(dims_global, nproc_dir, L_global, &
                  BC_x, BC_y, BC_z, &
                  use_2decomp=use_2decomp)

#ifdef CUDA
    cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
    allocator => cuda_allocator

    omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
    host_allocator => omp_allocator

    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend
#else
    omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
    allocator => omp_allocator
    host_allocator => omp_allocator

    omp_backend = omp_backend_t(mesh, allocator)
    backend => omp_backend
#endif

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
    x_periodic = (trim(BC_x(1)) == 'periodic')
    y_periodic = (trim(BC_y(1)) == 'periodic')
    z_periodic = (trim(BC_z(1)) == 'periodic')

    ! Run all 8 cosine tests
    test_types = [TEST_COS_X, TEST_COS_Y, TEST_COS_XY, TEST_COS_XYZ]
    test_ns = [2, 3]
    idx = 0

    do n = 1, NUM_NS
      do t = 1, NUM_TYPES
        idx = idx + 1

        all_xfail(idx, config_id) = is_expected_fail( &
          test_ns(n), test_types(t), x_periodic, y_periodic, z_periodic)

        call run_single_test(backend, host_allocator, mesh, &
                             xdirps, ydirps, zdirps, vector_calculus, &
                             test_ns(n), test_types(t), passed, &
                             all_poisson_errs(idx, config_id), &
                             all_divgrad_errs(idx, config_id))
        all_results(idx, config_id) = passed
      end do
    end do

  end subroutine run_config

  ! ================================================================
  ! Print results for one config in the grand summary
  ! ================================================================
  subroutine print_config_results(config_id)
    integer, intent(in) :: config_id
    integer :: n, t, idx
    integer :: test_types(NUM_TYPES), test_ns(NUM_NS)
    character(len=4) :: result_str, expected_str
    character(len=2) :: verdict_str
    logical :: passed, xfail

    test_types = [TEST_COS_X, TEST_COS_Y, TEST_COS_XY, TEST_COS_XYZ]
    test_ns = [2, 3]
    idx = 0
    do n = 1, NUM_NS
      do t = 1, NUM_TYPES
        idx = idx + 1
        passed = all_results(idx, config_id)
        xfail = all_xfail(idx, config_id)

        result_str = merge('PASS', 'FAIL', passed)
        expected_str = merge('FAIL', 'PASS', xfail)

        if (passed .neqv. xfail) then
          verdict_str = 'OK'
        else
          verdict_str = '!!'
        end if

        write (stderr, '(2X,A5,2X,A10,I4,ES14.4,ES14.4,2X,A4,4X,A4,4X,A)') &
          config_labels(config_id), &
          test_type_name(test_types(t)), test_ns(n), &
          all_poisson_errs(idx, config_id), &
          all_divgrad_errs(idx, config_id), &
          result_str, expected_str, trim(verdict_str)
      end do
    end do
  end subroutine print_config_results

  ! ================================================================
  ! Helper: test type name
  ! ================================================================
  pure function test_type_name(test_type) result(name)
    integer, intent(in) :: test_type
    character(len=10) :: name

    select case (test_type)
    case (TEST_COS_X); name = 'COS_X     '
    case (TEST_COS_Y); name = 'COS_Y     '
    case (TEST_COS_XY); name = 'COS_XY    '
    case (TEST_COS_XYZ); name = 'COS_XYZ   '
    case default; name = 'UNKNOWN   '
    end select
  end function test_type_name

  pure function is_expected_fail(n_wave, test_type, &
                                 x_periodic, y_periodic, z_periodic) &
                                 result(xfail)
    !! Determine if a test is expected to fail.
    !!
    !! n=3 on even-sized periodic grids (64 cells) does not resolve
    !! cos(3*pi*x) cleanly due to aliasing. A test is XFAIL when n=3
    !! AND any direction involved in the test function uses periodic BCs.
    integer, intent(in) :: n_wave, test_type
    logical, intent(in) :: x_periodic, y_periodic, z_periodic
    logical :: xfail

    xfail = .false.
    if (n_wave /= 3) return

    select case (test_type)
    case (TEST_COS_X)
      xfail = x_periodic
    case (TEST_COS_Y)
      xfail = y_periodic
    case (TEST_COS_XY)
      xfail = x_periodic .or. y_periodic
    case (TEST_COS_XYZ)
      xfail = x_periodic .or. y_periodic .or. z_periodic
    end select
  end function is_expected_fail

  ! ================================================================
  ! Create cosine test field
  ! ================================================================
  subroutine create_cosine_field(mesh, host_field, n_wave, test_type)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n_wave
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), n_pi

    dims = mesh%get_dims(CELL)
    n_pi = real(n_wave, dp)*pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = mesh%get_coordinates(i, j, k, CELL)
          select case (test_type)
          case (TEST_COS_X)
            host_field%data(i, j, k) = cos(n_pi*coords(1))
          case (TEST_COS_Y)
            host_field%data(i, j, k) = cos(n_pi*coords(2))
          case (TEST_COS_XY)
            host_field%data(i, j, k) = cos(n_pi*coords(1)) &
                                       *cos(n_pi*coords(2))
          case (TEST_COS_XYZ)
            host_field%data(i, j, k) = cos(n_pi*coords(1)) &
                                       *cos(n_pi*coords(2)) &
                                       *cos(n_pi*coords(3))
          end select
        end do
      end do
    end do
  end subroutine create_cosine_field

  ! ================================================================
  ! Create analytical Poisson solution
  ! ================================================================
  subroutine create_analytical_solution(mesh, host_field, n_wave, test_type)
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n_wave
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), n_pi, n_pi_sq

    dims = mesh%get_dims(CELL)
    n_pi = real(n_wave, dp)*pi
    n_pi_sq = n_pi*n_pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = mesh%get_coordinates(i, j, k, CELL)
          select case (test_type)
          case (TEST_COS_X)
            host_field%data(i, j, k) = -cos(n_pi*coords(1))/n_pi_sq
          case (TEST_COS_Y)
            host_field%data(i, j, k) = -cos(n_pi*coords(2))/n_pi_sq
          case (TEST_COS_XY)
            host_field%data(i, j, k) = -cos(n_pi*coords(1)) &
                                       *cos(n_pi*coords(2)) &
                                       /(2.0_dp*n_pi_sq)
          case (TEST_COS_XYZ)
            host_field%data(i, j, k) = -cos(n_pi*coords(1)) &
                                       *cos(n_pi*coords(2)) &
                                       *cos(n_pi*coords(3)) &
                                       /(3.0_dp*n_pi_sq)
          end select
        end do
      end do
    end do
  end subroutine create_analytical_solution

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
                             n_wave, test_type, test_passed, &
                             poisson_err_out, divgrad_err_out)
    class(base_backend_t), pointer, intent(in) :: backend
    type(allocator_t), pointer, intent(in) :: host_allocator
    type(mesh_t), intent(in) :: mesh
    type(dirps_t), pointer, intent(in) :: xdirps, ydirps, zdirps
    type(vector_calculus_t), intent(in) :: vector_calculus
    integer, intent(in) :: n_wave
    integer, intent(in) :: test_type
    logical, intent(out) :: test_passed
    real(dp), intent(out) :: poisson_err_out, divgrad_err_out

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, host_analytical, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3)
    real(dp) :: poisson_error_norm, div_grad_error_norm
    logical :: poisson_passed, div_grad_passed

    dims = mesh%get_dims(CELL)

    if (mesh%par%is_root()) then
      write (stderr, '(4X,A,A,A,I1)') &
        'Running: ', trim(test_type_name(test_type)), '  n = ', n_wave
    end if

    ! Allocate fields
    f_device => backend%allocator%get_block(DIR_C, CELL)
    f_reference => backend%allocator%get_block(DIR_X)
    host_field => host_allocator%get_block(DIR_C)

    ! Create test function on host and transfer to device
    call create_cosine_field(mesh, host_field, n_wave, test_type)
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
    call create_analytical_solution(mesh, host_analytical, n_wave, test_type)

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

    div_grad_passed = (div_grad_error_norm <= ERROR_TOLERANCE)

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