module m_case_cos2pix
  !! Poisson Solver Validation Test Case
  !!
  !! Validates the Poisson solver using the pressure_correction pattern:
  !!   1. Create test function: f = cos(n*pi*x/L) along specified axes
  !!   2. Solve Poisson equation: laplacian(p) = f
  !!   3. Check numerical solution matches analytical solution (L2 norm)
  !!   4. Compute gradient of solution: grad(p)
  !!   5. Compute divergence: div(grad(p))
  !!   6. Verify: div(grad(p)) ≈ f (within tolerance)
  !!
  !! Runs 8 test cases:
  !!   TEST_COS_X   -> cos(n*pi*x)                          for n = 2, 3
  !!   TEST_COS_Y   -> cos(n*pi*y)                          for n = 2, 3
  !!   TEST_COS_XY  -> cos(n*pi*x)*cos(n*pi*y)              for n = 2, 3
  !!   TEST_COS_XYZ -> cos(n*pi*x)*cos(n*pi*y)*cos(n*pi*z)  for n = 2, 3
  !!
  !! Analytical Poisson solutions (laplacian(p) = f):
  !!   COS_X   : p = -cos(n*pi*x) / (n*pi)^2
  !!   COS_Y   : p = -cos(n*pi*y) / (n*pi)^2
  !!   COS_XY  : p = -cos(n*pi*x)*cos(n*pi*y) / (2*(n*pi)^2)
  !!   COS_XYZ : p = -cos(n*pi*x)*cos(n*pi*y)*cos(n*pi*z) / (3*(n*pi)^2)
  !!
  !! Verbosity control via environment variable:
  !!   VERBOSE=1  -> print pointwise field comparisons (cosine & analytical)
  !!   (default)  -> only print error norms and pass/fail status
  !!
  !! Usage:
  !!   mpirun -n 1 ./build-gpu/bin/xcompact input.x3d              (quiet)
  !!   VERBOSE=1 mpirun -n 1 ./build-gpu/bin/xcompact input.x3d    (verbose)
  !!
  !! NOTE: For Dirichlet BCs, dims_global along that axis must be ODD (e.g., 65)

  use m_allocator
  use m_base_backend
  use m_base_case, only: base_case_t
  use m_common, only: pi, dp, DIR_C, DIR_X, DIR_Z, CELL, RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_field, only: field_t
  use m_mesh

  implicit none
  private

  public :: case_cos2pix_t

  real(dp), parameter :: ERROR_TOLERANCE = 1.0e-11_dp

  ! Test type identifiers
  integer, parameter :: TEST_COS_X = 1
  integer, parameter :: TEST_COS_Y = 2
  integer, parameter :: TEST_COS_XY = 3
  integer, parameter :: TEST_COS_XYZ = 4

  ! Number of test types and wavenumbers
  integer, parameter :: NUM_TYPES = 4
  integer, parameter :: NUM_NS = 2
  integer, parameter :: NUM_TESTS = NUM_TYPES*NUM_NS

  type, extends(base_case_t) :: case_cos2pix_t
    private
    logical :: verbose = .false.
  contains
    procedure :: boundary_conditions => boundary_conditions_cos2pix
    procedure :: initial_conditions => initial_conditions_cos2pix
    procedure :: forcings => forcings_cos2pix
    procedure :: pre_correction => pre_correction_cos2pix
    procedure :: postprocess => postprocess_cos2pix
    procedure :: run => run_cos2pix
    procedure, private :: create_cosine_field
    procedure, private :: create_analytical_solution
    procedure, private :: compute_error_norm
    procedure, private :: print_field_comparison
    procedure, private :: run_single_test
    procedure, private :: detect_verbose
  end type case_cos2pix_t

  interface case_cos2pix_t
    module procedure case_cos2pix_init
  end interface case_cos2pix_t

contains

  function case_cos2pix_init(backend, mesh, host_allocator) result(flow_case)
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_cos2pix_t) :: flow_case

    call flow_case%case_init(backend, mesh, host_allocator)
  end function case_cos2pix_init

  subroutine detect_verbose(self)
    !! Read the VERBOSE environment variable to set debug output level.
    !!   VERBOSE=1  -> verbose output (pointwise field comparisons)
    !!   (unset/0)  -> quiet output   (error norms and status only)
    class(case_cos2pix_t), intent(inout) :: self

    character(len=16) :: env_val
    integer :: stat

    call get_environment_variable('VERBOSE', env_val, status=stat)
    self%verbose = (stat == 0 .and. trim(env_val) == '1')

    if (self%solver%mesh%par%is_root()) then
      if (self%verbose) then
        write (*, '(4X,A)') 'Verbose   : ON  (VERBOSE=1 detected)'
      else
        write (*, '(4X,A)') 'Verbose   : OFF (set VERBOSE=1 for field dumps)'
      end if
    end if
  end subroutine detect_verbose

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

  pure function test_rhs_formula(test_type, n) result(formula)
    integer, intent(in) :: test_type, n
    character(len=60) :: formula
    character(len=1) :: cn

    write (cn, '(I1)') n
    select case (test_type)
    case (TEST_COS_X)
      formula = 'cos('//cn//'*pi*x)'
    case (TEST_COS_Y)
      formula = 'cos('//cn//'*pi*y)'
    case (TEST_COS_XY)
      formula = 'cos('//cn//'*pi*x)*cos('//cn//'*pi*y)'
    case (TEST_COS_XYZ)
      formula = 'cos('//cn//'*pi*x)*cos('//cn//'*pi*y)*cos('//cn//'*pi*z)'
    case default
      formula = '???'
    end select
  end function test_rhs_formula

  subroutine create_cosine_field(self, host_field, n, test_type)
    !! Create test field based on test_type:
    !!   TEST_COS_X   -> cos(n*pi*x)
    !!   TEST_COS_Y   -> cos(n*pi*y)
    !!   TEST_COS_XY  -> cos(n*pi*x) * cos(n*pi*y)
    !!   TEST_COS_XYZ -> cos(n*pi*x) * cos(n*pi*y) * cos(n*pi*z)
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3), x
    real(dp) :: coords(3), n_pi

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp)*pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
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

  subroutine create_analytical_solution(self, host_field, n, test_type)
    !! Create the analytical Poisson solution for laplacian(p) = f:
    !!   TEST_COS_X   : p = -cos(n*pi*x) / (n*pi)^2
    !!   TEST_COS_Y   : p = -cos(n*pi*y) / (n*pi)^2
    !!   TEST_COS_XY  : p = -cos(n*pi*x)*cos(n*pi*y) / (2*(n*pi)^2)
    !!   TEST_COS_XYZ : p = -cos(n*pi*x)*cos(n*pi*y)*cos(n*pi*z) / (3*(n*pi)^2)
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), n_pi, n_pi_sq

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp)*pi
    n_pi_sq = n_pi*n_pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
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

  function compute_error_norm(self, field) result(error_norm)
    !! Compute normalized L2 error norm
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: field
    real(dp) :: error_norm

    integer :: dims(3)

    dims = self%solver%mesh%get_dims(CELL)
  error_norm = norm2(field%data(1:dims(1), 1:dims(2), 1:dims(3)))/product(dims)
  end function compute_error_norm

  subroutine print_field_comparison(self, host_field, host_analytical, &
                                    label, n, test_type)
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: host_field
    class(field_t), intent(in) :: host_analytical
    character(len=*), intent(in) :: label
    integer, intent(in) :: n
    integer, intent(in) :: test_type

    integer :: ix, iy, iz, ii, dims(3)
    real(dp) :: n_pi
    real(dp) :: coord_x, coord_y, coord_z

    if (.not. self%verbose) return
    if (.not. self%solver%mesh%par%is_root()) return

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp)*pi

    write (*, '(A)') ''
    write (*, '(4X,A)') label
    write (*, '(4X,A,A,A,I1)') &
      'Test type : ', trim(test_type_name(test_type)), '   n = ', n
    write (*, '(4X,A,ES12.5)') 'n*pi      = ', n_pi

    select case (test_type)

    case (TEST_COS_X)
      write (*, '(A)') ''
      write (*, '(6X,A)') &
        '  ix    x-coord         Numerical           Analytical'
      write (*, '(A)') ''
      iy = dims(2)/2
      iz = dims(3)/2
      do ix = 1, dims(1)
        coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
        write (*, '(6X,I4,F12.6,2ES20.12)') ix, coord_x, &
          host_field%data(ix, iy, iz), host_analytical%data(ix, iy, iz)
      end do

    case (TEST_COS_Y)
      write (*, '(A)') ''
      write (*, '(6X,A)') &
        '  iy    y-coord         Numerical           Analytical'
      write (*, '(A)') ''
      ix = dims(1)/2
      iz = dims(3)/2
      do iy = 1, dims(2)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        write (*, '(6X,I4,F12.6,2ES20.12)') iy, coord_y, &
          host_field%data(ix, iy, iz), host_analytical%data(ix, iy, iz)
      end do

    case (TEST_COS_XY)
      write (*, '(A)') ''
      write (*, '(6X,A)') &
        '  ix    x-coord      y-coord         Numerical           Analytical'
      write (*, '(A)') ''
      iz = dims(3)/2
      do ii = 1, min(dims(1), dims(2))
        coord_x = self%solver%mesh%geo%midp_coords(ii, 1)
        coord_y = self%solver%mesh%geo%midp_coords(ii, 2)
        write (*, '(6X,I4,2F12.6,2ES20.12)') ii, coord_x, coord_y, &
          host_field%data(ii, ii, iz), host_analytical%data(ii, ii, iz)
      end do

    case (TEST_COS_XYZ)
      write (*, '(A)') ''
      write (*, '(6X,A)') &
        '  ii    x-coord      y-coord      z-coord' &
        //'         Numerical           Analytical'
      write (*, '(A)') ''
      do ii = 1, min(dims(1), dims(2), dims(3))
        coord_x = self%solver%mesh%geo%midp_coords(ii, 1)
        coord_y = self%solver%mesh%geo%midp_coords(ii, 2)
        coord_z = self%solver%mesh%geo%midp_coords(ii, 3)
        write (*, '(6X,I4,3F12.6,2ES20.12)') ii, coord_x, coord_y, coord_z, &
          host_field%data(ii, ii, ii), host_analytical%data(ii, ii, ii)
      end do

    end select

  end subroutine print_field_comparison

  subroutine run_single_test(self, n, test_type, test_passed, &
                             poisson_err_out, divgrad_err_out)
    !! Run a single test with the given test_type and wavenumber n.
    !! Two checks are performed:
    !!   1. Poisson solution vs analytical (L2 norm of numerical - analytical)
    !!   2. div(grad(p)) vs original RHS f (L2 norm)
    !! Both must pass for the test to pass.
    !! Error norms are returned for the summary table.
    class(case_cos2pix_t), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in) :: test_type
    logical, intent(out) :: test_passed
    real(dp), intent(out) :: poisson_err_out, divgrad_err_out

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, host_analytical, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3), x
    real(dp) :: poisson_error_norm, div_grad_error_norm
    logical :: poisson_passed, div_grad_passed

    dims = self%solver%mesh%get_dims(CELL)

    ! Print test header
    if (self%solver%mesh%par%is_root()) then
      write (*, '(A)') ''
      write (*, '(4X,A,A,A,I1)') &
        'Test  ', trim(test_type_name(test_type)), '   n = ', n
      write (*, '(4X,A,A)') &
        'RHS   f = ', trim(test_rhs_formula(test_type, n))
      write (*, '(4X,A)') &
        '========================================='
      write (*, '(A)') ''
    end if

    ! Allocate fields
    f_device => self%solver%backend%allocator%get_block(DIR_C, CELL)
    f_reference => self%solver%backend%allocator%get_block(DIR_X)
    host_field => self%solver%host_allocator%get_block(DIR_C)

    ! Create test function
    call self%create_cosine_field(host_field, n, test_type)

    ! Transfer to device
    call self%solver%backend%set_field_data(f_device, host_field%data, DIR_C)
    call f_device%set_data_loc(CELL)
    call self%solver%host_allocator%release_block(host_field)

    ! Store reference copy for later comparison
    call self%solver%backend%reorder(f_reference, f_device, RDR_C2X)

    ! Solve Poisson equation
    temp => self%solver%backend%allocator%get_block(DIR_C)
    call self%solver%backend%poisson_fft%solve_poisson(f_device, temp)
    call self%solver%backend%allocator%release_block(temp)

! ---- Check 1: Poisson solution vs analytical ----
    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)

    ! Remove arbitrary constant from numerical solution
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_field%data(1, 1, 1)

    host_analytical => self%solver%host_allocator%get_block(DIR_C)
    call self%create_analytical_solution(host_analytical, n, test_type)

    ! Remove same constant from analytical solution
    host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1, 1, 1)

    ! Verbose: print the Poisson solution vs analytical (both shifted)
    call self%print_field_comparison(host_field, host_analytical, &
             'Poisson solution (numerical vs analytical, constant removed):', &
                                     n, test_type)

    ! Compute pointwise difference
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

    poisson_error_norm = self%compute_error_norm(host_field)

    call self%solver%host_allocator%release_block(host_analytical)
    call self%solver%host_allocator%release_block(host_field)

    poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE)

    if (self%solver%mesh%par%is_root()) then
      write (*, '(6X,A)') 'Check 1  solve(f) = φ  vs  analytical'
      write (*, '(6X,A,ES14.6)') '  L2 error  : ', poisson_error_norm
      write (*, '(6X,A,ES14.6)') '  Tolerance : ', ERROR_TOLERANCE
      write (*, '(6X,A,A)') '  Status    : ', &
        merge('PASSED', 'FAILED', poisson_passed)
      write (*, '(A)') ''
    end if

    ! ---- Check 2: div(grad(p)) vs original RHS ----
    gradient_input => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%backend%reorder(gradient_input, f_device, RDR_C2Z)
    call self%solver%backend%allocator%release_block(f_device)

    dpdx => self%solver%backend%allocator%get_block(DIR_X)
    dpdy => self%solver%backend%allocator%get_block(DIR_X)
    dpdz => self%solver%backend%allocator%get_block(DIR_X)

    call self%solver%gradient_p2v(dpdx, dpdy, dpdz, gradient_input)
    call self%solver%backend%allocator%release_block(gradient_input)

    f_result => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%divergence_v2p(f_result, dpdx, dpdy, dpdz)

    call self%solver%backend%allocator%release_block(dpdx)
    call self%solver%backend%allocator%release_block(dpdy)
    call self%solver%backend%allocator%release_block(dpdz)

    f_device => self%solver%backend%allocator%get_block(DIR_X)
    call self%solver%backend%reorder(f_device, f_result, RDR_Z2X)
    call self%solver%backend%allocator%release_block(f_result)

    ! Compute error: div(grad(p)) - f
    call self%solver%backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)
    div_grad_error_norm = self%compute_error_norm(host_field)

    ! Cleanup
    call self%solver%backend%allocator%release_block(f_device)
    call self%solver%backend%allocator%release_block(f_reference)
    call self%solver%host_allocator%release_block(host_field)

    div_grad_passed = (div_grad_error_norm <= ERROR_TOLERANCE)

    if (self%solver%mesh%par%is_root()) then
      write (*, '(6X,A)') 'Check 2  ∇·(∇φ) = f  (round-trip)'
      write (*, '(6X,A,ES14.6)') '  L2 error  : ', div_grad_error_norm
      write (*, '(6X,A,ES14.6)') '  Tolerance : ', ERROR_TOLERANCE
      write (*, '(6X,A,A)') '  Status    : ', &
        merge('PASSED', 'FAILED', div_grad_passed)
      write (*, '(A)') ''
    end if

    ! Both checks must pass
    test_passed = poisson_passed .and. div_grad_passed
    poisson_err_out = poisson_error_norm
    divgrad_err_out = div_grad_error_norm

  end subroutine run_single_test

  subroutine run_cos2pix(self)
    class(case_cos2pix_t), intent(inout) :: self

    integer :: dims(3), n, t, idx
    logical :: passed, all_passed
    logical :: results(NUM_TESTS)
    integer :: test_types(NUM_TYPES), test_ns(NUM_NS)
    character(len=10) :: names(NUM_TESTS)
    real(dp) :: poisson_errs(NUM_TESTS), divgrad_errs(NUM_TESTS)

    dims = self%solver%mesh%get_dims(CELL)

    ! Detect verbosity from environment
    call self%detect_verbose()

    ! Print banner and schematic
    if (self%solver%mesh%par%is_root()) then
      write (*, '(A)') ''
      write (*, '(A)') &
        '  ========================================================='
      write (*, '(A)') &
        '            POISSON SOLVER VALIDATION TEST                  '
      write (*, '(A)') &
        '  ========================================================='
      write (*, '(A)') ''
      write (*, '(4X,A)') 'Pipeline under test:'
      write (*, '(A)') ''
      write (*, '(4X,A)') &
        'f  -->  solve(f) = φ  -->  ∇φ  -->  ∇·(∇φ) = ∇²φ = f'
      write (*, '(4X,A)') &
        '          |                  |              |'
      write (*, '(4X,A)') &
        '    "double integral"   1st derivative  1st derivative'
      write (*, '(4X,A)') &
        '    (inverse of ∇²)                  (together = ∇²'
      write (*, '(4X,A)') &
        '                                        = Laplacian)'
      write (*, '(A)') ''
      write (*, '(4X,A)') 'Verification checks (both must pass):'
      write (*, '(6X,A)') 'Check 1 : φ matches known analytical solution'
      write (*, '(6X,A)') 'Check 2 : ∇·(∇φ) recovers original f'
      write (*, '(A)') ''
      write (*, '(4X,A,I4,A,I4,A,I4)') &
        'Grid      : ', dims(1), ' x', dims(2), ' x', dims(3)
      write (*, '(4X,A,ES10.3)') &
        'Tolerance : ', ERROR_TOLERANCE
      write (*, '(A)') ''

      ! Test matrix as borderless table
      write (*, '(4X,A)') 'Test matrix (8 cases):'
      write (*, '(A)') ''
      write (*, '(6X,A10,A6,A3,A)') &
        'Type      ', '  n   ', '   ', 'RHS f(x,y,z)'
      write (*, '(A)') ''
      write (*, '(6X,A10,I6,A3,A)') 'COS_X     ', 2, '   ', 'cos(2πx)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_Y     ', 2, '   ', 'cos(2πy)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_XY    ', 2, '   ', &
        'cos(2πx)·cos(2πy)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_XYZ   ', 2, '   ', &
        'cos(2πx)·cos(2πy)·cos(2πz)'
      write (*, '(A)') ''
      write (*, '(6X,A10,I6,A3,A)') 'COS_X     ', 3, '   ', 'cos(3πx)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_Y     ', 3, '   ', 'cos(3πy)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_XY    ', 3, '   ', &
        'cos(3πx)·cos(3πy)'
      write (*, '(6X,A10,I6,A3,A)') 'COS_XYZ   ', 3, '   ', &
        'cos(3πx)·cos(3πy)·cos(3πz)'
      write (*, '(A)') ''
    end if

    test_types = [TEST_COS_X, TEST_COS_Y, TEST_COS_XY, TEST_COS_XYZ]
    test_ns = [2, 3]

    idx = 0
    all_passed = .true.

    ! Run all 8 tests: for each n, run COS_X, COS_Y, COS_XY, COS_XYZ
    do n = 1, NUM_NS
      do t = 1, NUM_TYPES
        idx = idx + 1
        call self%run_single_test(test_ns(n), test_types(t), passed, &
                                  poisson_errs(idx), divgrad_errs(idx))
        results(idx) = passed
        names(idx) = test_type_name(test_types(t))
        all_passed = all_passed .and. passed
      end do
    end do

    ! Summary results table (borderless)
    if (self%solver%mesh%par%is_root()) then
      write (*, '(A)') ''
      write (*, '(A)') &
        '  ========================================================='
      write (*, '(A)') &
        '                       TEST SUMMARY                        '
      write (*, '(A)') &
        '  ========================================================='
      write (*, '(A)') ''

      ! Column headers
      write (*, '(6X,A10,A4,A16,A16,A10)') &
        'Type      ', ' n  ', '  Poisson L2    ', '  ∇·(∇φ) L2     ', &
        '  Status  '
      write (*, '(A)') ''

      ! Data rows
      idx = 0
      do n = 1, NUM_NS
        do t = 1, NUM_TYPES
          idx = idx + 1
          write (*, '(6X,A10,I4,ES16.6,ES16.6,A4,A6)') &
            names(idx), test_ns(n), &
            poisson_errs(idx), divgrad_errs(idx), &
            '    ', merge('PASSED', 'FAILED', results(idx))
        end do
        ! Blank line between n=2 and n=3 groups
        if (n < NUM_NS) write (*, '(A)') ''
      end do

      write (*, '(A)') ''
      write (*, '(6X,A,ES10.3)') 'Tolerance : ', ERROR_TOLERANCE
      write (*, '(A)') ''
    end if

    ! Final verdict
    if (.not. all_passed) then
      if (self%solver%mesh%par%is_root()) then
        write (*, '(4X,A)') '!! ONE OR MORE TESTS FAILED !!'
        write (*, '(A)') ''
      end if
      error stop 'TEST FAILED: One or more tests did not pass'
    else
      if (self%solver%mesh%par%is_root()) then
        write (*, '(A)') &
          '  ========================================================='
        write (*, '(A)') &
          '                    ALL TESTS PASSED                        '
        write (*, '(A)') &
          '  ========================================================='
        write (*, '(A)') ''
      end if
    end if

    call self%case_finalise()
  end subroutine run_cos2pix

  ! Required interface implementations (empty for this test case)

  subroutine initial_conditions_cos2pix(self)
    class(case_cos2pix_t) :: self
  end subroutine initial_conditions_cos2pix

  subroutine boundary_conditions_cos2pix(self)
    class(case_cos2pix_t) :: self
  end subroutine boundary_conditions_cos2pix

  subroutine forcings_cos2pix(self, du, dv, dw, iter)
    class(case_cos2pix_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter
  end subroutine forcings_cos2pix

  subroutine pre_correction_cos2pix(self, u, v, w)
    class(case_cos2pix_t) :: self
    class(field_t), intent(inout) :: u, v, w
  end subroutine pre_correction_cos2pix

  subroutine postprocess_cos2pix(self, iter, t)
    class(case_cos2pix_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t
  end subroutine postprocess_cos2pix

end module m_case_cos2pix
