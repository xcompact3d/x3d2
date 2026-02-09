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
  !! Runs 6 test cases:
  !!   TEST_COS_X  -> cos(n*pi*x)             for n = 2, 3
  !!   TEST_COS_Y  -> cos(n*pi*y)             for n = 2, 3
  !!   TEST_COS_XY -> cos(n*pi*x)*cos(n*pi*y) for n = 2, 3
  !!
  !! Analytical Poisson solutions (laplacian(p) = f):
  !!   COS_X  : p = -cos(n*pi*x) / (n*pi)^2
  !!   COS_Y  : p = -cos(n*pi*y) / (n*pi)^2
  !!   COS_XY : p = -cos(n*pi*x)*cos(n*pi*y) / (2*(n*pi)^2)
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
  integer, parameter :: TEST_COS_X  = 1
  integer, parameter :: TEST_COS_Y  = 2
  integer, parameter :: TEST_COS_XY = 3

  type, extends(base_case_t) :: case_cos2pix_t
    private
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

  pure function test_type_name(test_type) result(name)
    integer, intent(in) :: test_type
    character(len=10) :: name
    select case (test_type)
    case (TEST_COS_X);  name = 'COS_X'
    case (TEST_COS_Y);  name = 'COS_Y'
    case (TEST_COS_XY); name = 'COS_XY'
    case default;        name = 'UNKNOWN'
    end select
  end function test_type_name

  subroutine create_cosine_field(self, host_field, n, test_type)
    !! Create test field based on test_type:
    !!   TEST_COS_X  -> cos(n*pi*x)
    !!   TEST_COS_Y  -> cos(n*pi*y)
    !!   TEST_COS_XY -> cos(n*pi*x) * cos(n*pi*y)
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), n_pi

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp) * pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
          select case (test_type)
          case (TEST_COS_X)
            host_field%data(i, j, k) = cos(n_pi * coords(1))
          case (TEST_COS_Y)
            host_field%data(i, j, k) = cos(n_pi * coords(2))
          case (TEST_COS_XY)
            host_field%data(i, j, k) = cos(n_pi * coords(1)) &
                                      * cos(n_pi * coords(2))
          end select
        end do
      end do
    end do
  end subroutine create_cosine_field

  subroutine create_analytical_solution(self, host_field, n, test_type)
    !! Create the analytical Poisson solution for laplacian(p) = f:
    !!   TEST_COS_X  : p = -cos(n*pi*x) / (n*pi)^2
    !!   TEST_COS_Y  : p = -cos(n*pi*y) / (n*pi)^2
    !!   TEST_COS_XY : p = -cos(n*pi*x)*cos(n*pi*y) / (2*(n*pi)^2)
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n
    integer, intent(in) :: test_type

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), n_pi, n_pi_sq

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp) * pi
    n_pi_sq = n_pi * n_pi

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
          select case (test_type)
          case (TEST_COS_X)
            host_field%data(i, j, k) = -cos(n_pi * coords(1)) / n_pi_sq
          case (TEST_COS_Y)
            host_field%data(i, j, k) = -cos(n_pi * coords(2)) / n_pi_sq
          case (TEST_COS_XY)
            host_field%data(i, j, k) = -cos(n_pi * coords(1)) &
                                      * cos(n_pi * coords(2)) &
                                      / (n_pi_sq + n_pi_sq)
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
    error_norm = norm2(field%data(1:dims(1), 1:dims(2), 1:dims(3))) / product(dims)
  end function compute_error_norm

  subroutine print_field_comparison(self, host_field, label, n, test_type, &
                                    is_poisson_solution)
    !! Print field values vs analytical solution for debugging
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: host_field
    character(len=*), intent(in) :: label
    integer, intent(in) :: n
    integer, intent(in) :: test_type
    logical, intent(in) :: is_poisson_solution

    integer :: ix, iy, iz, dims(3)
    real(dp) :: n_pi
    real(dp) :: coord_x, coord_y, numerical, analytical

    if (.not. self%solver%mesh%par%is_root()) return

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp) * pi

    print *, ''
    print *, label
    print *, 'Test type: ', test_type_name(test_type), '  n =', n
    print *, '  n_pi =', n_pi

    select case (test_type)

    case (TEST_COS_X)
      print *, 'ix | x-coord | Numerical | Analytical'
      iy = dims(2) / 2
      iz = dims(3) / 2
      do ix = 1, dims(1)
        coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
        numerical = host_field%data(ix, iy, iz)
        if (is_poisson_solution) then
          analytical = -cos(n_pi * coord_x) / (n_pi * n_pi)
        else
          analytical = cos(n_pi * coord_x)
        end if
        print *, ix, coord_x, numerical, analytical
      end do

    case (TEST_COS_Y)
      print *, 'iy | y-coord | Numerical | Analytical'
      ix = dims(1) / 2
      iz = dims(3) / 2
      do iy = 1, dims(2)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        numerical = host_field%data(ix, iy, iz)
        if (is_poisson_solution) then
          analytical = -cos(n_pi * coord_y) / (n_pi * n_pi)
        else
          analytical = cos(n_pi * coord_y)
        end if
        print *, iy, coord_y, numerical, analytical
      end do

    case (TEST_COS_XY)
      print *, 'ix | x-coord | y-coord | Numerical | Analytical'
      iz = dims(3) / 2
      iy = 0
      do ix = 1, dims(1)
        iy = iy + 1
        coord_x = self%solver%mesh%geo%midp_coords(iy, 1)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        numerical = host_field%data(iy, iy, iz)
        if (is_poisson_solution) then
          analytical = -cos(n_pi * coord_x) * cos(n_pi * coord_y) &
                      / (n_pi*n_pi + n_pi*n_pi)
        else
          analytical = cos(n_pi * coord_x) * cos(n_pi * coord_y)
        end if
        print *, iy, coord_x, coord_y, numerical, analytical
      end do

    end select

  end subroutine print_field_comparison

  subroutine run_single_test(self, n, test_type, test_passed)
    !! Run a single test with the given test_type and wavenumber n.
    !! Two checks are performed:
    !!   1. Poisson solution vs analytical (L2 norm of numerical - analytical)
    !!   2. div(grad(p)) vs original RHS f (L2 norm)
    !! Both must pass for the test to pass.
    class(case_cos2pix_t), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in) :: test_type
    logical, intent(out) :: test_passed

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, host_analytical, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3)
    real(dp) :: poisson_error_norm, div_grad_error_norm
    logical :: poisson_passed, div_grad_passed

    dims = self%solver%mesh%get_dims(CELL)

    ! Print test info
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '-----------------------------------------'
      print *, '  Testing ', test_type_name(test_type), ' with n =', n
      print *, '-----------------------------------------'
      print *, ''
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
    ! Get numerical solution to host
    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)

    ! Create analytical solution on host
    host_analytical => self%solver%host_allocator%get_block(DIR_C)
    call self%create_analytical_solution(host_analytical, n, test_type)

    ! Compute pointwise difference: host_field = numerical - analytical
    host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
      host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

    poisson_error_norm = self%compute_error_norm(host_field)

    call self%solver%host_allocator%release_block(host_analytical)
    call self%solver%host_allocator%release_block(host_field)

    poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE)

    if (self%solver%mesh%par%is_root()) then
      print *, '[Poisson solve] Error norm (L2):', poisson_error_norm
      print *, '[Poisson solve] Tolerance:      ', ERROR_TOLERANCE
      if (poisson_passed) then
        print *, '[Poisson solve] PASSED'
      else
        print *, '[Poisson solve] FAILED'
      end if
      print *, ''
    end if

    ! ---- Check 2: div(grad(p)) vs original RHS ----
    ! Compute gradient of pressure
    gradient_input => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%backend%reorder(gradient_input, f_device, RDR_C2Z)
    call self%solver%backend%allocator%release_block(f_device)

    dpdx => self%solver%backend%allocator%get_block(DIR_X)
    dpdy => self%solver%backend%allocator%get_block(DIR_X)
    dpdz => self%solver%backend%allocator%get_block(DIR_X)

    call self%solver%gradient_p2v(dpdx, dpdy, dpdz, gradient_input)
    call self%solver%backend%allocator%release_block(gradient_input)

    ! Compute divergence of gradient
    f_result => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%divergence_v2p(f_result, dpdx, dpdy, dpdz)

    call self%solver%backend%allocator%release_block(dpdx)
    call self%solver%backend%allocator%release_block(dpdy)
    call self%solver%backend%allocator%release_block(dpdz)

    ! Reorder result for comparison
    f_device => self%solver%backend%allocator%get_block(DIR_X)
    call self%solver%backend%reorder(f_device, f_result, RDR_Z2X)
    call self%solver%backend%allocator%release_block(f_result)

    ! Compute error: div(grad(p)) - f
    call self%solver%backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

    ! Get error field to host and compute norm
    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)
    div_grad_error_norm = self%compute_error_norm(host_field)

    ! Cleanup
    call self%solver%backend%allocator%release_block(f_device)
    call self%solver%backend%allocator%release_block(f_reference)
    call self%solver%host_allocator%release_block(host_field)

    div_grad_passed = (div_grad_error_norm <= ERROR_TOLERANCE)

    if (self%solver%mesh%par%is_root()) then
      print *, '[div(grad(p))] Error norm (L2):', div_grad_error_norm
      print *, '[div(grad(p))] Tolerance:      ', ERROR_TOLERANCE
      if (div_grad_passed) then
        print *, '[div(grad(p))] PASSED'
      else
        print *, '[div(grad(p))] FAILED'
      end if
      print *, ''
    end if

    ! Both checks must pass
    test_passed = poisson_passed .and. div_grad_passed

    if (self%solver%mesh%par%is_root()) then
      if (test_passed) then
        print *, 'TEST PASSED: ', test_type_name(test_type), ' n =', n
      else
        print *, 'TEST FAILED: ', test_type_name(test_type), ' n =', n
      end if
    end if

  end subroutine run_single_test

  subroutine run_cos2pix(self)
    class(case_cos2pix_t), intent(inout) :: self

    integer :: dims(3), n, t
    logical :: passed, all_passed
    logical :: results(6)
    integer :: test_types(3), test_ns(2)
    character(len=10) :: names(6)
    integer :: idx

    dims = self%solver%mesh%get_dims(CELL)

    ! Print test info
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '========================================='
      print *, '    POISSON SOLVER VALIDATION TEST'
      print *, '========================================='
      print *, ''
      print *, 'Test: poisson_solve(f) ≈ analytical'
      print *, '      div(grad(poisson_solve(f))) ≈ f'
      print *, ''
      print *, 'Grid dimensions:', dims
      print *, ''
    end if

    test_types = [TEST_COS_X, TEST_COS_Y, TEST_COS_XY]
    test_ns = [2, 3]

    idx = 0
    all_passed = .true.

    ! Run all 6 tests: for each n, run COS_X, COS_Y, COS_XY
    do n = 1, 2
      do t = 1, 3
        idx = idx + 1
        call self%run_single_test(test_ns(n), test_types(t), passed)
        results(idx) = passed
        names(idx) = test_type_name(test_types(t))
        all_passed = all_passed .and. passed
      end do
    end do

    ! Summary
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '========================================='
      print *, '           TEST SUMMARY'
      print *, '========================================='
      print *, ''

      idx = 0
      do n = 1, 2
        do t = 1, 3
          idx = idx + 1
          print *, names(idx), ' n=', test_ns(n), ': ', &
                   merge('PASSED', 'FAILED', results(idx))
        end do
      end do

      print *, ''
    end if

    ! Check overall pass/fail
    if (.not. all_passed) then
      error stop 'TEST FAILED: One or more tests did not pass'
    else
      if (self%solver%mesh%par%is_root()) then
        print *, '========================================='
        print *, '          ALL TESTS PASSED'
        print *, '========================================='
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