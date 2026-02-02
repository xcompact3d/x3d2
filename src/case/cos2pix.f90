module m_case_cos2pix
  !! Poisson Solver Validation Test Case
  !!
  !! Validates the Poisson solver using the pressure_correction pattern:
  !!   1. Create test function: f = cos(n*pi*x/L) along Dirichlet axis
  !!   2. Solve Poisson equation: laplacian(p) = f
  !!   3. Compute gradient of solution: grad(p)
  !!   4. Compute divergence: div(grad(p))
  !!   5. Verify: div(grad(p)) ≈ f (within tolerance)
  !!
  !! Tests both cos(2*pi*x) and cos(3*pi*x)
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

  type, extends(base_case_t) :: case_cos2pix_t
    private
    integer :: dirichlet_axis(3) = 0
    real(dp) :: domain_length = 0.0_dp
  contains
    procedure :: boundary_conditions => boundary_conditions_cos2pix
    procedure :: initial_conditions => initial_conditions_cos2pix
    procedure :: forcings => forcings_cos2pix
    procedure :: pre_correction => pre_correction_cos2pix
    procedure :: postprocess => postprocess_cos2pix
    procedure :: run => run_cos2pix
    procedure, private :: detect_dirichlet_axis
    procedure, private :: create_cosine_field
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

  subroutine detect_dirichlet_axis(self)
    !! Detect which axis has Dirichlet boundary conditions
    class(case_cos2pix_t), intent(inout) :: self
    associate(bc => self%solver%mesh%grid%periodic_BC, L => self%solver%mesh%geo%L)
      if (.not. bc(1)) then
        self%dirichlet_axis(1)= 1
        self%domain_length = L(1)
      end if
      if (.not. bc(2)) then
        self%dirichlet_axis(2) = 1
        self%domain_length = L(2)
      end if
      if (.not. bc(3)) then
        self%dirichlet_axis(3)= 1
        self%domain_length = L(3)
      end if
    end associate
  end subroutine detect_dirichlet_axis

  subroutine create_cosine_field(self, host_field, n)
    !! Create test field: f = cos(n*pi*x/L) along Dirichlet axis
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field
    integer, intent(in) :: n

    integer :: i, j, k, dims(3), axis
    real(dp) :: coords(3), arg, result

    dims = self%solver%mesh%get_dims(CELL)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
          result = 1.0_dp
          do axis = 1, 3
            if (self%dirichlet_axis(axis) == 1) then
              arg = real(n, dp) * pi * coords(axis) / self%solver%mesh%geo%L(axis)
              result = result * cos(arg)
            end if
          end do
          host_field%data(i, j, k) = result
        end do
      end do
    end do
  end subroutine create_cosine_field

  function compute_error_norm(self, field) result(error_norm)
    !! Compute normalized error norm
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: field
    real(dp) :: error_norm

    integer :: dims(3)

    dims = self%solver%mesh%get_dims(CELL)
    error_norm = norm2(field%data(1:dims(1), 1:dims(2), 1:dims(3))) / product(dims)
  end function compute_error_norm

  subroutine print_field_comparison(self, host_field, label, n)
    !! Print field values vs analytical solution for debugging
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: host_field
    character(len=*), intent(in) :: label
    integer, intent(in) :: n

    integer :: idx, dims(3)
    integer :: Lx, Ly, Lz
    integer :: bc_pattern
    real(dp) :: coord, numerical, analytical, coord_x, coord_y, coord_z
    real(dp) :: n_pi  ! n*pi factor

    if (.not. self%solver%mesh%par%is_root()) return

    dims = self%solver%mesh%get_dims(CELL)
    n_pi = real(n, dp) * pi

    print *, ''
    print *, label
    print *, 'Testing cos(', n, '*pi*x/L)'
    print *, 'Coordinate | Numerical | Analytical'
    ! Convert to bit pattern: z + y*2 + x*4 gives values 0-7
    ! This way [1,0,0] (x only) = 4, [0,1,0] (y only) = 2, [0,0,1] (z only) = 1
    bc_pattern = self%dirichlet_axis(3) + self%dirichlet_axis(2)*2 + self%dirichlet_axis(1)*4
    print *, bc_pattern
    select case (bc_pattern)
      case (2) ! [0,1,0] - y only
        do idx = 1, dims(2)
          coord = self%solver%mesh%geo%midp_coords(idx, 2)
          numerical = host_field%data(4, idx, 8)
          analytical = cos(n_pi * coord) / (-n_pi * n_pi)
          print *, coord, numerical, analytical
        end do
      case (4) ![1,0,0] - x only
        do idx = 1, dims(1)
          coord = self%solver%mesh%geo%midp_coords(idx, 1)
          numerical = host_field%data(idx, 4, 8)
          analytical = cos(n_pi * coord) / (-n_pi * n_pi)
          print *, coord, numerical, analytical
        end do
      case (6) ! [1,1,0] - x and y
        Lx = 1
        Ly = 1
        
        print *, 'Printing along x-axis (at j=4, k=8):'
        do idx = 1, dims(1)
          coord_x = self%solver%mesh%geo%midp_coords(idx, 1)
          coord_y = self%solver%mesh%geo%midp_coords(4, 2)
          numerical = host_field%data(idx, 4, 8)
          analytical = cos(n_pi * coord_x / Lx) * cos(n_pi * coord_y / Ly) &
                       / (-n_pi * n_pi * (1.0_dp/(Lx*Lx) + 1.0_dp/(Ly*Ly)))
          print *, coord_x, numerical, analytical
        end do
        
        print *, ''
        print *, 'Printing along y-axis (at i=4, k=8):'
        do idx = 1, dims(2)
          coord_x = self%solver%mesh%geo%midp_coords(4, 1)
          coord_y = self%solver%mesh%geo%midp_coords(idx, 2)
          numerical = host_field%data(4, idx, 8)
          analytical = cos(n_pi * coord_x / Lx) * cos(n_pi * coord_y / Ly) &
                       / (-n_pi * n_pi * (1.0_dp/(Lx*Lx) + 1.0_dp/(Ly*Ly)))
          print *, coord_y, numerical, analytical
        end do
      case (3)
        error stop "Z-axis comparison not implemented"
    end select
  end subroutine print_field_comparison

  subroutine run_single_test(self, n, test_passed)
    !! Run a single test with cos(n*pi*x/L)
    class(case_cos2pix_t), intent(inout) :: self
    integer, intent(in) :: n
    logical, intent(out) :: test_passed

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3)
    real(dp) :: error_norm

    dims = self%solver%mesh%get_dims(CELL)

    ! Print test info
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '-----------------------------------------'
      print *, '  Testing cos(', n, '*pi*x/L)'
      print *, '-----------------------------------------'
      print *, ''
    end if

    ! Allocate fields
    f_device => self%solver%backend%allocator%get_block(DIR_C, CELL)
    f_reference => self%solver%backend%allocator%get_block(DIR_X)
    host_field => self%solver%host_allocator%get_block(DIR_C)

    ! Create test function f = cos(n*pi*x/L)
    call self%create_cosine_field(host_field, n)
    call self%print_field_comparison(host_field, 'Initial cosine field:', n)

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

    ! Debug: print Poisson solution
    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)
    call self%print_field_comparison(host_field, 'Poisson solution:', n)
    call self%solver%host_allocator%release_block(host_field)

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

    ! Compute error: result - reference
    call self%solver%backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

    ! Get error field to host and compute norm
    host_field => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_field%data, f_device)
    error_norm = self%compute_error_norm(host_field)

    ! Cleanup
    call self%solver%backend%allocator%release_block(f_device)
    call self%solver%backend%allocator%release_block(f_reference)
    call self%solver%host_allocator%release_block(host_field)

    ! Report results
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, 'Error:', error_norm
      print *, 'Tolerance:      ', ERROR_TOLERANCE
      print *, ''
    end if

    ! Check pass/fail
    test_passed = (error_norm <= ERROR_TOLERANCE)

    if (.not. test_passed) then
      if (self%solver%mesh%par%is_root()) then
        print *, 'TEST FAILED for cos(', n, '*pi*x/L)'
      end if
    else
      if (self%solver%mesh%par%is_root()) then
        print *, 'TEST PASSED for cos(', n, '*pi*x/L)'
      end if
    end if

  end subroutine run_single_test

  subroutine run_cos2pix(self)
    class(case_cos2pix_t), intent(inout) :: self

    integer :: dims(3)
    logical :: test_passed_2, test_passed_3

    ! Setup
    call self%detect_dirichlet_axis()
    dims = self%solver%mesh%get_dims(CELL)

    ! Print test info
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '========================================='
      print *, '    POISSON SOLVER VALIDATION TEST'
      print *, '========================================='
      print *, ''
      print *, 'Test: div(grad(poisson_solve(f))) ≈ f'
      print *, ''
      print *, 'Grid dimensions:', dims
      print *, 'Dirichlet axis: ', self%dirichlet_axis
      print *, 'Domain length:  ', self%domain_length
      print *, ''
    end if

    ! Run test for cos(2*pi*x)
    call self%run_single_test(2, test_passed_2)

    ! Run test for cos(3*pi*x)
    call self%run_single_test(3, test_passed_3)

    ! Summary
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '========================================='
      print *, '           TEST SUMMARY'
      print *, '========================================='
      print *, ''
      print *, 'cos(2*pi*x/L) test:', merge('PASSED', 'FAILED', test_passed_2)
      print *, 'cos(3*pi*x/L) test:', merge('PASSED', 'FAILED', test_passed_3)
      print *, ''
    end if

    ! Check overall pass/fail
    if (.not. (test_passed_2 .and. test_passed_3)) then
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