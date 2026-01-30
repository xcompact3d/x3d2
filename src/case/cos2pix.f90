module m_case_cos2pix
  !! Poisson Solver Validation Test Case
  !!
  !! Validates the Poisson solver using the pressure_correction pattern:
  !!   1. Create test function: f = cos(2*pi*x/L) along Dirichlet axis
  !!   2. Solve Poisson equation: laplacian(p) = f
  !!   3. Compute gradient of solution: grad(p)
  !!   4. Compute divergence: div(grad(p))
  !!   5. Verify: div(grad(p)) ≈ f (within tolerance)
  !!
  !! Analytical solution: p = -cos(2*pi*x/L) / (4*pi^2)
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

  real(dp), parameter :: ERROR_TOLERANCE = 1.0e-9_dp

  type, extends(base_case_t) :: case_cos2pix_t
    private
    integer :: dirichlet_axis = 0
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

    associate(bc => self%solver%mesh%geo%BC, L => self%solver%mesh%geo%L)
      if (bc%x(1) == 'dirichlet' .and. bc%x(2) == 'dirichlet') then
        self%dirichlet_axis = 1
        self%domain_length = L(1)
      else if (bc%y(1) == 'dirichlet' .and. bc%y(2) == 'dirichlet') then
        self%dirichlet_axis = 2
        self%domain_length = L(2)
      else if (bc%z(1) == 'dirichlet' .and. bc%z(2) == 'dirichlet') then
        self%dirichlet_axis = 3
        self%domain_length = L(3)
      else
        error stop "No axis with Dirichlet BCs on both ends found"
      end if
    end associate
  end subroutine detect_dirichlet_axis

  subroutine create_cosine_field(self, host_field)
    !! Create test field: f = cos(2*pi*x/L) along Dirichlet axis
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(inout) :: host_field

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), arg

    dims = self%solver%mesh%get_dims(CELL)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k, CELL)
          arg = 2.0_dp * pi * coords(self%dirichlet_axis) / self%domain_length
          host_field%data(i, j, k) = cos(arg)
        end do
      end do
    end do
  end subroutine create_cosine_field

  function compute_error_norm(self, field) result(error_norm)
    !! Compute normalized L2 error norm
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: field
    real(dp) :: error_norm

    integer :: dims(3)

    dims = self%solver%mesh%get_dims(CELL)
    error_norm = norm2(field%data(1:dims(1), 1:dims(2), 1:dims(3))) / product(dims)
  end function compute_error_norm

  subroutine print_field_comparison(self, host_field, label)
    !! Print field values vs analytical solution for debugging
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: host_field
    character(len=*), intent(in) :: label

    integer :: idx, dims(3)
    real(dp) :: coord, numerical, analytical

    if (.not. self%solver%mesh%par%is_root()) return

    dims = self%solver%mesh%get_dims(CELL)
    print *, ''
    print *, label
    print *, 'Coordinate | Numerical | Analytical'

    select case (self%dirichlet_axis)
    case (1)
      do idx = 1, dims(1)
        coord = self%solver%mesh%geo%midp_coords(idx, 1)
        numerical = host_field%data(idx, 4, 8)
        analytical = cos(2.0_dp * pi * coord) / (-4.0_dp * pi * pi)
        print *, coord, numerical, analytical
      end do
    case (2)
      do idx = 1, dims(2)
        coord = self%solver%mesh%geo%midp_coords(idx, 2)
        numerical = host_field%data(4, idx, 8)
        analytical = cos(2.0_dp * pi * coord) / (-4.0_dp * pi * pi)
        print *, coord, numerical, analytical
      end do
    case (3)
      error stop "Z-axis comparison not implemented"
    end select
  end subroutine print_field_comparison

  subroutine run_cos2pix(self)
    class(case_cos2pix_t), intent(inout) :: self

    class(field_t), pointer :: f_device, f_reference, f_result
    class(field_t), pointer :: host_field, temp
    class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input
    integer :: dims(3)
    real(dp) :: error_norm
    character(len=1) :: axis_label

    ! Setup
    call self%detect_dirichlet_axis()
    dims = self%solver%mesh%get_dims(CELL)
    axis_label = char(ichar('w') + self%dirichlet_axis)  ! 'x', 'y', or 'z'

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
      print *, 'Dirichlet axis: ', axis_label
      print *, 'Domain length:  ', self%domain_length
      print *, ''
    end if

    ! Allocate fields
    f_device => self%solver%backend%allocator%get_block(DIR_C, CELL)
    f_reference => self%solver%backend%allocator%get_block(DIR_X)
    host_field => self%solver%host_allocator%get_block(DIR_C)

    ! Create test function f = cos(2*pi*x/L)
    call self%create_cosine_field(host_field)
    call self%print_field_comparison(host_field, 'Initial cosine field:')

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
    call self%print_field_comparison(host_field, 'Poisson solution:')
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
      print *, 'Error norm (L2):', error_norm
      print *, 'Tolerance:      ', ERROR_TOLERANCE
      print *, ''
    end if

    ! Check pass/fail
    if (error_norm > ERROR_TOLERANCE) then
      error stop 'TEST FAILED: div(grad(poisson_solve(f))) != f'
    else
      if (self%solver%mesh%par%is_root()) then
        print *, '========================================='
        print *, '              TEST PASSED'
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