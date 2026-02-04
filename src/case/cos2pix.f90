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
    !! Compute normalized L2 error norm
    class(case_cos2pix_t), intent(in) :: self
    class(field_t), intent(in) :: field
    real(dp) :: error_norm

    integer :: dims(3)

    dims = self%solver%mesh%get_dims(CELL)
    error_norm = norm2(field%data(1:dims(1), 1:dims(2), 1:dims(3))) / product(dims)
  end function compute_error_norm

subroutine print_field_comparison(self, host_field, label, n, is_poisson_solution)
  !! Print field values vs analytical solution for debugging
  class(case_cos2pix_t), intent(in) :: self
  class(field_t), intent(in) :: host_field
  character(len=*), intent(in) :: label
  integer, intent(in) :: n
  logical, intent(in) :: is_poisson_solution

  integer :: ix, iy, iz, dims(3)
  real(dp) :: Lx, Ly, Lz
  integer :: bc_pattern
  real(dp) :: coord_x, coord_y, coord_z
  real(dp) :: numerical, analytical
  real(dp) :: kx, ky, kz
  real(dp) :: val_z1, val_zk
  real(dp) :: phi_100_analytical, phi_010_analytical, phi_110_analytical, phi_product
  real(dp) :: error_sum, error_max, error_l2
  integer :: count

  if (.not. self%solver%mesh%par%is_root()) return

  dims = self%solver%mesh%get_dims(CELL)

  ! Get domain lengths
  Lx = self%solver%mesh%geo%L(1)
  Ly = self%solver%mesh%geo%L(2)
  Lz = self%solver%mesh%geo%L(3)

  ! Wavenumbers
  kx = real(n, dp) * pi / Lx
  ky = real(n, dp) * pi / Ly
  kz = real(n, dp) * pi / Lz

  print *, ''
  print *, '##################################################'
  print *, label
  print *, '##################################################'
  print *, 'Domain: Lx=', Lx, ' Ly=', Ly, ' Lz=', Lz
  print *, 'Grid: nx=', dims(1), ' ny=', dims(2), ' nz=', dims(3)
  print *, 'Mode number n =', n
  print *, 'kx =', kx, ' ky =', ky, ' kz =', kz

  ! Convert to bit pattern: z + y*2 + x*4
  bc_pattern = self%dirichlet_axis(3) + self%dirichlet_axis(2)*2 + self%dirichlet_axis(1)*4
  print *, 'BC pattern: ', bc_pattern

  select case (bc_pattern)

  case (2)  ! [0,1,0] - y only (010 case)
    print *, ''
    print *, '=========================================='
    print *, '010 Case: Dirichlet Y, Periodic X and Z'
    print *, 'RHS: f = cos(ky*y)'
    print *, 'Analytical: phi = -cos(ky*y) / ky^2'
    print *, '=========================================='
    print *, ''
    
    error_sum = 0._dp
    error_max = 0._dp
    count = 0
    
    ! Fixed x and z indices (middle of domain)
    ix = dims(1) / 2
    iz = dims(3) / 2
    
    print *, 'Sampling along Y at ix=', ix, ', iz=', iz
    print *, ''
    print '(A12, A16, A16, A16)', 'y-coord', 'Numerical', 'Analytical', 'Error'
    print '(A12, A16, A16, A16)', '-------', '---------', '----------', '-----'
    
    do iy = 1, dims(2)
      coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
      
      ! Array access: data(x_index, y_index, z_index)
      numerical = host_field%data(ix, iy, iz)
      
      if (is_poisson_solution) then
        analytical = -cos(ky * coord_y) / (ky * ky)
      else
        analytical = cos(ky * coord_y)
      end if
      
      error_sum = error_sum + (numerical - analytical)**2
      error_max = max(error_max, abs(numerical - analytical))
      count = count + 1
      
      print '(F12.6, 3ES16.6)', coord_y, numerical, analytical, abs(numerical - analytical)
    end do
    
    error_l2 = sqrt(error_sum / count)
    
    print *, ''
    print *, '010 Case Summary:'
    print *, 'L2 error:  ', error_l2
    print *, 'Max error: ', error_max
    print *, '=========================================='

  case (4)  ! [1,0,0] - x only (100 case)
    print *, ''
    print *, '=========================================='
    print *, '100 Case: Dirichlet X, Periodic Y and Z'
    print *, 'RHS: f = cos(kx*x)'
    print *, 'Analytical: phi = -cos(kx*x) / kx^2'
    print *, '=========================================='
    print *, ''
    
    error_sum = 0._dp
    error_max = 0._dp
    count = 0
    
    ! Fixed y and z indices (middle of domain)
    iy = dims(2) / 2
    iz = dims(3) / 2
    
    print *, 'Sampling along X at iy=', iy, ', iz=', iz
    print *, ''
    print '(A12, A16, A16, A16)', 'x-coord', 'Numerical', 'Analytical', 'Error'
    print '(A12, A16, A16, A16)', '-------', '---------', '----------', '-----'
    
    do ix = 1, dims(1)
      coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
      
      ! Array access: data(x_index, y_index, z_index)
      numerical = host_field%data(ix, iy, iz)
      
      if (is_poisson_solution) then
        analytical = -cos(kx * coord_x) / (kx * kx)
      else
        analytical = cos(kx * coord_x)
      end if
      
      error_sum = error_sum + (numerical - analytical)**2
      error_max = max(error_max, abs(numerical - analytical))
      count = count + 1
      
      print '(F12.6, 3ES16.6)', coord_x, numerical, analytical, abs(numerical - analytical)
    end do
    
    error_l2 = sqrt(error_sum / count)
    
    print *, ''
    print *, '100 Case Summary:'
    print *, 'L2 error:  ', error_l2
    print *, 'Max error: ', error_max
    print *, '=========================================='

  case (6)  ! [1,1,0] - x and y (110 case)
    print *, ''
    print *, '=========================================='
    print *, '110 Case: Dirichlet X and Y, Periodic Z'
    print *, 'RHS: f = cos(kx*x) * cos(ky*y)'
    print *, 'Analytical: phi = -cos(kx*x)*cos(ky*y) / (kx^2 + ky^2)'
    print *, '=========================================='
    print *, ''
    print *, 'kx^2 = ', kx*kx
    print *, 'ky^2 = ', ky*ky
    print *, 'kx^2 + ky^2 = ', kx*kx + ky*ky
    
    ! =========================================================================
    ! CHECK 1: Z-independence (solution should be constant in z)
    ! =========================================================================
    print *, ''
    print *, '-------------------------------------------'
    print *, 'CHECK 1: Z-independence'
    print *, 'Solution should be constant across z-planes'
    print *, '-------------------------------------------'
    print *, ''
    
    ix = dims(1) / 2
    iy = dims(2) / 2
    
    print *, 'Sampling at ix=', ix, ', iy=', iy
    print *, ''
    print '(A10, A20, A20)', 'z-plane', 'Value', 'Diff from z=1'
    print '(A10, A20, A20)', '-------', '-----', '-------------'
    
    val_z1 = host_field%data(ix, iy, 1)
    do iz = 1, dims(3)
      val_zk = host_field%data(ix, iy, iz)
      print '(I10, 2ES20.10)', iz, val_zk, abs(val_zk - val_z1)
    end do
    
    ! =========================================================================
    ! CHECK 2: Verify analytical product formula
    ! =========================================================================
    print *, ''
    print *, '-------------------------------------------'
    print *, 'CHECK 2: Analytical Product Formula'
    print *, '-------------------------------------------'
    print *, ''
    print *, 'Verifying that:'
    print *, '  phi_110 = phi_100 * phi_010 * (-1) * kx^2 * ky^2 / (kx^2 + ky^2)'
    print *, ''
    print *, 'where:'
    print *, '  phi_100(x) = -cos(kx*x) / kx^2'
    print *, '  phi_010(y) = -cos(ky*y) / ky^2'
    print *, '  phi_110(x,y) = -cos(kx*x)*cos(ky*y) / (kx^2 + ky^2)'
    print *, ''
    print '(A8, A8, A14, A14, A14, A14, A14)', &
          'x', 'y', 'phi_100', 'phi_010', 'Product', 'phi_110', 'Diff'
    print '(A8, A8, A14, A14, A14, A14, A14)', &
          '---', '---', '-------', '-------', '-------', '-------', '----'
    
    iz = dims(3) / 2
    do iy = 1, dims(2), max(1, dims(2)/8)  ! Sample ~8 points
      do ix = 1, dims(1), max(1, dims(1)/8)
        coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        
        ! Analytical 1D solutions at these coordinates
        phi_100_analytical = -cos(kx * coord_x) / (kx * kx)
        phi_010_analytical = -cos(ky * coord_y) / (ky * ky)
        
        ! Product formula
        phi_product = phi_100_analytical * phi_010_analytical * (-1._dp) &
                      * (kx*kx * ky*ky) / (kx*kx + ky*ky)
        
        ! Direct 110 analytical
        phi_110_analytical = -cos(kx * coord_x) * cos(ky * coord_y) / (kx*kx + ky*ky)
        
        print '(2F8.4, 5ES14.4)', coord_x, coord_y, &
              phi_100_analytical, phi_010_analytical, phi_product, &
              phi_110_analytical, abs(phi_product - phi_110_analytical)
      end do
    end do
    
    ! =========================================================================
    ! CHECK 3: Numerical vs Analytical (sampled)
    ! =========================================================================
    print *, ''
    print *, '-------------------------------------------'
    print *, 'CHECK 3: Numerical vs Analytical (sampled)'
    print *, '-------------------------------------------'
    print *, ''
    
    iz = dims(3) / 2
    print *, 'Sampling at z-plane iz=', iz
    print *, ''
    print '(A8, A8, A16, A16, A16)', 'x', 'y', 'Numerical', 'Analytical', 'Error'
    print '(A8, A8, A16, A16, A16)', '---', '---', '---------', '----------', '-----'
    
    do iy = 1, dims(2), max(1, dims(2)/8)
      do ix = 1, dims(1), max(1, dims(1)/8)
        coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        
        numerical = host_field%data(ix, iy, iz)
        
        if (is_poisson_solution) then
          analytical = -cos(kx * coord_x) * cos(ky * coord_y) / (kx*kx + ky*ky)
        else
          analytical = cos(kx * coord_x) * cos(ky * coord_y)
        end if
        
        print '(2F8.4, 3ES16.6)', coord_x, coord_y, numerical, analytical, &
              abs(numerical - analytical)
      end do
    end do
    
    ! =========================================================================
    ! CHECK 4: Full error statistics (single z-plane)
    ! =========================================================================
    print *, ''
    print *, '-------------------------------------------'
    print *, 'CHECK 4: Error Statistics (z-midplane)'
    print *, '-------------------------------------------'
    print *, ''
    
    error_sum = 0._dp
    error_max = 0._dp
    count = 0
    
    iz = dims(3) / 2
    do iy = 1, dims(2)
      do ix = 1, dims(1)
        coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
        coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
        
        numerical = host_field%data(ix, iy, iz)
        
        if (is_poisson_solution) then
          analytical = -cos(kx * coord_x) * cos(ky * coord_y) / (kx*kx + ky*ky)
        else
          analytical = cos(kx * coord_x) * cos(ky * coord_y)
        end if
        
        error_sum = error_sum + (numerical - analytical)**2
        error_max = max(error_max, abs(numerical - analytical))
        count = count + 1
      end do
    end do
    
    error_l2 = sqrt(error_sum / count)
    
    print *, 'Z-plane:', iz
    print *, 'Points evaluated:', count
    print *, 'L2 error:  ', error_l2
    print *, 'Max error: ', error_max
    
    ! =========================================================================
    ! CHECK 5: Global error (all z-planes)
    ! =========================================================================
    print *, ''
    print *, '-------------------------------------------'
    print *, 'CHECK 5: Global Error (all z-planes)'
    print *, '-------------------------------------------'
    print *, ''
    
    error_sum = 0._dp
    error_max = 0._dp
    count = 0
    
    do iz = 1, dims(3)
      do iy = 1, dims(2)
        do ix = 1, dims(1)
          coord_x = self%solver%mesh%geo%midp_coords(ix, 1)
          coord_y = self%solver%mesh%geo%midp_coords(iy, 2)
          
          numerical = host_field%data(ix, iy, iz)
          
          if (is_poisson_solution) then
            analytical = -cos(kx * coord_x) * cos(ky * coord_y) / (kx*kx + ky*ky)
          else
            analytical = cos(kx * coord_x) * cos(ky * coord_y)
          end if
          
          error_sum = error_sum + (numerical - analytical)**2
          error_max = max(error_max, abs(numerical - analytical))
          count = count + 1
        end do
      end do
    end do
    
    error_l2 = sqrt(error_sum / count)
    
    print *, 'Total points evaluated:', count
    print *, 'Global L2 error:  ', error_l2
    print *, 'Global Max error: ', error_max
    
    ! =========================================================================
    ! SUMMARY
    ! =========================================================================
    print *, ''
    print *, '=========================================='
    print *, '110 CASE SUMMARY'
    print *, '=========================================='
    print *, 'L2 error:  ', error_l2
    print *, 'Max error: ', error_max
    print *, ''
    print *, 'Expected: Similar to 100 and 010 cases (~0.02-0.03)'
    print *, '=========================================='

  case default
    print *, 'BC pattern', bc_pattern, 'not implemented'
    error stop "BC pattern not implemented"

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
    call self%print_field_comparison(host_field, 'Initial cosine field:', n, .false.)

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
    call self%print_field_comparison(host_field, 'Poisson solution:', n, .true.)
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
