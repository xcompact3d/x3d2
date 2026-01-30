module m_case_cos2pix
  !! Test case for validating Poisson solver using pressure_correction pattern
  !!
  !! This is the CORRECT way to test the Poisson solver because it uses
  !! the operators (divergence_v2p, poisson, gradient_p2v) as designed.
  !!
  !! Test:
  !! 1. Create divergent velocity: u=0, v=cos(2*pi*y/Ly), w=0
  !! 2. Compute div(u,v,w) using divergence_v2p
  !! 3. Solve Poisson for pressure
  !! 4. Compute pressure gradient
  !! 5. Correct velocity: u_new = u - grad(p)
  !! 6. Verify: div(u_new) ≈ 0
  !!
  !! NOTE: For y-Dirichlet BCs, dims_global(2) must be ODD (e.g., 65)

  use m_allocator
  use m_base_backend
  use m_base_case, only: base_case_t
  use m_common, only: pi, dp, get_argument, DIR_C, VERT, CELL, &
                      RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_config, only: domain_config_t, solver_config_t
  use m_field, only: field_t
  use m_mesh
  use m_solver, only: solver_t
  implicit none
  
  type, extends(base_case_t) :: case_cos2pix_t
  contains
    procedure :: boundary_conditions => boundary_conditions_cos2pix
    procedure :: initial_conditions => initial_conditions_cos2pix
    procedure :: forcings => forcings_cos2pix
    procedure :: pre_correction => pre_correction_cos2pix
    procedure :: postprocess => postprocess_cos2pix
    procedure :: run => run_cos2pix
  end type case_cos2pix_t
  
  interface case_cos2pix_t
    module procedure case_cos2pix_init
  end interface case_cos2pix_t
  
contains

  function case_cos2pix_init(backend, mesh, host_allocator) result(flow_case)
    implicit none
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_cos2pix_t) :: flow_case
    call flow_case%case_init(backend, mesh, host_allocator)
  end function case_cos2pix_init

  subroutine initial_conditions_cos2pix(self)
    implicit none
    class(case_cos2pix_t) :: self
  end subroutine initial_conditions_cos2pix

  subroutine boundary_conditions_cos2pix(self)
    implicit none
    class(case_cos2pix_t) :: self
  end subroutine boundary_conditions_cos2pix

  subroutine forcings_cos2pix(self, du, dv, dw, iter)
    implicit none
    class(case_cos2pix_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter
  end subroutine forcings_cos2pix

  subroutine pre_correction_cos2pix(self, u, v, w)
    implicit none
    class(case_cos2pix_t) :: self
    class(field_t), intent(inout) :: u, v, w
  end subroutine pre_correction_cos2pix

  subroutine postprocess_cos2pix(self, iter, t)
    implicit none
    class(case_cos2pix_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t
  end subroutine postprocess_cos2pix

  subroutine run_cos2pix(self)
    implicit none
    class(case_cos2pix_t), intent(inout) :: self
    class(field_t), pointer :: u, v, w
    class(field_t), pointer :: div_u, pressure, dpdx, dpdy, dpdz
    class(field_t), pointer :: host_field
    real(dp) :: div_max, div_mean, p_max, p_mean
    integer :: i, j, k, dims(3),axis_value
    real(dp) :: xloc(3),dbg_xloc(3), L, error_norm
    character(len=1) :: axis
    class(field_t), pointer :: f
    class(field_t), pointer :: gradient_result
    class(field_t), pointer :: f_final
    class(field_t), pointer :: f_check
    class(field_t), pointer :: f_final_check
    class(field_t), pointer :: temp


    !! #TODO: 2D support for boundary conditions check
    if(self%solver%mesh%geo%BC%x(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%x(2) == 'dirichlet') then
      L = self%solver%mesh%geo%L(1)
      axis = 'x'
      axis_value = 1
    else if(self%solver%mesh%geo%BC%y(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%y(2) == 'dirichlet') then  
      L = self%solver%mesh%geo%L(2)
      axis = 'y'
      axis_value = 2
    else if(self%solver%mesh%geo%BC%z(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%z(2) == 'dirichlet') then  
      L = self%solver%mesh%geo%L(3)
      axis = 'z'
      axis_value = 3
    end if
    print *, "[axis_value] = ", axis_value 
    dims = self%solver%mesh%get_dims(CELL)
    
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '=== POISSON SOLVER VALIDATION ==='
      print *, 'Using pressure_correction pattern'
      print *, ''
      print *, 'Grid: dims =', dims
      print *, 'Domain: L', axis, ' =', L
      print *, ''
      print *, 'Test: Create divergent velocity, correct it, verify div=0'
      print *, ''
    end if


! Device field
f => self%solver%backend%allocator%get_block(DIR_C, CELL)
f_check => self%solver%backend%allocator%get_block(DIR_X)
! Host field
host_field => self%solver%host_allocator%get_block(DIR_C)

    ! Create p = cos(2*pi*[axis]/L[axis])
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k, CELL)
          ! host_field%data(i, j, k) = cos(2* pi * xloc(axis_value) / L)
          host_field%data(i, j, k) = cos(2._dp*xloc(axis_value)*pi / L)
        end do
      end do
    end do

    ! Sanity check
    print*, 'cos\n'
    if( axis_value == 1) then
      do i = 1, dims(axis_value)
        print *,  self%solver%mesh%geo%midp_coords(i, axis_value), (host_field%data(i,4,8)), (cos(2._dp * pi * self%solver%mesh%geo%midp_coords(i, axis_value))/(-4._dp * pi *pi))
      end do
    else if(axis_value == 2) then
      do j = 1, dims(axis_value)
        print *,  self%solver%mesh%geo%midp_coords(j, axis_value), (host_field%data(4,j,8)), (cos(2._dp * pi * self%solver%mesh%geo%midp_coords(j, axis_value))/(-4._dp * pi *pi))
      end do
    else
      error stop "not implemented yet!!"
    end if
    ! Assign host_field to f and delete host_field
    print *, "1"
    call self%solver%backend%set_field_data(f, host_field%data, DIR_C)
    print *, "2"
    call f%set_data_loc(CELL)
    print *, "3"
    call self%solver%host_allocator%release_block(host_field)
print *, "4"
    call self%solver%backend%reorder(f_check, f, RDR_C2X)
print *, "5"
    ! Step 1: Solve Poisson  
    temp => self%solver%backend%allocator%get_block(DIR_C)
    print *, "6"
    ! solve poisson equation with FFT based approach
    call self%solver%backend%poisson_fft%solve_poisson(f, temp)
    print *, "7"
    call self%solver%backend%allocator%release_block(temp)
    print *, "8"
    
    host_field => self%solver%host_allocator%get_block(DIR_C)
    print *, "9"
    call self%solver%backend%get_field_data(host_field%data, f)
    print*, 'last\n'    

    if( axis_value == 1) then
      do i = 1, dims(axis_value)
        print *,  self%solver%mesh%geo%midp_coords(i, axis_value), (host_field%data(i,4,8)), (cos(2._dp * pi * self%solver%mesh%geo%midp_coords(i, axis_value))/(-4._dp * pi *pi))
      end do
    else if(axis_value == 2) then
      do j = 1, dims(axis_value)
        print *,  self%solver%mesh%geo%midp_coords(j, axis_value), (host_field%data(4,j,8)), (cos(2._dp * pi * self%solver%mesh%geo%midp_coords(j, axis_value))/(-4._dp * pi *pi))
      end do
    else
      error stop "not implemented yet!!"
    end if
    call self%solver%host_allocator%release_block(host_field)
    ! Step 2: Compute gradient of pressure
    gradient_result => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%backend%reorder(gradient_result, f, RDR_C2Z)
    call self%solver%backend%allocator%release_block(f)

    dpdx  => self%solver%backend%allocator%get_block(DIR_X)
    dpdy  => self%solver%backend%allocator%get_block(DIR_X)
    dpdz  => self%solver%backend%allocator%get_block(DIR_X)

    
    call self%solver%gradient_p2v(dpdx, dpdy, dpdz,gradient_result)
    
    call self%solver%backend%allocator%release_block(gradient_result)

    f_final => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%divergence_v2p(f_final, dpdx, dpdy, dpdz)

    call self%solver%backend%allocator%release_block(dpdx)
    call self%solver%backend%allocator%release_block(dpdy)
    call self%solver%backend%allocator%release_block(dpdz)
    

  f_final_check => self%solver%backend%allocator%get_block(DIR_X)
  call self%solver%backend%reorder(f_final_check, f_final, RDR_Z2X)
  call self%solver%backend%allocator%release_block(f_final)

  call self%solver%backend%vecadd(-1._dp, f_check, 1._dp, f_final_check)

  host_field => self%solver%host_allocator%get_block(DIR_C)

  call self%solver%backend%get_field_data(host_field%data, f_final_check)
  call self%solver%backend%allocator%release_block(f_final_check)
  dims = self%solver%mesh%get_dims(CELL)
  error_norm = norm2(host_field%data(1:dims(1), 1:dims(2), 1:dims(3)))/product(dims)
  print*, 'norm', error_norm

    call self%solver%backend%allocator%release_block(f_check)
  call self%solver%host_allocator%release_block(host_field)


  if (error_norm .gt. 1d-9) then
    error stop 'div(grad(fft())) test failed!'
  else
    if (self%solver%mesh%par%is_root()) print *, "TEST PASS"
  end if


    call self%case_finalise()
    
  end subroutine run_cos2pix

end module m_case_cos2pix
