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

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, VERT, DIR_C, DIR_X, DIR_Z, pi
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init
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
    integer :: i, j, k, dims(3)
    real(dp) :: xloc(3),dbg_xloc(3), axis_value, L
    character(len=1) :: axis

    !! #TODO: 2D support for boundary conditions check
    if(self%solver%mesh%geo%BC%x(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%x(2) == 'dirichlet') then
      L = self%solver%mesh%geo%L(1)
      axis = 'x'
      axis_value = 0
    else if(self%solver%mesh%geo%BC%y(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%y(2) == 'dirichlet') then  
      L = self%solver%mesh%geo%L(2)
      axis = 'y'
      axis_value = 1
    else if(self%solver%mesh%geo%BC%z(1) == 'dirichlet' .and. self%solver%mesh%geo%BC%z(2) == 'dirichlet') then  
      L = self%solver%mesh%geo%L(3)
      axis = 'z'
      axis_value = 2
    end if

    dims = self%solver%mesh%get_dims(VERT)
    
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
    
    host_field => self%solver%host_allocator%get_block(DIR_C)
    
    ! Create u = 0
    u => self%solver%backend%allocator%get_block(DIR_X)
    call u%fill(0._dp)


    ! Create v = cos(2*pi*[axis]/L[axis])
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          host_field%data(i, j, k) = cos(2* pi * xloc(axis_value) / L)
        end do
      end do
    end do
    v => self%solver%backend%allocator%get_block(DIR_X)
    call self%solver%backend%set_field_data(v, host_field%data)
    
    ! Create w = 0
    w => self%solver%backend%allocator%get_block(DIR_X)
    call w%fill(0._dp)
    
    call self%solver%host_allocator%release_block(host_field)
    
    ! Step 1: Compute divergence
    div_u => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%divergence_v2p(div_u, u, v, w)
    
    call div_u%set_data_loc(VERT)
    call self%solver%backend%field_max_mean(div_max, div_mean, div_u)
    if (self%solver%mesh%par%is_root()) then
      print *, 'Step 1: Initial divergence'
      print *, '  max|div(u,v,w)| =', div_max
      print *, '  mean =', div_mean
      print *, '  (expected max: pi/L', axis,' =', pi/L, ')'
      print *, ''
    end if
    
    ! Step 2: Solve Poisson
    pressure => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%poisson(pressure, div_u)
    call self%solver%backend%allocator%release_block(div_u)
    
    call pressure%set_data_loc(VERT)
    call self%solver%backend%field_max_mean(p_max, p_mean, pressure)
    if (self%solver%mesh%par%is_root()) then
      print *, 'Step 2: Pressure from Poisson solve'
      print *, '  max|p| =', p_max
      print *, '  mean(p) =', p_mean
      print *, ''
    end if
    
    ! Step 3: Compute gradient
    dpdx => self%solver%backend%allocator%get_block(DIR_X)
    dpdy => self%solver%backend%allocator%get_block(DIR_X)
    dpdz => self%solver%backend%allocator%get_block(DIR_X)
    call self%solver%gradient_p2v(dpdx, dpdy, dpdz, pressure)
    call self%solver%backend%allocator%release_block(pressure)
    
    if (self%solver%mesh%par%is_root()) then
      print *, 'Step 3: Computed pressure gradient'
    end if
    
    ! Step 4: Correct velocity
    call self%solver%backend%vecadd(-1._dp, dpdx, 1._dp, u)
    call self%solver%backend%vecadd(-1._dp, dpdy, 1._dp, v)
    call self%solver%backend%vecadd(-1._dp, dpdz, 1._dp, w)
    
    call self%solver%backend%allocator%release_block(dpdx)
    call self%solver%backend%allocator%release_block(dpdy)
    call self%solver%backend%allocator%release_block(dpdz)
    
    if (self%solver%mesh%par%is_root()) then
      print *, 'Step 4: Corrected velocity'
    end if
    
    ! Step 5: Compute divergence of corrected velocity
    div_u => self%solver%backend%allocator%get_block(DIR_Z)
    call self%solver%divergence_v2p(div_u, u, v, w)
    
    call div_u%set_data_loc(VERT)
    call self%solver%backend%field_max_mean(div_max, div_mean, div_u)
    
    if (self%solver%mesh%par%is_root()) then
      print *, ''
      print *, '=== RESULTS ==='
      print *, 'Divergence of corrected velocity:'
      print *, '  max|div| =', div_max
      print *, '  mean|div| =', div_mean
      print *, ''
      if (div_max < 1.0e-10_dp) then
        print *, 'STATUS: TEST PASS! (div < 1e-10)'
      else
        print *, 'STATUS: TEST FAILED'
      end if
      print *, ''
    end if
    
    call self%solver%backend%allocator%release_block(u)
    call self%solver%backend%allocator%release_block(v)
    call self%solver%backend%allocator%release_block(w)
    call self%solver%backend%allocator%release_block(div_u)
    call self%case_finalise()
    
  end subroutine run_cos2pix

end module m_case_cos2pix
