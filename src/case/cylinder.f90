module m_case_cylinder
  !! Cylinder flow case using IBM (Immersed Boundary Method).
  !!
  !! Algorithm (projection / fractional-step method):
  !!   1. Define initial conditions for velocity (uniform + noise)
  !!   TIME LOOP (U^m -> U^{m+1}):
  !!     2. Store BCs in 2D fields (faces) using velocity U^m
  !!     3. Convective-diffusive (transport) terms computed via transeq
  !!     4. Forcings (IBM handled by solver when ibm_on=.true.)
  !!     5. Time integration  U^m -> U*
  !!     6. Apply BC via pre-correction on velocity U* (using faces)
  !!     7. IBM body forcing on U* (handled by solver)
  !!     8. Divergence of U* (intermediate velocity)
  !!     9. Poisson equation to get pressure: nabla^2 p = div(U*)
  !!    10. Pressure gradients (cell -> face)
  !!    11. Correction: U^{m+1} = U* - nabla p
  !!    12. Check divergence ~ zero (machine precision)
  !!     - Post-processing
  !!   END TIME LOOP
  !!
  !! Outflow BC uses the advective (convective) condition:
  !!   du/dt + Uc * du/dx = 0
  !! Discretised as:
  !!   u_N^{n+1} = u_N^n - Uc * dt/dx * (u_N^n - u_{N-1}^n)
  !!
  !! Note: BCs on U* must be applied BEFORE solving the Poisson equation.
  !! The cylinder is represented via the solver's built-in IBM (ibm_on=.true.).

  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, MPI_X3D2_DP, get_argument, DIR_C, VERT, CELL, X_FACE
  use m_field, only: field_t
  use m_config, only: cylinder_config_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  real(dp), parameter :: U_conv = 1._dp  !! Convection velocity for outflow BC

  type, extends(base_case_t) :: case_cylinder_t
    type(cylinder_config_t) :: cylinder_cfg
  contains
    procedure :: boundary_conditions => boundary_conditions_cylinder
    procedure :: initial_conditions => initial_conditions_cylinder
    procedure :: forcings => forcings_cylinder
    procedure :: pre_correction => pre_correction_cylinder
    procedure :: postprocess => postprocess_cylinder
    procedure :: apply_outflow_convective
  end type case_cylinder_t

  interface case_cylinder_t
    module procedure case_cylinder_init
  end interface case_cylinder_t

contains

  ! ==========================================================================
  ! Initialization
  ! ==========================================================================
  function case_cylinder_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_cylinder_t) :: flow_case

    call flow_case%cylinder_cfg%read(nml_file=get_argument(1))

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_cylinder_init

  ! ==========================================================================
  ! Initial Conditions: uniform flow u=1, v=0, w=0 with localized noise
  ! to help trigger vortex shedding behind the cylinder.
  ! ==========================================================================
  subroutine initial_conditions_cylinder(self)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), pointer :: u_init, v_init, w_init

    integer :: i, j, k, dims(3)
    real(dp) :: xloc(3), x, noise(3), um

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

    ! Fill with random numbers first
    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(v_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(w_init%data(1:dims(1), 1:dims(2), 1:dims(3)))

    noise = self%cylinder_cfg%init_noise
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          ! Gaussian envelope centred at domain midpoint in x
          x = xloc(1) - self%solver%mesh%geo%L(1)/2._dp
          um = exp(-0.2_dp*x*x)

          ! Uniform flow + localized noise perturbation (per component)
          u_init%data(i, j, k) = 1._dp &
                                 + noise(1)*um*(2*u_init%data(i, j, k) - 1._dp)
          v_init%data(i, j, k) = noise(2)*um*(2*v_init%data(i, j, k) - 1._dp)
          w_init%data(i, j, k) = noise(3)*um*(2*w_init%data(i, j, k) - 1._dp)
        end do
      end do
    end do

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_cylinder

  ! ==========================================================================
  ! Convective outflow BC for a single field:
  !   du/dt + Uc * du/dx = 0
  ! Discretised (backward Euler in time, backward difference in space):
  !   u_N^{n+1} = u_N^n - Uc * dt/dx * (u_N^n - u_{N-1}^n)
  !
  ! Also sets the inflow (left) face to the specified value.
  !
  ! This copies the field to host, updates both faces, then copies back.
  ! ==========================================================================
  subroutine apply_outflow_convective(self, fld, inflow_val)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: fld
    real(dp), intent(in) :: inflow_val

    class(field_t), pointer :: host_fld
    integer :: dims(3), nx
    real(dp) :: cfl  ! Uc * dt / dx

    dims = self%solver%mesh%get_dims(VERT)
    nx = dims(1)

    ! Courant number for the outflow advection
    cfl = U_conv*self%solver%dt/self%solver%mesh%geo%d(1)

    ! Copy field from device to host
    host_fld => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_fld%data, fld)

    ! Inflow (left face): set to prescribed value
    host_fld%data(1, :, :) = inflow_val

    ! Outflow (right face): advective BC
    ! u_N^{n+1} = u_N^n - cfl * (u_N^n - u_{N-1}^n)
    host_fld%data(nx, :, :) = host_fld%data(nx, :, :) &
                  - cfl*(host_fld%data(nx, :, :) - host_fld%data(nx - 1, :, :))

    ! Copy back to device
    call self%solver%backend%set_field_data(fld, host_fld%data)

    call self%solver%host_allocator%release_block(host_fld)

  end subroutine apply_outflow_convective

  ! ==========================================================================
  ! Boundary Conditions: applied to U^m at the start of each substep.
  ! Inflow (i=1): u=1, v=0, w=0  (Dirichlet)
  ! Outflow (i=nx): convective BC for u, fixed v=0, w=0
  ! ==========================================================================
  subroutine boundary_conditions_cylinder(self)
    implicit none
    class(case_cylinder_t) :: self
    ! ! u: inflow = 1, outflow = convective
    ! call self%apply_outflow_convective(self%solver%u, 1._dp)
    ! u, v and w: fixed on both faces
    call self%solver%backend%field_set_face(self%solver%u, 1._dp, 1._dp, X_FACE)
    call self%solver%backend%field_set_face(self%solver%v, 0._dp, 0._dp, X_FACE)
    call self%solver%backend%field_set_face(self%solver%w, 0._dp, 0._dp, X_FACE)

  end subroutine boundary_conditions_cylinder

  ! ==========================================================================
  ! Pre-correction: enforce BCs on U* (intermediate velocity) at faces
  ! before the Poisson pressure solve.
  ! This is critical — BCs on U* must be applied BEFORE nabla^2 p = div(U*).
  ! ==========================================================================
  subroutine pre_correction_cylinder(self, u, v, w)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: u, v, w


    ! ! u: inflow = 1, outflow = convective
    ! call self%apply_outflow_convective(u, 1._dp)

    ! v and w: fixed on both faces
    call self%solver%backend%field_set_face(u, 1._dp, 1._dp, X_FACE)
    call self%solver%backend%field_set_face(v, 0._dp, 0._dp, X_FACE)
    call self%solver%backend%field_set_face(w, 0._dp, 0._dp, X_FACE)

  end subroutine pre_correction_cylinder

  ! ==========================================================================
  ! Forcings: empty — cylinder forcing is handled by the solver's IBM
  ! (self%solver%ibm%body is called in the run loop when ibm_on=.true.)
  ! ==========================================================================
  subroutine forcings_cylinder(self, du, dv, dw, iter)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter
    ! No additional forcing terms — IBM is handled by solver
  end subroutine forcings_cylinder

  ! ==========================================================================
  ! Post-processing: report enstrophy and divergence
  ! ==========================================================================
  subroutine postprocess_cylinder(self, iter, t)
    implicit none
    class(case_cylinder_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_cylinder

end module m_case_cylinder
