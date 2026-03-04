module m_case_cylinder

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
    procedure :: apply_outflow
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
  ! This copies the field to host, updates both faces, then copies back.
  ! ==========================================================================
subroutine apply_outflow(self, fld, inflow_val)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: fld
    real(dp), intent(in) :: inflow_val

    class(field_t), pointer :: host_fld, host_u
    integer :: dims(3), nx, j, k, ierr
    real(dp) :: cfl, uxmin, uxmax, Uc
    real(dp) :: dx, dt

    dims = self%solver%mesh%get_dims(VERT)
    nx = dims(1)
    dx = self%solver%mesh%geo%d(1)
    dt = self%solver%dt

    ! --- Get u field on host to compute convection velocity ---
    host_u => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_u%data, self%solver%u)

    ! --- Compute min/max of u at nx-1 (like Xcompact3d) ---
    uxmax = -huge(1._dp)
    uxmin =  huge(1._dp)
    do k = 1, dims(3)
      do j = 1, dims(2)
        if (host_u%data(nx - 1, j, k) > uxmax) uxmax = host_u%data(nx - 1, j, k)
        if (host_u%data(nx - 1, j, k) < uxmin) uxmin = host_u%data(nx - 1, j, k)
      end do
    end do

    call self%solver%host_allocator%release_block(host_u)

    ! Global reduction across all MPI ranks
    call MPI_ALLREDUCE(MPI_IN_PLACE, uxmax, 1, MPI_X3D2_DP, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, uxmin, 1, MPI_X3D2_DP, MPI_MIN, &
                       MPI_COMM_WORLD, ierr)

    ! Convection velocity: use max velocity (safest for stability)
    Uc = uxmax

    ! CFL for the outflow advection
    cfl = Uc*dt/dx

    ! --- Copy field from device to host ---
    host_fld => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_fld%data, fld)

    ! Inflow (left face): Dirichlet
    host_fld%data(1, :, :) = inflow_val

    ! Outflow (right face):  BC
    host_fld%data(nx, :, :) = host_fld%data(nx, :, :) &
        - cfl*(host_fld%data(nx, :, :) - host_fld%data(nx - 1, :, :))

    ! Copy back to device
    call self%solver%backend%set_field_data(fld, host_fld%data)

    call self%solver%host_allocator%release_block(host_fld)

  end subroutine apply_outflow

  ! ==========================================================================
  ! Boundary Conditions: applied to U^m at the start of each substep.
  ! Inflow (i=1): u=1, v=0, w=0  (Dirichlet)
  ! Outflow (i=nx): for u, fixed v=0, w=0
  ! ==========================================================================
  subroutine boundary_conditions_cylinder(self)
    implicit none
    class(case_cylinder_t) :: self
    ! ! u: inflow = 1
    ! u, v and w: fixed on both faces
    call self%apply_outflow(self%solver%u, 1._dp)
    call self%apply_outflow(self%solver%v, 0._dp)
    call self%apply_outflow(self%solver%w, 0._dp)
    ! call self%solver%backend%field_set_face(self%solver%v, 0._dp, 0._dp, X_FACE)
    ! call self%solver%backend%field_set_face(self%solver%w, 0._dp, 0._dp, X_FACE)

  end subroutine boundary_conditions_cylinder

  ! ==========================================================================
  ! Pre-correction: enforce BCs on U* (intermediate velocity) at faces
  ! before the Poisson pressure solve.
  ! ==========================================================================
  subroutine pre_correction_cylinder(self, u, v, w)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: u, v, w


    ! ! u: inflow = 1, outflow = convective
    call self%apply_outflow(u, 1._dp)
    call self%apply_outflow(v, 0._dp)
    call self%apply_outflow(w, 0._dp)

    ! v and w: fixed on both faces
    ! call self%solver%backend%field_set_face(u, 1._dp, 1._dp, X_FACE)
    ! call self%solver%backend%field_set_face(v, 0._dp, 0._dp, X_FACE)
    ! call self%solver%backend%field_set_face(w, 0._dp, 0._dp, X_FACE)

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
