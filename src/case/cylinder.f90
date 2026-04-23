module m_case_cylinder
  use iso_fortran_env, only: stderr => error_unit
  use mpi
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, MPI_X3D2_DP, get_argument, DIR_C, DIR_X, &
                       VERT, CELL, X_FACE, BC_DIRICHLET
  use m_field, only: field_t
  use m_config, only: cylinder_config_t
  use m_mesh, only: mesh_t
  use m_solver, only: init
  implicit none

  type, extends(base_case_t) :: case_cylinder_t
    type(cylinder_config_t) :: cylinder_cfg
  contains
    procedure :: boundary_conditions => boundary_conditions_cylinder
    procedure :: initial_conditions => initial_conditions_cylinder
    procedure :: forcings => forcings_cylinder
    procedure :: pre_correction => pre_correction_cylinder
    procedure :: postprocess => postprocess_cylinder
    procedure :: compute_outflow_params
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

    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(v_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(w_init%data(1:dims(1), 1:dims(2), 1:dims(3)))

    noise = self%cylinder_cfg%init_noise

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          x = xloc(1) - self%solver%mesh%geo%L(1)/2._dp
          um = exp(-0.2_dp*x*x)

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
  ! Compute outflow velocity number.
  !
  ! This is the ONE place we still touch the host: we copy only the u field
  ! to extract max(u) at the penultimate x-station, then do an MPI reduction.
  ! The result is a single scalar velocity passed to the GPU kernels.
  !
  ! TODO: Replace with a device-side reduction (cub::DeviceReduce or a
  !       custom kernel over a single x-slice) to eliminate this transfer.
  ! ==========================================================================
  subroutine compute_outflow_params(self, out_vel)
    implicit none
    class(case_cylinder_t) :: self
    real(dp) :: out_vel

    class(field_t), pointer :: host_u
    integer :: dims(3), nx, j, k, ierr
    real(dp) :: uxmax, dx, dt

    dims = self%solver%mesh%get_dims(VERT)
    nx   = dims(1)
    dx   = self%solver%mesh%geo%d(1)
    dt   = self%solver%dt

    ! Copy u field to host
    host_u => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(host_u%data, self%solver%u)

    ! Local max of u at x = nx-1
    uxmax = -huge(1._dp)
    do k = 1, dims(3)
      do j = 1, dims(2)
        if (host_u%data(nx - 1, j, k) > uxmax) &
          uxmax = host_u%data(nx - 1, j, k)
      end do
    end do

    call self%solver%host_allocator%release_block(host_u)

    ! Global max across all MPI ranks
    call MPI_ALLREDUCE(MPI_IN_PLACE, uxmax, 1, MPI_X3D2_DP, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)

    out_vel = uxmax * dt / dx
  end subroutine compute_outflow_params

  ! ==========================================================================
  ! Boundary Conditions: applied to U^m at the start of each substep.
  !
  ! Inflow  (left,  i=1):  Dirichlet  u=1, v=0, w=0
  ! Outflow (right, i=nx): Convective du/dt + Uc*du/dx = 0
  !
  ! Everything runs on GPU via field_set_face — no per-field host round-trips.
  ! ==========================================================================
  subroutine boundary_conditions_cylinder(self)
    implicit none
    class(case_cylinder_t) :: self
    real(dp) :: out_vel

    call self%compute_outflow_params(out_vel)

    call self%solver%backend%field_set_face( &
        self%solver%u, 1._dp, out_vel, X_FACE)

    call self%solver%backend%field_set_face( &
        self%solver%v, 0._dp, out_vel, X_FACE)

    call self%solver%backend%field_set_face( &
        self%solver%w, 0._dp, out_vel, X_FACE) 
  end subroutine boundary_conditions_cylinder

  ! ==========================================================================
  ! Pre-correction: enforce BCs on U* before the pressure Poisson solve.
  ! Same logic as boundary_conditions — Dirichlet inflow, convective outflow.
  ! ==========================================================================
  subroutine pre_correction_cylinder(self, u, v, w)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: u, v, w
    real(dp) :: out_vel

    call self%compute_outflow_params(out_vel)

    call self%solver%backend%field_set_face( &
        u, 1._dp, out_vel, X_FACE)

    call self%solver%backend%field_set_face( &
        v, 0._dp, out_vel, X_FACE)

    call self%solver%backend%field_set_face( &
        w, 0._dp, out_vel, X_FACE) 
  end subroutine pre_correction_cylinder

  ! ==========================================================================
  ! Forcings: empty — cylinder forcing is handled by the solver's IBM
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