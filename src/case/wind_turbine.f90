module m_case_wind_turbine
  !! Case for wind-turbine simulations (the `two_turbines` setup):
  !! uniform inflow and convective outflow in x, periodic in y and z.
  !!
  !! The case owns a polymorphic turbine forcing model (turbine_model_t) and is
  !! the only place the `iturbine` switch appears (at allocation time). It drives
  !! the model from its forcings()/postprocess() hooks.
  !!
  !! The inflow/outflow boundary treatment mirrors m_case_cylinder. BCs are
  !! dirichlet/periodic/periodic for the same reason cylinder uses them: the FFT
  !! Poisson solver requires z periodic and only supports a single rank for
  !! non-periodic directions. The two_turbines case uses free-slip walls
  !! in y and z (wind-tunnel walls), which the FFT solver cannot currently do.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
  use m_adm, only: adm_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, MPI_X3D2_DP, get_argument, DIR_C, DIR_X, &
                      VERT, X_FACE, BC_DIRICHLET
  use m_config, only: wind_turbine_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init
  use m_turbine_model, only: turbine_model_t
  use m_turbine_dummy, only: turbine_dummy_t

  implicit none

  type, extends(base_case_t) :: case_wind_turbine_t
    type(wind_turbine_config_t) :: wt_cfg
    class(turbine_model_t), pointer :: turbine => null()
    real(dp) :: out_vel_cached = 0._dp
    real(dp) :: flow_rate_diff_cached = 0._dp
  contains
    procedure :: define_BC => define_BC_wind_turbine
    procedure :: initial_conditions => initial_conditions_wind_turbine
    procedure :: forcings => forcings_wind_turbine
    procedure :: apply_BC => apply_BC_wind_turbine
    procedure :: postprocess => postprocess_wind_turbine
    procedure :: finalise_case_specific => finalise_wind_turbine
    procedure :: compute_outflow_params
  end type case_wind_turbine_t

  interface case_wind_turbine_t
    module procedure case_wind_turbine_init
  end interface case_wind_turbine_t

contains

  function case_wind_turbine_init(backend, mesh, host_allocator) &
    result(flow_case)
    implicit none
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_wind_turbine_t) :: flow_case

    call flow_case%wt_cfg%read(nml_file=get_argument(1))
    call flow_case%case_init(backend, mesh, host_allocator)

    ! Select and initialise the turbine forcing model. This is the only place
    ! the iturbine switch lives; downstream the case calls the polymorphic
    ! turbine object with no branching.
    select case (flow_case%wt_cfg%iturbine)
    case (1)
      ! Actuator line model - falls back to dummy until landed.
      allocate (turbine_dummy_t :: flow_case%turbine)
    case (2)
      allocate (adm_t :: flow_case%turbine)
      select type (turbine => flow_case%turbine)
      type is (adm_t)
        turbine%coords_file = flow_case%wt_cfg%adm_coords
        turbine%rho_air = flow_case%wt_cfg%rho_air
        turbine%T_relax = flow_case%wt_cfg%T_relax
        turbine%stats_start = max(1, flow_case%io_mgr%stats_mgr%config%initstat)
        turbine%stats_freq = max(1, flow_case%io_mgr%stats_mgr%config%istatfreq)
      end select
    case (0)
      allocate (turbine_dummy_t :: flow_case%turbine)
    case default
      error stop 'Unsupported wind turbine model selected.'
    end select

    call flow_case%turbine%init(flow_case%solver%backend, &
                                flow_case%solver%mesh, &
                                flow_case%solver%host_allocator, &
                                flow_case%solver%dt)
    call flow_case%io_mgr%register_checkpoint_state(flow_case%turbine)
    call flow_case%turbine%setup_output( &
      flow_case%solver%mesh%par%is_root(), flow_case%io_mgr%is_restart() &
      )
  end function case_wind_turbine_init

  ! Initial Conditions: uniform inflow with localised noise
  subroutine initial_conditions_wind_turbine(self)
    implicit none
    class(case_wind_turbine_t) :: self

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

    noise = self%wt_cfg%init_noise

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          x = xloc(1) - self%solver%mesh%geo%L(1)/2._dp
          um = exp(-0.2_dp*x*x)

          u_init%data(i, j, k) = self%wt_cfg%bc_start_u &
                                 + noise(1)*um*(2*u_init%data(i, j, k) - 1._dp)
          v_init%data(i, j, k) = self%wt_cfg%bc_start_v &
                                 + noise(2)*um*(2*v_init%data(i, j, k) - 1._dp)
          w_init%data(i, j, k) = self%wt_cfg%bc_start_w &
                                 + noise(3)*um*(2*w_init%data(i, j, k) - 1._dp)
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
  end subroutine initial_conditions_wind_turbine

  ! Convective outflow parameters (identical treatment to m_case_cylinder).
  subroutine compute_outflow_params(self, out_vel, flow_rate_diff)
    implicit none
    class(case_wind_turbine_t) :: self
    real(dp), intent(out) :: out_vel, flow_rate_diff

    integer :: dims(3), nx, ierr
    real(dp) :: uxmax, uxmax_discard
    real(dp) :: flow_rate_in, flow_rate_out
    real(dp) :: flow_rate_in_max_discard, flow_rate_out_max_discard
    real(dp) :: fl_sums(2), ny_nz
    real(dp) :: dx, gdt

    dims = self%solver%mesh%get_dims(VERT)
    nx = dims(1)
    dx = self%solver%mesh%geo%d(1)
    ny_nz = real(dims(2)*dims(3), dp)

    gdt = self%solver%time_integrator%gdt

    call self%solver%backend%slice_max_sum( &
      uxmax, uxmax_discard, self%solver%u, nx - 1)
    call self%solver%backend%slice_max_sum( &
      flow_rate_in_max_discard, flow_rate_in, self%solver%u, 1)
    call self%solver%backend%slice_max_sum( &
      flow_rate_out_max_discard, flow_rate_out, self%solver%u, nx)

    call MPI_Allreduce(MPI_IN_PLACE, uxmax, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    fl_sums(1) = flow_rate_in
    fl_sums(2) = flow_rate_out
    call MPI_Allreduce(MPI_IN_PLACE, fl_sums, 2, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    flow_rate_in = fl_sums(1)
    flow_rate_out = fl_sums(2)

    flow_rate_in = flow_rate_in/ny_nz
    flow_rate_out = flow_rate_out/ny_nz

    out_vel = uxmax*gdt/dx
    flow_rate_diff = flow_rate_in - flow_rate_out

  end subroutine compute_outflow_params

  subroutine define_BC_wind_turbine(self)
    !! Boundary Conditions hook (called per substep before transeq).
    !! Mirrors initial_conditions structurally, but at a single x-plane:
    !!   u  = bc_start_u + inlet_noise(1) * um * (2r - 1)
    !!   v  = bc_start_v + inlet_noise(2) * um * (2r - 1)
    !!   w  = bc_start_w + inlet_noise(3) * um * (2r - 1)
    implicit none
    class(case_wind_turbine_t) :: self

    class(field_t), pointer :: hu, hv, hw
    integer :: j, k, dims(3)
    real(dp) :: noise(3), um, half_L

    dims = self%solver%mesh%get_dims(VERT)
    noise = self%wt_cfg%inlet_noise
    half_L = self%solver%mesh%geo%L(1)/2._dp
    um = exp(-0.2_dp*half_L*half_L)

    ! Sample outflow params from solver%u once per substep
    call self%compute_outflow_params(self%out_vel_cached, &
                                     self%flow_rate_diff_cached)

    ! Allocate persistent device BC fields on first call
    if (.not. associated(self%bc_start_u_x)) then
      self%bc_start_u_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_v_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_w_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      call self%bc_start_u_x%set_data_loc(VERT)
      call self%bc_start_v_x%set_data_loc(VERT)
      call self%bc_start_w_x%set_data_loc(VERT)
    end if

    ! Build the inflow profile in DIR_C host buffers, then upload
    hu => self%solver%host_allocator%get_block(DIR_C)
    hv => self%solver%host_allocator%get_block(DIR_C)
    hw => self%solver%host_allocator%get_block(DIR_C)

    ! Fill the inlet plane with random numbers in [0, 1); the loop below
    ! maps these onto noise in [-1, 1)
    call random_number(hu%data(1, 1:dims(2), 1:dims(3)))
    call random_number(hv%data(1, 1:dims(2), 1:dims(3)))
    call random_number(hw%data(1, 1:dims(2), 1:dims(3)))

    do k = 1, dims(3)
      do j = 1, dims(2)
        hu%data(1, j, k) = self%wt_cfg%bc_start_u &
                         + noise(1)*um*(2._dp*hu%data(1, j, k) - 1._dp)
        hv%data(1, j, k) = self%wt_cfg%bc_start_v &
                         + noise(2)*um*(2._dp*hv%data(1, j, k) - 1._dp)
        hw%data(1, j, k) = self%wt_cfg%bc_start_w &
                         + noise(3)*um*(2._dp*hw%data(1, j, k) - 1._dp)
      end do
    end do

    call self%solver%backend%set_field_data(self%bc_start_u_x, hu%data)
    call self%solver%backend%set_field_data(self%bc_start_v_x, hv%data)
    call self%solver%backend%set_field_data(self%bc_start_w_x, hw%data)

    call self%solver%host_allocator%release_block(hu)
    call self%solver%host_allocator%release_block(hv)
    call self%solver%host_allocator%release_block(hw)
  end subroutine define_BC_wind_turbine

  subroutine apply_BC_wind_turbine(self, u, v, w)
    !! Pre-correction (called per substep after the integrator step):
    !! Enforce inflow Dirichlet from bc_start_u_x/v/w on the inlet plane and the
    !! convective outflow update on the right face.
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: u, v, w

    call self%solver%backend%field_set_face_from_field( &
      u, self%bc_start_u_x, self%out_vel_cached, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff_cached)
    call self%solver%backend%field_set_face_from_field( &
      v, self%bc_start_v_x, self%out_vel_cached, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff_cached)
    call self%solver%backend%field_set_face_from_field( &
      w, self%bc_start_w_x, self%out_vel_cached, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff_cached)
  end subroutine apply_BC_wind_turbine

  ! Forcings: advance the turbine model and accumulate its momentum source.
  subroutine forcings_wind_turbine(self, du, dv, dw, iter)
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

    real(dp) :: t

    t = real(iter, dp)*self%solver%dt
    call self%turbine%update(t, self%solver%dt)
    call self%turbine%project_forces(du, dv, dw, &
                                     self%solver%u, self%solver%v, &
                                     self%solver%w)
  end subroutine forcings_wind_turbine

  ! Post-processing: turbine diagnostics + flow diagnostics.
  subroutine postprocess_wind_turbine(self, iter, t)
    implicit none
    class(case_wind_turbine_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%wt_cfg%iturboutput > 0) then
      if (iter > 0 .and. mod(iter, self%wt_cfg%iturboutput) == 0) then
        call self%turbine%write_output(iter, self%solver%mesh%par%is_root())
      end if
    end if

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if
    call self%monitoring%write_step( &
      self%solver, t, self%solver%u, self%solver%v, self%solver%w)
  end subroutine postprocess_wind_turbine

  subroutine finalise_wind_turbine(self)
    class(case_wind_turbine_t) :: self

    call self%turbine%finalise()
    call self%io_mgr%unregister_checkpoint_state()
    deallocate (self%turbine)
  end subroutine finalise_wind_turbine

end module m_case_wind_turbine
