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
  !! non-periodic directions. The legacy two_turbines case uses free-slip walls
  !! in y and z (wind-tunnel walls), which the FFT solver cannot currently do.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
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
    class(turbine_model_t), allocatable :: turbine
    real(dp) :: out_vel_cached = 0._dp
    real(dp) :: flow_rate_diff_cached = 0._dp
    logical :: outflow_params_valid = .false.
  contains
    procedure :: define_BC => define_BC_wind_turbine
    procedure :: initial_conditions => initial_conditions_wind_turbine
    procedure :: forcings => forcings_wind_turbine
    procedure :: apply_BC => apply_BC_wind_turbine
    procedure :: postprocess => postprocess_wind_turbine
    procedure :: compute_outflow_params
    procedure :: apply_outflow_bc
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
      ! Actuator disc model - falls back to dummy until landed.
      allocate (turbine_dummy_t :: flow_case%turbine)
    case default
      ! iturbine = 0: no turbine forcing.
      allocate (turbine_dummy_t :: flow_case%turbine)
    end select

    call flow_case%turbine%init(flow_case%solver%backend, &
                                flow_case%solver%mesh, &
                                flow_case%solver%host_allocator, &
                                flow_case%solver%dt)
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

  subroutine apply_outflow_bc(self, u, v, w)
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: u, v, w
    real(dp) :: out_vel, flow_rate_diff

    if (self%outflow_params_valid) then
      out_vel = self%out_vel_cached
      flow_rate_diff = self%flow_rate_diff_cached
      self%outflow_params_valid = .false.
    else
      call self%compute_outflow_params(out_vel, flow_rate_diff)
      self%out_vel_cached = out_vel
      self%flow_rate_diff_cached = flow_rate_diff
    end if
    associate (cfg => self%wt_cfg)
      call self%solver%backend%field_set_face( &
        u, cfg%bc_start_u, out_vel, X_FACE, &
        bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
        flow_rate_diff=flow_rate_diff)
      call self%solver%backend%field_set_face( &
        v, cfg%bc_start_v, out_vel, X_FACE, &
        bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
        flow_rate_diff=flow_rate_diff)
      call self%solver%backend%field_set_face( &
        w, cfg%bc_start_w, out_vel, X_FACE, &
        bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
        flow_rate_diff=flow_rate_diff)
    end associate
  end subroutine apply_outflow_bc

  subroutine define_BC_wind_turbine(self)
    implicit none
    class(case_wind_turbine_t) :: self
    call self%apply_outflow_bc(self%solver%u, self%solver%v, self%solver%w)
    self%outflow_params_valid = .true.
  end subroutine define_BC_wind_turbine

  subroutine apply_BC_wind_turbine(self, u, v, w)
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: u, v, w
    call self%apply_outflow_bc(u, v, w)
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

    if (mod(iter, self%wt_cfg%iturboutput) == 0) then
      call self%turbine%write_output(iter, self%solver%mesh%par%is_root())
    end if

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if
    call self%monitoring%write_step( &
      self%solver, t, self%solver%u, self%solver%v, self%solver%w)
  end subroutine postprocess_wind_turbine

end module m_case_wind_turbine
