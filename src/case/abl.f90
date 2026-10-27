module m_case_abl
  use m_allocator, only: allocator_t
  use m_abl, only: abl_t
  use m_adm, only: adm_t
  use m_abl_diagnostics, only: abl_diagnostics_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, MPI_X3D2_DP, DIR_X, VERT, &
                      X_FACE, Y_FACE, BC_DIRICHLET, BC_NEUMANN, BC_PERIODIC
  use m_config, only: abl_config_t, inflow_config_t, wind_turbine_config_t
  use m_field, only: field_t
  use m_inflow_planes, only: inflow_plane_writer_t, inflow_plane_reader_t
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce
  use m_logger, only: log_info, abort_run, str
  use m_turbine_model, only: turbine_model_t

  implicit none

  type, extends(base_case_t) :: case_abl_t
    type(abl_config_t) :: abl_cfg
    type(abl_t) :: abl
    type(abl_diagnostics_t) :: diagnostics
    type(inflow_config_t) :: inflow_cfg
    !> Precursor output: inflow planes for a successor run
    type(inflow_plane_writer_t) :: inflow_writer
    !> Successor input: inflow planes from a precursor run
    type(inflow_plane_reader_t) :: inflow_reader
    !> Iteration the inlet plane in bc_start_*_x was built for
    integer :: inflow_iter = -1
    real(dp) :: out_vel = 0._dp
    !> y weights of the inlet and outlet flow rates (outflow correction)
    real(dp), allocatable :: flow_rate_y_weights(:)
    !> Turbine settings from wind_turbine_nml, which is optional here
    type(wind_turbine_config_t) :: turbine_cfg
    !> Turbine forcing model; only allocated when wind_turbine_nml asks for
    !> one
    class(turbine_model_t), pointer :: turbine => null()
  contains
    procedure :: define_BC => define_BC_abl
    procedure :: initial_conditions => initial_conditions_abl
    procedure :: forcings => forcings_abl
    procedure :: apply_BC => apply_BC_abl
    procedure :: postprocess => postprocess_abl
    procedure :: post_step => post_step_abl
    procedure :: finalise_case_specific => finalise_abl
    procedure :: update_inflow_plane
    procedure :: init_turbine
  end type case_abl_t

  interface case_abl_t
    module procedure case_abl_init
  end interface case_abl_t

contains

  function case_abl_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_abl_t) :: flow_case

    integer :: bcs(3, 2), bc_x, dims(3)

    call flow_case%abl_cfg%read(nml_file=get_argument(1))
    call flow_case%inflow_cfg%read(nml_file=get_argument(1))
    ! Turbines are set in wind_turbine_nml, as for the wind_turbine case,
    ! and are optional here. Its uniform-inflow keys belong to that case, so
    ! setting one is an error rather than silently ignored.
    flow_case%turbine_cfg%required = .false.
    call flow_case%turbine_cfg%read(nml_file=get_argument(1))
    if (any(flow_case%turbine_cfg%init_noise /= 0._dp) &
        .or. any(flow_case%turbine_cfg%inlet_noise /= 0._dp) &
        .or. flow_case%turbine_cfg%bc_start_u /= 1._dp &
        .or. flow_case%turbine_cfg%bc_start_v /= 0._dp &
        .or. flow_case%turbine_cfg%bc_start_w /= 0._dp) &
      call abort_run('ABL: wind_turbine_nml init_noise, inlet_noise and &
                     &bc_start_* belong to the wind_turbine case; the ABL &
                     &takes its inflow from abl_nml and inflow_nml.', &
                     collective=.true.)

    ! apply_BC_abl and the wall model hard-code these boundaries, so the
    ! mesh must agree with them for the pressure solve and precorrect_walls.
    ! A successor reading precursor planes has an x inflow and outflow.
    bc_x = BC_PERIODIC
    if (flow_case%inflow_cfg%read_planes) bc_x = BC_DIRICHLET
    bcs = mesh%grid%BCs_global
    if (any(bcs(1, :) /= bc_x) .or. any(bcs(3, :) /= BC_PERIODIC) &
        .or. bcs(2, 1) /= BC_DIRICHLET .or. bcs(2, 2) /= BC_NEUMANN) &
      call abort_run('ABL requires BC_x periodic (dirichlet when read_planes), &
                     &periodic BC_z, and BC_y = ''dirichlet'', ''neumann''.', &
                     collective=.true.)

    ! Damping is a sponge rather than a driving mechanism; a successor is
    ! driven by its inflow.
    if (.not. (flow_case%abl_cfg%pressure_gradient .or. &
               flow_case%abl_cfg%coriolis .or. &
               flow_case%abl_cfg%mass_conserve .or. &
               flow_case%inflow_cfg%read_planes)) &
      call abort_run('ABL config error: enable pressure_gradient, coriolis, &
                     &or mass_conserve, or read inflow planes.', &
                     collective=.true.)
    if (flow_case%abl_cfg%shifted_periodic .and. &
        flow_case%inflow_cfg%read_planes) &
      call abort_run('ABL shifted_periodic needs a periodic x direction; it &
                     &cannot be used with read_planes.', collective=.true.)
    flow_case%abl = abl_t(backend, mesh, host_allocator, flow_case%abl_cfg)
    if (flow_case%abl_cfg%shifted_periodic) &
      call flow_case%abl%configure_shifted_periodic()

    call flow_case%case_init(backend, mesh, host_allocator)
    flow_case%solver%keep_wall_correction = .true.
    call flow_case%abl%configure_wall_boundary_correction(flow_case%solver%les)
    flow_case%diagnostics = abl_diagnostics_t(backend, mesh, &
                                              flow_case%abl_cfg)
    flow_case%inflow_writer = inflow_plane_writer_t( &
                              backend, mesh, flow_case%inflow_cfg, &
                              flow_case%solver%dt)
    flow_case%inflow_reader = inflow_plane_reader_t(mesh, flow_case%inflow_cfg)

    ! Flow rate through an x face as the projection sees it: the y
    ! interpolation gives the floor and lid planes half weight. The floor
    ! is then left out, since apply_BC_abl stamps it to zero on both faces
    ! after the outflow correction. y is not decomposed (the wall model
    ! requires it), so the local planes are the global ones.
    if (flow_case%inflow_cfg%read_planes) then
      dims = mesh%get_dims(VERT)
      allocate (flow_case%flow_rate_y_weights(dims(2)))
      flow_case%flow_rate_y_weights = 1._dp
      flow_case%flow_rate_y_weights(1) = 0._dp
      flow_case%flow_rate_y_weights(dims(2)) = 0.5_dp
    end if

    if (flow_case%turbine_cfg%iturbine /= 0) call flow_case%init_turbine()

  end function case_abl_init

  subroutine init_turbine(self)
    !! Turbines in the boundary layer, typically in a successor driven by a
    !! precursor's inflow (Incompact3d iturbine with an ABL itype). They are
    !! set in wind_turbine_nml, as for the wind_turbine case.
    implicit none

    class(case_abl_t) :: self

    select case (self%turbine_cfg%iturbine)
    case (2)
      allocate (adm_t :: self%turbine)
      select type (turbine => self%turbine)
      type is (adm_t)
        turbine%coords_file = self%turbine_cfg%adm_coords
        turbine%rho_air = self%turbine_cfg%rho_air
        turbine%T_relax = self%turbine_cfg%T_relax
        turbine%stats_start = max(1, self%io_mgr%stats_mgr%config%initstat)
        turbine%stats_freq = max(1, self%io_mgr%stats_mgr%config%istatfreq)
      end select
    case default
      call abort_run('ABL: wind_turbine_nml iturbine must be 0 (none) or 2 &
                     &(actuator disc).', collective=.true.)
    end select

    call self%turbine%init(self%solver%backend, self%solver%mesh, &
                           self%solver%host_allocator, self%solver%dt)
    call self%io_mgr%register_checkpoint_state(self%turbine)
    call self%turbine%setup_output(self%solver%mesh%par%is_root(), &
                                   self%io_mgr%is_restart())
  end subroutine init_turbine

  subroutine initial_conditions_abl(self)
    implicit none

    class(case_abl_t) :: self

    call self%abl%initialise(self%solver%u, self%solver%v, self%solver%w)

  end subroutine initial_conditions_abl

  subroutine define_BC_abl(self)
    implicit none

    class(case_abl_t) :: self

    real(dp) :: ub, target_mean, can, ly
    real(dp), allocatable :: sums(:)
    integer :: dims(3), global_dims(3), j, ierr

    ! Constant-flow-rate correction (Incompact3d forceabl); mirrors the channel
    ! bulk-velocity shift, targeting the log-law flow rate. The wall stress
    ! is not applied here: the wall model supplies it to the SGS stress.
    if (self%abl_cfg%mass_conserve) then
      ly = self%solver%mesh%geo%L(2)
      ! Bulk velocity as forceabl computes it: the plane mean of u, integrated
      ! over y with the trapezoidal rule on the vertices. A plain vertex sum
      ! would give the free-slip lid full weight. y is not decomposed (see
      ! configure_wall_boundary_correction), so each rank holds every plane.
      dims = self%solver%mesh%get_dims(VERT)
      global_dims = self%solver%mesh%get_global_dims(VERT)
      allocate (sums(dims(2)))
      call self%solver%backend%field_plane_sums(sums, self%solver%u)
      call MPI_Allreduce(MPI_IN_PLACE, sums, dims(2), MPI_X3D2_DP, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      sums = sums/real(global_dims(1)*global_dims(3), dp)
      ub = 0._dp
      do j = 1, dims(2) - 1
        ub = ub + 0.5_dp*(sums(j) + sums(j + 1)) &
             *(self%solver%mesh%geo%vert_coords(j + 1, 2) &
               - self%solver%mesh%geo%vert_coords(j, 2))
      end do
      ub = ub/ly
      if (self%abl_cfg%u_bulk > 0._dp) then
        target_mean = self%abl_cfg%u_bulk
      else
        target_mean = self%abl_cfg%u_star/self%abl_cfg%kappa &
                      *(ly*log(self%abl_cfg%delta/self%abl_cfg%z0) &
                        - self%abl_cfg%delta)/ly
      end if
      can = target_mean - ub
      call self%solver%backend%field_shift(self%solver%u, can)
    end if

    ! After the flow-rate correction, as in Incompact3d
    ! boundary_conditions_abl
    if (self%abl_cfg%shifted_periodic) &
      call self%abl%apply_shifted_periodic(self%solver%u, self%solver%v, &
                                           self%solver%w)

    ! The flow-rate correction is left to apply_BC, which sees the field
    ! the stamp is applied to
    if (self%inflow_cfg%read_planes) then
      call self%compute_outflow_params(self%out_vel)
      call self%update_inflow_plane()
    end if

  end subroutine define_BC_abl

  subroutine update_inflow_plane(self)
    !! Fill bc_start_*_x with the precursor plane for the iteration being
    !! computed. The plane is held for every substep of the iteration, as
    !! Incompact3d does, so it is uploaded once per iteration.
    implicit none

    class(case_abl_t) :: self

    real(dp), allocatable :: pu(:, :), pv(:, :), pw(:, :)
    integer :: iter, dims(3)

    iter = self%solver%current_iter + 1
    if (iter == self%inflow_iter) return
    self%inflow_iter = iter

    if (.not. associated(self%bc_start_u_x)) then
      self%bc_start_u_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_v_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_w_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      call self%bc_start_u_x%set_data_loc(VERT)
      call self%bc_start_v_x%set_data_loc(VERT)
      call self%bc_start_w_x%set_data_loc(VERT)
    end if

    ! Only the x = 1 plane of the BC fields is read, by the face stamp and
    ! the outflow correction, so only that plane is uploaded
    dims = self%solver%mesh%get_dims(VERT)
    allocate (pu(dims(2), dims(3)), pv(dims(2), dims(3)), pw(dims(2), dims(3)))
    call self%inflow_reader%planes_at(iter*self%solver%dt, pu, pv, pw)
    call self%solver%backend%field_set_x_inlet_plane(self%bc_start_u_x, pu)
    call self%solver%backend%field_set_x_inlet_plane(self%bc_start_v_x, pv)
    call self%solver%backend%field_set_x_inlet_plane(self%bc_start_w_x, pw)
  end subroutine update_inflow_plane

  subroutine forcings_abl(self, du, dv, dw, iter)
    implicit none

    class(case_abl_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

    call self%abl%apply_forcing(du, dv, dw, &
                                self%solver%u, self%solver%v, self%solver%w)

    if (associated(self%turbine)) then
      call self%turbine%update(real(iter, dp)*self%solver%dt, self%solver%dt)
      call self%turbine%project_forces(du, dv, dw, self%solver%u, &
                                       self%solver%v, self%solver%w)
    end if

  end subroutine forcings_abl

  subroutine apply_BC_abl(self, u, v, w)
    implicit none

    class(case_abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    real(dp) :: flow_rate_diff

    ! Successor: precursor plane at the inlet, convective outflow at the
    ! outlet (Incompact3d inflow/outflow). The wall stamps below then win
    ! on the edges where the faces meet.
    if (self%inflow_cfg%read_planes) then
      ! Measured on this u against the plane being imposed, so the flow
      ! rate out matches the flow rate in exactly, whatever the
      ! mass_conserve shift or the time step did to the boundary planes.
      ! An imbalance is a uniform divergence the Neumann-x projection
      ! cannot remove.
      flow_rate_diff = self%outflow_flow_rate_correction( &
                       u, self%out_vel, self%flow_rate_y_weights)
      call self%solver%backend%field_set_face_from_field( &
        u, self%bc_start_u_x, self%out_vel, X_FACE, &
        flow_rate_diff=flow_rate_diff)
      ! The flow-rate correction balances the streamwise flux, so only u
      ! takes it
      call self%solver%backend%field_set_face_from_field( &
        v, self%bc_start_v_x, self%out_vel, X_FACE)
      call self%solver%backend%field_set_face_from_field( &
        w, self%bc_start_w_x, self%out_vel, X_FACE)
    end if

    ! No-slip floor, free-slip lid (Incompact3d ncly1=2, nclyn=1). The
    ! tangential components are pinned only at the floor; the lid is left to
    ! the even-symmetry closure. The wall stress does not come from this
    ! resolved gradient -- the wall model supplies it to the SGS stress.
    call self%solver%backend%field_set_face(u, 0._dp, 0._dp, Y_FACE, &
                                            bc_start=BC_DIRICHLET, &
                                            bc_end=BC_NEUMANN)
    call self%solver%backend%field_set_face(w, 0._dp, 0._dp, Y_FACE, &
                                            bc_start=BC_DIRICHLET, &
                                            bc_end=BC_NEUMANN)
    ! Impermeability at both faces. The Neumann pressure BC gives a zero
    ! projection correction on the boundary planes, so without this stamp the
    ! wall-normal velocity there feeds an unstable uniform-divergence mode.
    call self%solver%backend%field_set_face(v, 0._dp, 0._dp, Y_FACE)

  end subroutine apply_BC_abl

  subroutine postprocess_abl(self, iter, t)
    implicit none

    class(case_abl_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    call log_info('time = '//str(t)//' iteration = '//str(iter))

    call self%monitoring%write_step( &
      self%solver, t, self%solver%u, self%solver%v, self%solver%w)
    call self%diagnostics%sample( &
      iter, t, self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_abl

  subroutine post_step_abl(self, iter, t)
    implicit none

    class(case_abl_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    call self%inflow_writer%sample( &
      iter, t, self%solver%u, self%solver%v, self%solver%w)

    ! Every iturboutput steps, independent of n_output, as in Incompact3d
    if (associated(self%turbine) .and. self%turbine_cfg%iturboutput > 0) then
      if (mod(iter, self%turbine_cfg%iturboutput) == 0) &
        call self%turbine%write_output(iter, self%solver%mesh%par%is_root())
    end if

  end subroutine post_step_abl

  subroutine finalise_abl(self)
    implicit none

    class(case_abl_t) :: self

    if (.not. associated(self%turbine)) return
    call self%turbine%finalise()
    call self%io_mgr%unregister_checkpoint_state()
    deallocate (self%turbine)
  end subroutine finalise_abl

end module m_case_abl
