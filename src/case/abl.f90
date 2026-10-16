module m_case_abl
  use m_allocator, only: allocator_t
  use m_abl, only: abl_t
  use m_abl_diagnostics, only: abl_diagnostics_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, MPI_X3D2_DP, DIR_C, DIR_X, VERT, &
                      X_FACE, Y_FACE, BC_DIRICHLET, BC_NEUMANN, BC_PERIODIC
  use m_config, only: abl_config_t, inflow_config_t
  use m_field, only: field_t
  use m_inflow_planes, only: inflow_plane_writer_t, inflow_plane_reader_t
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce
  use m_logger, only: log_info, abort_run, str

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
    real(dp) :: flow_rate_diff = 0._dp
  contains
    procedure :: define_BC => define_BC_abl
    procedure :: initial_conditions => initial_conditions_abl
    procedure :: forcings => forcings_abl
    procedure :: apply_BC => apply_BC_abl
    procedure :: postprocess => postprocess_abl
    procedure :: post_step => post_step_abl
    procedure :: update_inflow_plane
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

    integer :: bcs(3, 2), bc_x

    call flow_case%abl_cfg%read(nml_file=get_argument(1))
    call flow_case%inflow_cfg%read(nml_file=get_argument(1))

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
    flow_case%abl = abl_t(backend, mesh, host_allocator, flow_case%abl_cfg)

    call flow_case%case_init(backend, mesh, host_allocator)
    flow_case%solver%keep_wall_correction = .true.
    call flow_case%abl%configure_wall_boundary_correction(flow_case%solver%les)
    flow_case%diagnostics = abl_diagnostics_t(backend, mesh, &
                                              flow_case%abl_cfg)
    flow_case%inflow_writer = inflow_plane_writer_t( &
                              backend, mesh, flow_case%inflow_cfg, &
                              flow_case%solver%dt)
    flow_case%inflow_reader = inflow_plane_reader_t(mesh, flow_case%inflow_cfg)

  end function case_abl_init

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

    if (self%inflow_cfg%read_planes) then
      call self%compute_outflow_params(self%out_vel, self%flow_rate_diff)
      call self%update_inflow_plane()
    end if

  end subroutine define_BC_abl

  subroutine update_inflow_plane(self)
    !! Fill bc_start_*_x with the precursor plane for the iteration being
    !! computed. The plane is held for every substep of the iteration, as
    !! Incompact3d does, so it is uploaded once per iteration.
    implicit none

    class(case_abl_t) :: self

    class(field_t), pointer :: hu, hv, hw
    real(dp), allocatable :: pu(:, :), pv(:, :), pw(:, :)
    integer :: iter, dims(3), j, k

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

    dims = self%solver%mesh%get_dims(VERT)
    hu => self%solver%host_allocator%get_block(DIR_C)
    hv => self%solver%host_allocator%get_block(DIR_C)
    hw => self%solver%host_allocator%get_block(DIR_C)

    ! Only the x = 1 plane of the BC fields is read by the face stamp
    allocate (pu(dims(2), dims(3)), pv(dims(2), dims(3)), pw(dims(2), dims(3)))
    call self%inflow_reader%planes_at(iter*self%solver%dt, pu, pv, pw)
    do k = 1, dims(3)
      do j = 1, dims(2)
        hu%data(1, j, k) = pu(j, k)
        hv%data(1, j, k) = pv(j, k)
        hw%data(1, j, k) = pw(j, k)
      end do
    end do

    call self%solver%backend%set_field_data(self%bc_start_u_x, hu%data)
    call self%solver%backend%set_field_data(self%bc_start_v_x, hv%data)
    call self%solver%backend%set_field_data(self%bc_start_w_x, hw%data)

    call self%solver%host_allocator%release_block(hu)
    call self%solver%host_allocator%release_block(hv)
    call self%solver%host_allocator%release_block(hw)
  end subroutine update_inflow_plane

  subroutine forcings_abl(self, du, dv, dw, iter)
    implicit none

    class(case_abl_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

    call self%abl%apply_forcing(du, dv, dw, &
                                self%solver%u, self%solver%v, self%solver%w)

  end subroutine forcings_abl

  subroutine apply_BC_abl(self, u, v, w)
    implicit none

    class(case_abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    ! Successor: precursor plane at the inlet, convective outflow at the
    ! outlet (Incompact3d inflow/outflow). The wall stamps below then win
    ! on the edges where the faces meet.
    if (self%inflow_cfg%read_planes) then
      call self%solver%backend%field_set_face_from_field( &
        u, self%bc_start_u_x, self%out_vel, X_FACE, &
        flow_rate_diff=self%flow_rate_diff)
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

  end subroutine post_step_abl

end module m_case_abl
