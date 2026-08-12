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
    class(field_t), pointer :: sponge_sigma => null() !! outflow relaxation rate
  contains
    procedure :: define_BC => define_BC_wind_turbine
    procedure :: initial_conditions => initial_conditions_wind_turbine
    procedure :: forcings => forcings_wind_turbine
    procedure :: apply_BC => apply_BC_wind_turbine
    procedure :: postprocess => postprocess_wind_turbine
    procedure :: finalise_case_specific => finalise_wind_turbine
    procedure :: compute_outflow_params
    procedure :: init_sponge
    procedure :: apply_sponge
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

    ! Seed the RNG deterministically, offset by rank, so inlet-noise runs are
    ! reproducible for a fixed decomposition while ranks stay decorrelated.
    call seed_rng(flow_case%solver%mesh%par%nrank)

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

    ! Optional outflow sponge/relaxation zone (off by default).
    call flow_case%init_sponge()
  end function case_wind_turbine_init

  subroutine seed_rng(nrank)
    !! Deterministic, rank-offset seeding of the intrinsic RNG.
    integer, intent(in) :: nrank
    ! Odd multiplier that spreads the seed array over distinct values.
    integer, parameter :: seed_stride = 37
    integer :: seed_size, s
    integer, allocatable :: seed(:)

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = [(seed_stride*s + 1, s=1, seed_size)] + nrank
    call random_seed(put=seed)
  end subroutine seed_rng

  ! Initial Conditions: uniform inflow with localised noise
  subroutine initial_conditions_wind_turbine(self)
    implicit none
    class(case_wind_turbine_t) :: self

    class(field_t), pointer :: u_init, v_init, w_init
    integer :: i, j, k, dims(3)
    real(dp) :: noise(3)

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
          u_init%data(i, j, k) = self%wt_cfg%bc_start_u &
                                 + noise(1)*(2*u_init%data(i, j, k) - 1._dp)
          v_init%data(i, j, k) = self%wt_cfg%bc_start_v &
                                 + noise(2)*(2*v_init%data(i, j, k) - 1._dp)
          w_init%data(i, j, k) = self%wt_cfg%bc_start_w &
                                 + noise(3)*(2*w_init%data(i, j, k) - 1._dp)
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

  ! Convective outflow velocity (identical treatment to Incompact3d's
  ! outflow subroutine):  cx = 0.5*(uxmax + uxmin)*gdt/dx, where uxmax and
  ! uxmin are the global max/min of the streamwise velocity on the x = nx-1
  ! plane of the post-integration field.
  subroutine compute_outflow_params(self, out_vel)
    implicit none
    class(case_wind_turbine_t) :: self
    real(dp), intent(out) :: out_vel

    integer :: dims(3), nx, ierr
    real(dp) :: uxmax, uxmin, sum_discard
    real(dp) :: dx, gdt

    dims = self%solver%mesh%get_dims(VERT)
    nx = dims(1)
    dx = self%solver%mesh%geo%d(1)
    gdt = self%solver%time_integrator%gdt

    call self%solver%backend%slice_max_sum( &
      uxmax, sum_discard, self%solver%u, nx - 1, min_val=uxmin)

    call MPI_Allreduce(MPI_IN_PLACE, uxmax, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, uxmin, 1, MPI_X3D2_DP, &
                       MPI_MIN, MPI_COMM_WORLD, ierr)

    out_vel = 0.5_dp*(uxmax + uxmin)*gdt/dx

  end subroutine compute_outflow_params

  subroutine define_BC_wind_turbine(self)
    !! Boundary Conditions hook (called per substep before transeq).
    !! Builds the inlet Dirichlet profile only; outflow parameters are
    !! computed in apply_BC from the post-integration velocity.
    !!   u  = bc_start_u + inlet_noise(1) * um * (2r - 1)
    !!   v  = bc_start_v + inlet_noise(2) * um * (2r - 1)
    !!   w  = bc_start_w + inlet_noise(3) * um * (2r - 1)
    implicit none
    class(case_wind_turbine_t) :: self

    class(field_t), pointer :: hu, hv, hw
    integer :: j, k, dims(3)
    real(dp) :: noise(3), um
    logical :: first_call

    dims = self%solver%mesh%get_dims(VERT)
    noise = self%wt_cfg%inlet_noise
    um = abs(self%wt_cfg%bc_start_u)

    ! Allocate persistent device BC fields on first call
    first_call = .not. associated(self%bc_start_u_x)
    if (first_call) then
      self%bc_start_u_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_v_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_w_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      call self%bc_start_u_x%set_data_loc(VERT)
      call self%bc_start_v_x%set_data_loc(VERT)
      call self%bc_start_w_x%set_data_loc(VERT)
    end if

    ! A noise-free inlet is constant in time: build it once and reuse. Only a
    ! noisy inlet needs refreshing every substep.
    if (.not. first_call .and. all(noise == 0._dp)) return

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
    !! enforce the inflow Dirichlet and convective outflow BCs, then apply
    !! an exact uniform mass correction to the outlet, reproducing the
    !! Incompact3d treatment (outflow + pre_correc).
    !!
    !! Steps (mirroring legacy exactly):
    !!  1. cx = 0.5*(uxmax + uxmin)*gdt/dx from the post-integration field.
    !!  2. Convective outflow at x=nx and Dirichlet inflow at x=1 for u,v,w
    !!     (no correction yet):  u(nx) = u(nx) - cx*(u(nx) - u(nx-1)).
    !!  3. Measure the inlet mean (from the prescribed Dirichlet plane) and
    !!     the outlet mean AFTER the convective update, then shift the outlet
    !!     plane of u uniformly so mean(u_outlet) = mean(u_inlet) exactly
    !!     (legacy: bxxn = bxxn - ut + ut1). Applied to u only.
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: u, v, w

    integer :: dims(3), nx, ierr
    real(dp) :: out_vel, ny_nz, udif, max_discard
    real(dp) :: fl_sums(2)
    integer :: gdims(3)

    dims = self%solver%mesh%get_dims(VERT)
    nx = dims(1)
    ! Use global y-z plane size so the mass-flux average is correct
    ! regardless of MPI decomposition in y or z.
    gdims = self%solver%mesh%get_global_dims(VERT)
    ny_nz = real(gdims(2)*gdims(3), dp)

    ! Convective speed from the already-integrated velocity field.
    call self%compute_outflow_params(out_vel)

    ! Dirichlet inflow + convective outflow (no mass correction yet).
    call self%solver%backend%field_set_face_from_field( &
      u, self%bc_start_u_x, out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
    call self%solver%backend%field_set_face_from_field( &
      v, self%bc_start_v_x, out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
    call self%solver%backend%field_set_face_from_field( &
      w, self%bc_start_w_x, out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

    ! Exact mass conservation (legacy bxxn = bxxn - ut + ut1):
    ! ut1 = mean of the prescribed inlet plane (bc_start_u_x at x=1),
    ! ut  = mean of the outlet plane AFTER the convective update (u at x=nx).
    call self%solver%backend%slice_max_sum( &
      max_discard, fl_sums(1), self%bc_start_u_x, 1)
    call self%solver%backend%slice_max_sum( &
      max_discard, fl_sums(2), u, nx)
    call MPI_Allreduce(MPI_IN_PLACE, fl_sums, 2, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    udif = (fl_sums(1) - fl_sums(2))/ny_nz

    call self%solver%backend%field_add_const_x_face(u, udif, at_end=.true.)
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

    ! Outflow sponge: relax (u,v,w) toward the freestream near the outlet to
    ! damp the turbine wake before it reaches the convective outflow boundary
    ! and reflects. Adds a source term -sigma(x)*(q - q_ref) to each momentum
    ! RHS; the subsequent pressure correction restores divergence-free flow.
    if (self%wt_cfg%sponge_on) then
      call self%apply_sponge(du, self%solver%u, self%wt_cfg%bc_start_u)
      call self%apply_sponge(dv, self%solver%v, self%wt_cfg%bc_start_v)
      call self%apply_sponge(dw, self%solver%w, self%wt_cfg%bc_start_w)
    end if
  end subroutine forcings_wind_turbine

  ! Build the outflow-sponge relaxation-rate field sigma(x). The rate is zero
  ! upstream of sponge_start and ramps quadratically to sponge_strength at the
  ! outlet (x = L_x), giving a smooth (C1) onset that avoids reflections.
  subroutine init_sponge(self)
    implicit none
    class(case_wind_turbine_t) :: self

    class(field_t), pointer :: sigma_host
    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), x0, xL, s, sigma_max

    if (.not. self%wt_cfg%sponge_on) return

    x0 = self%wt_cfg%sponge_start
    xL = self%solver%mesh%geo%L(1)
    sigma_max = self%wt_cfg%sponge_strength

    if (xL <= x0) error stop 'wind_turbine sponge: sponge_start >= L_x.'

    dims = self%solver%mesh%get_dims(VERT)
    sigma_host => self%solver%host_allocator%get_block(DIR_C)
    sigma_host%data(:, :, :) = 0._dp

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k)
          if (coords(1) > x0) then
            s = (coords(1) - x0)/(xL - x0)   ! normalised 0..1
            sigma_host%data(i, j, k) = sigma_max*s*s
          end if
        end do
      end do
    end do

    self%sponge_sigma => self%solver%backend%allocator%get_block(DIR_X, VERT)
    call self%sponge_sigma%set_data_loc(VERT)
    call self%solver%backend%set_field_data(self%sponge_sigma, sigma_host%data)

    call self%solver%host_allocator%release_block(sigma_host)

    if (self%solver%mesh%par%is_root()) then
      print *, 'Outflow sponge enabled: start =', x0, &
        ' strength =', sigma_max
    end if
  end subroutine init_sponge

  ! Apply one component of the sponge relaxation to a momentum RHS:
  !   dq = dq - sigma(x)*(q - q_ref)
  subroutine apply_sponge(self, dq, q, q_ref)
    implicit none
    class(case_wind_turbine_t) :: self
    class(field_t), intent(inout) :: dq
    class(field_t), intent(in) :: q
    real(dp), intent(in) :: q_ref

    class(field_t), pointer :: tmp

    tmp => self%solver%backend%allocator%get_block(DIR_X, VERT)
    call self%solver%backend%veccopy(tmp, q)        ! tmp = q
    if (q_ref /= 0._dp) call self%solver%backend%field_shift(tmp, -q_ref)
    call self%solver%backend%vecmult(tmp, self%sponge_sigma) ! tmp = sigma*(q-q_ref)
    call self%solver%backend%vecadd(-1._dp, tmp, 1._dp, dq)  ! dq -= sigma*(q-q_ref)
    call self%solver%backend%allocator%release_block(tmp)
  end subroutine apply_sponge

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

    if (associated(self%sponge_sigma)) then
      call self%solver%backend%allocator%release_block(self%sponge_sigma)
      self%sponge_sigma => null()
    end if

    call self%turbine%finalise()
    call self%io_mgr%unregister_checkpoint_state()
    deallocate (self%turbine)
  end subroutine finalise_wind_turbine

end module m_case_wind_turbine
