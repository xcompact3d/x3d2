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
    real(dp) :: out_vel = 0._dp
    real(dp) :: flow_rate_diff = 0._dp
    ! The persistent device BC fields (DIR_X, VERT) for the inlet
    ! plane (bc_start_u/v/w_x) live on base_case_t. They are allocated
    ! on the first call to define_BC_cylinder, then refilled every
    ! substep (cheap upload, no allocator churn). Released only at
    ! program end.
  contains
    procedure :: define_BC => define_BC_cylinder
    procedure :: initial_conditions => initial_conditions_cylinder
    procedure :: forcings => forcings_cylinder
    procedure :: apply_BC => apply_BC_cylinder
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
  ! Compute outflow velocity number and inlet/outlet flow-rate imbalance.
  ! Three slice reductions on the GPU, two batched MPI_Allreduces (MAX
  ! scalar, SUM 2-vector). No D2H copy of the full field.
  ! ==========================================================================
  subroutine compute_outflow_params(self, out_vel, flow_rate_diff)
    implicit none
    class(case_cylinder_t) :: self
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
    ! NOTE: preserved from original -- uses local ny*nz, not global. If the
    ! y-z plane is decomposed, this is not the true per-cell flow rate.
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
    flow_rate_in = fl_sums(1)/ny_nz
    flow_rate_out = fl_sums(2)/ny_nz

    out_vel = uxmax*gdt/dx
    flow_rate_diff = flow_rate_in - flow_rate_out
  end subroutine compute_outflow_params

  ! ==========================================================================
  ! Boundary Conditions hook (called per substep before transeq).
  ! Sets up / refreshes the inflow profile fields. No writes to
  ! solver%u/v/w here.
  !
  ! Mirrors initial_conditions structurally, but at a single x-plane:
  !   u  = 1 + noise(1) * um * (2r - 1)
  !   v  = 0 + noise(2) * um * (2r - 1)
  !   w  = 0 + noise(3) * um * (2r - 1)
  ! ==========================================================================
  subroutine define_BC_cylinder(self)
    implicit none
    class(case_cylinder_t) :: self

    class(field_t), pointer :: hu, hv, hw
    integer :: j, k, dims(3)
    real(dp) :: noise(3), um, half_L

    dims = self%solver%mesh%get_dims(VERT)
    noise = self%cylinder_cfg%inlet_noise
    half_L = self%solver%mesh%geo%L(1)/2._dp
    um = exp(-0.2_dp*half_L*half_L)

    ! Sample outflow params from solver%u once per substep, stored on self
    ! so apply_BC (later in the substep) and postprocess can both read them.
    ! Note: these are sampled pre-step here, whereas the old layout sampled
    ! them post-step inside apply_BC; the resulting BC parameters drift by
    ! one substep relative to the stamp.

    ! TODO: remove cache data when the PR#300 merges in or add a commit
    !      to fix it here too accordingly
    call self%compute_outflow_params(self%out_vel, self%flow_rate_diff)

    ! Allocate persistent device BC fields on first call.
    if (.not. associated(self%bc_start_u_x)) then
      self%bc_start_u_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_v_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      self%bc_start_w_x => self%solver%backend%allocator%get_block(DIR_X, VERT)
      call self%bc_start_u_x%set_data_loc(VERT)
      call self%bc_start_v_x%set_data_loc(VERT)
      call self%bc_start_w_x%set_data_loc(VERT)
    end if

    ! Build the inflow profile in DIR_C host buffers, then upload.
    hu => self%solver%host_allocator%get_block(DIR_C)
    hv => self%solver%host_allocator%get_block(DIR_C)
    hw => self%solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        hu%data(1, j, k) = 1._dp + noise(1)*um*(2._dp*hu%data(1, j, k) - 1._dp)
        hv%data(1, j, k) = noise(2)*um*(2._dp*hv%data(1, j, k) - 1._dp)
        hw%data(1, j, k) = noise(3)*um*(2._dp*hw%data(1, j, k) - 1._dp)
      end do
    end do

    call self%solver%backend%set_field_data(self%bc_start_u_x, hu%data)
    call self%solver%backend%set_field_data(self%bc_start_v_x, hv%data)
    call self%solver%backend%set_field_data(self%bc_start_w_x, hw%data)

    call self%solver%host_allocator%release_block(hu)
    call self%solver%host_allocator%release_block(hv)
    call self%solver%host_allocator%release_block(hw)
  end subroutine define_BC_cylinder

  ! ==========================================================================
  ! Pre-correction (called per substep after the integrator step):
  ! enforce inflow Dirichlet from bc_start_u_x/v/w on the inlet plane and the
  ! convective outflow update on the right face.
  ! ==========================================================================
  subroutine apply_BC_cylinder(self, u, v, w)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: u, v, w

    call self%solver%backend%field_set_face_from_field( &
      u, self%bc_start_u_x, self%out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff)
    call self%solver%backend%field_set_face_from_field( &
      v, self%bc_start_v_x, self%out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff)
    call self%solver%backend%field_set_face_from_field( &
      w, self%bc_start_w_x, self%out_vel, X_FACE, &
      bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
      flow_rate_diff=self%flow_rate_diff)
  end subroutine apply_BC_cylinder

  ! ==========================================================================
  ! Forcings: empty -- cylinder forcing is handled by the solver's IBM
  ! ==========================================================================
  subroutine forcings_cylinder(self, du, dv, dw, iter)
    implicit none
    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter
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
      print '(A, ES12.5, A, ES12.5)', &
        ' out_vel = ', self%out_vel, &
        '  flow_rate_diff = ', self%flow_rate_diff
    end if
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)
  end subroutine postprocess_cylinder

end module m_case_cylinder
