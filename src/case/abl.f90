module m_case_abl
  use m_allocator, only: allocator_t
  use m_abl, only: abl_t
  use m_abl_diagnostics, only: abl_diagnostics_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, MPI_X3D2_DP, VERT, Y_FACE, &
                      BC_DIRICHLET, BC_NEUMANN, BC_PERIODIC
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce

  implicit none

  type, extends(base_case_t) :: case_abl_t
    type(abl_config_t) :: abl_cfg
    type(abl_t) :: abl
    type(abl_diagnostics_t) :: diagnostics
  contains
    procedure :: define_BC => define_BC_abl
    procedure :: initial_conditions => initial_conditions_abl
    procedure :: forcings => forcings_abl
    procedure :: apply_BC => apply_BC_abl
    procedure :: postprocess => postprocess_abl
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

    integer :: bcs(3, 2)

    ! apply_BC_abl and the wall model hard-code these boundaries, so the
    ! mesh must agree with them for the pressure solve and precorrect_walls.
    bcs = mesh%grid%BCs_global
    if (any(bcs(1, :) /= BC_PERIODIC) .or. any(bcs(3, :) /= BC_PERIODIC) &
        .or. bcs(2, 1) /= BC_DIRICHLET .or. bcs(2, 2) /= BC_NEUMANN) &
      error stop 'ABL requires periodic BC_x and BC_z, and &
                 &BC_y = ''dirichlet'', ''neumann''.'

    call flow_case%abl_cfg%read(nml_file=get_argument(1))
    flow_case%abl = abl_t(backend, mesh, host_allocator, flow_case%abl_cfg)

    call flow_case%case_init(backend, mesh, host_allocator)
    flow_case%solver%keep_wall_correction = .true.
    call flow_case%abl%configure_wall_boundary_correction(flow_case%solver%les)
    flow_case%diagnostics = abl_diagnostics_t(backend, mesh, &
                                              flow_case%abl_cfg)

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

  end subroutine define_BC_abl

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

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%monitoring%write_step( &
      self%solver, t, self%solver%u, self%solver%v, self%solver%w)
    call self%diagnostics%sample( &
      iter, t, self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_abl

end module m_case_abl
