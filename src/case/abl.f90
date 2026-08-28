module m_case_abl
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
  use m_abl, only: abl_t
  use m_abl_diagnostics, only: abl_diagnostics_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, MPI_X3D2_DP, CELL, VERT, Y_FACE
  use m_config, only: abl_config_t, solver_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

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

    type(solver_config_t) :: solver_cfg

    call flow_case%abl_cfg%read(nml_file=get_argument(1))
    ! dt is needed by the driver but only known via the solver namelist.
    call solver_cfg%read(nml_file=get_argument(1))
    flow_case%abl = abl_t(backend, mesh, host_allocator, &
                          flow_case%abl_cfg, solver_cfg%dt)

    call flow_case%case_init(backend, mesh, host_allocator)
    call flow_case%abl%configure_wall_boundary_correction(flow_case%solver%les)
    flow_case%diagnostics = abl_diagnostics_t( &
                            backend, mesh, host_allocator, flow_case%abl_cfg)

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
    integer :: ierr

    ! Constant-flow-rate correction (Incompact3d forceabl); mirrors the channel
    ! bulk-velocity shift, targeting the log-law flow rate. Free-slip y-walls
    ! are otherwise handled by the Poisson solver; the wall model wires in
    ! here in issue #317 commit 3.
    if (self%abl_cfg%mass_conserve) then
      ly = self%solver%mesh%geo%L(2)
      ub = self%solver%backend%field_volume_integral(self%solver%u)
      ub = ub/product(self%solver%mesh%get_global_dims(CELL))
      call MPI_Allreduce(MPI_IN_PLACE, ub, 1, MPI_X3D2_DP, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
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

    ! Slip walls: u and w are free, but impermeability must be enforced
    ! explicitly. The Neumann pressure BC gives a zero projection correction
    ! on the boundary planes, so without this stamp the wall-normal velocity
    ! there integrates freely and feeds an unstable uniform-divergence mode.
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
