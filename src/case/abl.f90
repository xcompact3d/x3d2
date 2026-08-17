module m_case_abl
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
  use m_abl, only: abl_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, VERT
  use m_config, only: abl_config_t, solver_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_abl_t
    type(abl_config_t) :: abl_cfg
    type(abl_t) :: abl
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

  end function case_abl_init

  subroutine initial_conditions_abl(self)
    implicit none

    class(case_abl_t) :: self

    call self%abl%initialise(self%solver%u, self%solver%v, self%solver%w)

  end subroutine initial_conditions_abl

  subroutine define_BC_abl(self)
    implicit none

    class(case_abl_t) :: self

    ! Free-slip y-walls are handled by the Poisson solver via the mesh BCs.
    ! The wall model (and mass-conserve correction) wire in here in later
    ! issue #317 commits.

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

    ! Free-slip; nothing to stamp onto the velocity fields.

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

  end subroutine postprocess_abl

end module m_case_abl
