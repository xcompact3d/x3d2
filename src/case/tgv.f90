module m_case_tgv
  !! Taylor-Green vortex (TGV) case for validation and benchmarking.
  !!
  !! The Taylor-Green vortex is a canonical test case for incompressible
  !! Navier-Stokes solvers. It features an analytically-defined initial
  !! condition that transitions from laminar to turbulent flow, providing
  !! a rigorous test of:
  !!
  !! - Spatial discretisation accuracy
  !! - Time integration stability
  !! - Energy conservation properties
  !! - Transition to turbulence physics
  !!
  !! **Initial Conditions:**
  !!
  !! \[ u = \sin(x) \cos(y) \cos(z) \]
  !! \[ v = -\cos(x) \sin(y) \cos(z) \]
  !! \[ w = 0 \]
  !!
  !! This satisfies incompressibility (\( \nabla \cdot \mathbf{u} = 0 \)) exactly and is periodic
  !! in all three directions.
  !!
  !! **Domain:**
  !!
  !! Typically \( [0, 2\pi]^3 \) with periodic boundary conditions in all directions.
  !!
  !! **Validation Metrics:**
  !!
  !! - Kinetic energy decay rate
  !! - Enstrophy evolution
  !! - Dissipation rate
  !! - Vorticity dynamics
  !!
  !! **Reference:**
  !!
  !! Taylor, G. I., & Green, A. E. (1937). Mechanism of the production of
  !! small eddies from large ones. Proc. R. Soc. Lond. A, 158(895), 499-521.
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, VERT, DIR_C
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_tgv_t
    !! Taylor-Green vortex case (no additional state needed beyond base).
  contains
    procedure :: boundary_conditions => boundary_conditions_tgv !! No action (periodic BCs)
    procedure :: initial_conditions => initial_conditions_tgv   !! Set TGV velocity field
    procedure :: forcings => forcings_tgv                       !! No forcing
    procedure :: pre_correction => pre_correction_tgv           !! No correction
    procedure :: postprocess => postprocess_tgv                 !! Compute diagnostics
  end type case_tgv_t

  interface case_tgv_t
    module procedure case_tgv_init
  end interface case_tgv_t

contains

  function case_tgv_init(backend, mesh, host_allocator) result(flow_case)
    !! Initialise Taylor-Green vortex case.
    implicit none

    class(base_backend_t), target, intent(inout) :: backend         !! Computational backend
    type(mesh_t), target, intent(inout) :: mesh                     !! Mesh with decomposition
    type(allocator_t), target, intent(inout) :: host_allocator      !! Host memory allocator
    type(case_tgv_t) :: flow_case                                   !! Initialised TGV case

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_tgv_init

  subroutine initial_conditions_tgv(self)
    !! Set Taylor-Green vortex initial velocity field.
    !!
    !! Initialises the three velocity components according to the TGV
    !! analytical solution. The field is exactly divergence-free and
    !! periodic, making it ideal for testing solver accuracy.
    implicit none

    class(case_tgv_t) :: self !! TGV case instance

    call self%set_init(self%solver%u, u_func)
    call self%set_init(self%solver%v, v_func)
    call self%solver%w%fill(0._dp)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_tgv

  pure function u_func(coords) result(r)
    !! Compute x-velocity component of TGV at given coordinates.
    !!
    !! \[ u = \sin(x) \cos(y) \cos(z) \]
    implicit none

    real(dp), intent(in) :: coords(3) !! Position [x, y, z]
    real(dp) :: r                     !! Velocity component u

    r = sin(coords(1))*cos(coords(2))*cos(coords(3))
  end function u_func

  pure function v_func(coords) result(r)
    !! Compute y-velocity component of TGV at given coordinates.
    !!
    !! \[ v = -\cos(x) \sin(y) \cos(z) \]
    implicit none

    real(dp), intent(in) :: coords(3) !! Position [x, y, z]
    real(dp) :: r                     !! Velocity component v

    r = -cos(coords(1))*sin(coords(2))*cos(coords(3))
  end function v_func

  subroutine boundary_conditions_tgv(self)
    implicit none

    class(case_tgv_t) :: self

    ! do nothing for TGV case
  end subroutine boundary_conditions_tgv

  subroutine forcings_tgv(self, du, dv, dw, iter)
    implicit none

    class(case_tgv_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

    ! do nothing for TGV case
  end subroutine forcings_tgv

  subroutine pre_correction_tgv(self, u, v, w)
    implicit none

    class(case_tgv_t) :: self
    class(field_t), intent(inout) :: u, v, w

    ! do nothing for TGV case
  end subroutine pre_correction_tgv

  subroutine postprocess_tgv(self, iter, t)
    implicit none

    class(case_tgv_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_tgv

end module m_case_tgv
