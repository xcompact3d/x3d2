module m_case_generic
  !! Generic freestream flow case for general-purpose simulations.
  !!
  !! This module provides a minimal template for setting up custom flow
  !! cases. It implements a simple uniform freestream flow (u=1, v=0, w=0)
  !! with no forcing or boundary corrections.
  !!
  !! **Use Cases:**
  !! - Starting point for implementing new flow cases
  !! - Testing solver functionality with simple initial conditions
  !! - Freestream simulations with immersed boundaries (add IBM via forcings)
  !! - Custom flow setups requiring minimal default behaviour
  !!
  !! **Default Configuration:**
  !! - Initial condition: Uniform flow u=1, v=w=0
  !! - No boundary condition corrections
  !! - No forcing terms
  !! - No pre-correction
  !! - Minimal postprocessing
  !!
  !! **Customisation:**
  !! Users can extend this case or modify the procedures directly to implement
  !! specific flow physics, boundary conditions, or forcing terms.
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_generic_t
    !! Generic case with minimal default behaviour.
  contains
    procedure :: boundary_conditions => boundary_conditions_generic !! No action (use domain BCs)
    procedure :: initial_conditions => initial_conditions_generic   !! Uniform freestream
    procedure :: forcings => forcings_generic                       !! No forcing
    procedure :: pre_correction => pre_correction_generic           !! No correction
    procedure :: postprocess => postprocess_generic                 !! Minimal diagnostics
  end type case_generic_t

  interface case_generic_t
    module procedure case_generic_init
  end interface case_generic_t

contains

  function case_generic_init(backend, mesh, host_allocator) result(flow_case)
    !! Initialise generic flow case.
    implicit none

    class(base_backend_t), target, intent(inout) :: backend         !! Computational backend
    type(mesh_t), target, intent(inout) :: mesh                     !! Mesh with decomposition
    type(allocator_t), target, intent(inout) :: host_allocator      !! Host memory allocator
    type(case_generic_t) :: flow_case                               !! Initialised generic case

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_generic_init

  subroutine boundary_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

  end subroutine boundary_conditions_generic

  subroutine initial_conditions_generic(self)
    !! Set initial velocity field for generic freestream case.
    !!
    !! Initialises a uniform flow field with:
    !! - \( u = 1 \) (streamwise velocity)
    !! - \( v = 0 \) (cross-stream velocity)
    !! - \( w = 0 \) (spanwise velocity)
    !!
    !! All velocity components are located at vertices (VERT).
    !! This simple uniform flow serves as a starting point that users
    !! can modify for their specific applications.
    implicit none

    class(case_generic_t) :: self !! Generic case instance

    call self%solver%u%fill(1._dp)
    call self%solver%v%fill(0._dp)
    call self%solver%w%fill(0._dp)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_generic

  subroutine forcings_generic(self, du, dv, dw, iter)
    implicit none

    class(case_generic_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

  end subroutine forcings_generic

  subroutine pre_correction_generic(self, u, v, w)
    implicit none

    class(case_generic_t) :: self
    class(field_t), intent(inout) :: u, v, w

  end subroutine pre_correction_generic

  subroutine postprocess_generic(self, iter, t)
    implicit none

    class(case_generic_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_generic

end module m_case_generic
