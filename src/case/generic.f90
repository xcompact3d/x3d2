module m_case_generic
  !! An example case set up to run and sustain a freestream flow.
  !! This is a good place to start for adding a new flow case.
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
  contains
    procedure :: boundary_conditions => boundary_conditions_generic
    procedure :: initial_conditions => initial_conditions_generic
    procedure :: forcings => forcings_generic
    procedure :: pre_correction => pre_correction_generic
    procedure :: postprocess => postprocess_generic
  end type case_generic_t

  interface case_generic_t
    module procedure case_generic_init
  end interface case_generic_t

contains

  function case_generic_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_generic_t) :: flow_case

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_generic_init

  subroutine boundary_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

  end subroutine boundary_conditions_generic

  subroutine initial_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

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
