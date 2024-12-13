module m_case_generic
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_generic_t
  contains
    procedure :: boundary_conditions => boundary_conditions_generic
    procedure :: initial_conditions => initial_conditions_generic
    procedure :: post_transeq => post_transeq_generic
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

    flow_case%solver = init(backend, mesh, host_allocator)
  end function case_generic_init

  subroutine boundary_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

  end subroutine boundary_conditions_generic

  subroutine initial_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

  end subroutine initial_conditions_generic

  subroutine postprocess_generic(self, t)
    implicit none

    class(case_generic_t) :: self
    real(dp), intent(in) :: t

  end subroutine postprocess_generic

  subroutine post_transeq_generic(self, du, dv, dw)
    implicit none

    class(case_generic_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

  end subroutine post_transeq_generic

end module m_case_generic
