module m_case_generic
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, init
  use m_tdsops, only: tdsops_t, dirps_t
  use m_time_integrator, only: time_intg_t
  use m_vector_calculus, only: vector_calculus_t

  implicit none

  type, extends(solver_t) :: case_generic_t
  contains
    procedure :: post_transeq => post_transeq_generic
    procedure :: boundary_conditions => boundary_conditions_generic
    procedure :: postprocess => postprocess_generic
  end type case_generic_t

  interface case_generic_t
    module procedure case_generic_init
  end interface case_generic_t

contains

  function case_generic_init(backend, mesh, host_allocator) result(solver)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_generic_t) :: solver

    solver%solver_t = init(backend, mesh, host_allocator)
  end function case_generic_init

  subroutine boundary_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

    call self%solver_t%boundary_conditions()

  end subroutine boundary_conditions_generic

  subroutine post_transeq_generic(self, du, dv, dw)
    implicit none

    class(case_generic_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

    ! first call the parent class
    call self%solver_t%post_transeq(du, dv, dw)

  end subroutine post_transeq_generic

  subroutine postprocess_generic(self, t)
    implicit none

    class(case_generic_t), intent(in) :: self
    real(dp), intent(in) :: t

    call self%solver_t%postprocess(t)

  end subroutine postprocess_generic

end module m_case_generic
