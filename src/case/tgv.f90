module m_case_tgv
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

  type, extends(solver_t) :: case_tgv_t
  contains
    procedure :: post_transeq => post_transeq_tgv
  end type case_tgv_t

  interface case_tgv_t
    module procedure case_tgv_init
  end interface case_tgv_t

contains


  function case_tgv_init(backend, mesh, host_allocator) result(solver)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_tgv_t) :: solver

    solver%solver_t = init(backend, mesh, host_allocator)
  end function case_tgv_init

  subroutine post_transeq_tgv(self, du, dv, dw)
    implicit none

    class(case_tgv_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

    ! first call the parent class
    call self%solver_t%post_transeq(du, dv, dw)

    print *, 'post_transeq for the tgv case'
  end subroutine post_transeq_tgv

end module m_case_tgv
