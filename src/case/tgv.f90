module m_case_tgv
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_tgv_t
  contains
    procedure :: boundary_conditions => boundary_conditions_tgv
    procedure :: initial_conditions => initial_conditions_tgv
    procedure :: post_transeq => post_transeq_tgv
    procedure :: postprocess => postprocess_tgv
  end type case_tgv_t

  interface case_tgv_t
    module procedure case_tgv_init
  end interface case_tgv_t

contains

  function case_tgv_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_tgv_t) :: flow_case

    flow_case%solver = init(backend, mesh, host_allocator)
  end function case_tgv_init

  subroutine initial_conditions_tgv(self)
    implicit none

    class(case_tgv_t) :: self

    ! set initial conditions for the TGV case
  end subroutine initial_conditions_tgv

  subroutine boundary_conditions_tgv(self)
    implicit none

    class(case_tgv_t) :: self

    ! do nothing for TGV case
  end subroutine boundary_conditions_tgv

  subroutine post_transeq_tgv(self, du, dv, dw)
    implicit none

    class(case_tgv_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

    ! do nothing for TGV case
  end subroutine post_transeq_tgv

  subroutine postprocess_tgv(self, t)
    implicit none

    class(case_tgv_t) :: self
    real(dp), intent(in) :: t

    call self%solver%print_enstrophy(self%solver%u, self%solver%v, &
                                     self%solver%w)
    call self%solver%print_div_max_mean(self%solver%u, self%solver%v, &
                                        self%solver%w)

  end subroutine postprocess_tgv

end module m_case_tgv
