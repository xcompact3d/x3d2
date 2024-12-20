module m_case_channel
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_channel_t
  contains
    procedure :: boundary_conditions => boundary_conditions_channel
    procedure :: initial_conditions => initial_conditions_channel
    procedure :: post_transeq => post_transeq_channel
    procedure :: postprocess => postprocess_channel
  end type case_channel_t

  interface case_channel_t
    module procedure case_channel_init
  end interface case_channel_t

contains

  function case_channel_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_channel_t) :: flow_case

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_channel_init

  subroutine boundary_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

  end subroutine boundary_conditions_channel

  subroutine initial_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    call self%solver%u%fill(1._dp)
    call self%solver%v%fill(0._dp)
    call self%solver%w%fill(0._dp)

  end subroutine initial_conditions_channel

  subroutine postprocess_channel(self, t)
    implicit none

    class(case_channel_t) :: self
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) print *, 'time =', t
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_channel

  subroutine post_transeq_channel(self, du, dv, dw)
    implicit none

    class(case_channel_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

  end subroutine post_transeq_channel

end module m_case_channel
