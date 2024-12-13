module m_case_generic
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, DIR_C, VERT
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

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_generic_init

  subroutine boundary_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

  end subroutine boundary_conditions_generic

  subroutine initial_conditions_generic(self)
    implicit none

    class(case_generic_t) :: self

    class(field_t), pointer :: u_init, v_init, w_init

    integer :: i, j, k, dims(3)

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          u_init%data(i, j, k) = 1_dp
          v_init%data(i, j, k) = 0
          w_init%data(i, j, k) = 0
        end do
      end do
    end do

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

  end subroutine initial_conditions_generic

  subroutine postprocess_generic(self, t)
    implicit none

    class(case_generic_t) :: self
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) print *, 'time =', t
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_generic

  subroutine post_transeq_generic(self, du, dv, dw)
    implicit none

    class(case_generic_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

  end subroutine post_transeq_generic

end module m_case_generic
