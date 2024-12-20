module m_case_tgv
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, VERT, DIR_C
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

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_tgv_init

  subroutine initial_conditions_tgv(self)
    implicit none

    class(case_tgv_t) :: self

    class(field_t), pointer :: u_init, v_init, w_init

    integer :: i, j, k, dims(3)
    real(dp) :: xloc(3), x, y, z

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          x = xloc(1)
          y = xloc(2)
          z = xloc(3)

          u_init%data(i, j, k) = sin(x)*cos(y)*cos(z)
          v_init%data(i, j, k) = -cos(x)*sin(y)*cos(z)
          w_init%data(i, j, k) = 0
        end do
      end do
    end do

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

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

    if (self%solver%mesh%par%is_root()) print *, 'time =', t
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_tgv

end module m_case_tgv
