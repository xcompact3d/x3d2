module m_case_tgv
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
  contains
    procedure :: boundary_conditions => boundary_conditions_tgv
    procedure :: initial_conditions => initial_conditions_tgv
    procedure :: forcings => forcings_tgv
    procedure :: pre_correction => pre_correction_tgv
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

    call self%set_init(self%solver%u, u_func)
    call self%set_init(self%solver%v, v_func)
    call self%solver%w%fill(0._dp)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_tgv

  pure function u_func(coords) result(r)
    implicit none

    real(dp), intent(in) :: coords(3)
    real(dp) :: r

    r = sin(coords(1))*cos(coords(2))*cos(coords(3))
  end function u_func

  pure function v_func(coords) result(r)
    implicit none

    real(dp), intent(in) :: coords(3)
    real(dp) :: r

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
