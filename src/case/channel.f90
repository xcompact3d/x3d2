module m_case_channel
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, DIR_C, VERT, CELL
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_channel_t
  contains
    procedure :: boundary_conditions => boundary_conditions_channel
    procedure :: initial_conditions => initial_conditions_channel
    procedure :: forcings => forcings_channel
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

    real(dp) :: can, ub, coeff, dy, L_y
    integer :: dims(3)

    ub = 0._dp
    dy = self%solver%mesh%geo%d(2)
    L_y = self%solver%mesh%geo%L(2)
    dims = self%solver%mesh%get_global_dims(CELL)
    coeff = dy/(L_y*dims(1)*dims(3))
    ub = ub*coeff
    can = ub - 2._dp/3._dp

    call self%solver%backend%field_shift(self%solver%u, can)

  end subroutine boundary_conditions_channel

  subroutine initial_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    class(field_t), pointer :: u_init

    integer :: i, j, k, dims(3)
    real(dp) :: xloc(3), y

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          y = xloc(2) - self%solver%mesh%geo%L(2)/2

          u_init%data(i, j, k) = 1._dp - y*y
        end do
      end do
    end do

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)

    call self%solver%host_allocator%release_block(u_init)

    call self%solver%v%fill(0._dp)
    call self%solver%w%fill(0._dp)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_channel

  subroutine postprocess_channel(self, i, t)
    implicit none

    class(case_channel_t) :: self
    integer, intent(in) :: i
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) print *, 'time =', t, 'iteration =', i
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_channel

  subroutine forcings_channel(self, du, dv, dw)
    implicit none

    class(case_channel_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

  end subroutine forcings_channel

end module m_case_channel
