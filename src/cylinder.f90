module m_case_cylinder
  !! An example case set up to run a cylinder flow.
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, DIR_C, VERT
  use m_config, only: cylinder_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_cylinder_t
    type(cylinder_config_t) :: cylinder_cfg
  contains
    procedure :: boundary_conditions => boundary_conditions_cylinder
    procedure :: initial_conditions => initial_conditions_cylinder
    procedure :: forcings => forcings_cylinder
    procedure :: pre_correction => pre_correction_cylinder
    procedure :: postprocess => postprocess_cylinder
  end type case_cylinder_t

  interface case_cylinder_t
    module procedure case_cylinder_init
  end interface case_cylinder_t

contains

  function case_cylinder_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_cylinder_t) :: flow_case

    call flow_case%cylinder_cfg%read(nml_file=get_argument(1))

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_cylinder_init

  subroutine boundary_conditions_cylinder(self)
    implicit none

    class(case_cylinder_t) :: self

  end subroutine boundary_conditions_cylinder

  subroutine initial_conditions_cylinder(self)
    implicit none

    class(case_cylinder_t) :: self

    real(dp) :: init_noise(3)
    integer :: dims(3)
    class(field_t), pointer :: u_init, v_init, w_init

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

    ! Initial value is 0.
    call u_init%fill(0._dp)
    call v_init%fill(0._dp)
    call w_init%fill(0._dp)

    ! Set initial noise in [0,1]
    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))

    ! Offset and rescale the noise
    init_noise = self%cylinder_cfg%init_noise
    u_init%data(:, :, :) = 1._dp &
                           + (u_init%data(:, :, :) - 0.5_dp)*init_noise(1)
    v_init%data(:, :, :) = (v_init%data(:, :, :) - 0.5_dp)*init_noise(2)
    w_init%data(:, :, :) = (w_init%data(:, :, :) - 0.5_dp)*init_noise(3)

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_cylinder

  subroutine forcings_cylinder(self, du, dv, dw, iter)
    implicit none

    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

  end subroutine forcings_cylinder

  subroutine pre_correction_cylinder(self, u, v, w)
    implicit none

    class(case_cylinder_t) :: self
    class(field_t), intent(inout) :: u, v, w

  end subroutine pre_correction_cylinder

  subroutine postprocess_cylinder(self, iter, t)
    implicit none

    class(case_cylinder_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_cylinder

end module m_case_cylinder
