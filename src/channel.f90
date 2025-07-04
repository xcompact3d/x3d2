module m_case_channel
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, get_argument, DIR_C, VERT, CELL, Y_FACE
  use m_config, only: channel_config_t
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_channel_t
    type(channel_config_t) :: channel_cfg
  contains
    procedure :: boundary_conditions => boundary_conditions_channel
    procedure :: initial_conditions => initial_conditions_channel
    procedure :: forcings => forcings_channel
    procedure :: pre_correction => pre_correction_channel
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

    call flow_case%channel_cfg%read(nml_file=get_argument(1))

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_channel_init

  subroutine boundary_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    real(dp) :: can, ub
    integer :: ierr

    ub = self%solver%backend%field_volume_integral(self%solver%u)

    ub = ub/(product(self%solver%mesh%get_global_dims(CELL)))
    call MPI_Allreduce(MPI_IN_PLACE, ub, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

    can = 2._dp/3._dp - ub

    call self%solver%backend%field_shift(self%solver%u, can)

  end subroutine boundary_conditions_channel

  subroutine initial_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    class(field_t), pointer :: u_init, v_init, w_init

    integer :: i, j, k, dims(3)
    real(dp) :: xloc(3), y, noise, um

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(v_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(w_init%data(1:dims(1), 1:dims(2), 1:dims(3)))

    noise = self%channel_cfg%noise
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          y = xloc(2) - self%solver%mesh%geo%L(2)/2._dp
          um = exp(-0.2_dp*y*y)

          u_init%data(i, j, k) = 1._dp - y*y &
                                 + noise*um*(2*u_init%data(i, j, k) - 1._dp)
          v_init%data(i, j, k) = noise*um*(2*v_init%data(i, j, k) - 1._dp)
          w_init%data(i, j, k) = noise*um*(2*w_init%data(i, j, k) - 1._dp)
        end do
      end do
    end do

    u_init%data(:, 1, :) = 0
    v_init%data(:, 1, :) = 0
    w_init%data(:, 1, :) = 0
    u_init%data(:, dims(2), :) = 0
    v_init%data(:, dims(2), :) = 0
    w_init%data(:, dims(2), :) = 0

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_channel

  subroutine forcings_channel(self, du, dv, dw, iter)
    implicit none

    class(case_channel_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: iter

    real(dp) :: rot

    if (self%channel_cfg%rotation .and. iter < self%channel_cfg%n_rotate) then
      rot = self%channel_cfg%omega_rot
      call self%solver%backend%vecadd(-rot, self%solver%v, 1._dp, du)
      call self%solver%backend%vecadd(rot, self%solver%u, 1._dp, dv)
    end if

  end subroutine forcings_channel

  subroutine pre_correction_channel(self, u, v, w)
    implicit none

    class(case_channel_t) :: self
    class(field_t), intent(inout) :: u, v, w

    call self%solver%backend%field_set_face(u, 0._dp, 0._dp, Y_FACE)
    call self%solver%backend%field_set_face(v, 0._dp, 0._dp, Y_FACE)
    call self%solver%backend%field_set_face(w, 0._dp, 0._dp, Y_FACE)

  end subroutine pre_correction_channel

  subroutine postprocess_channel(self, iter, t)
    implicit none

    class(case_channel_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    if (self%solver%mesh%par%is_root()) then
      print *, 'time =', t, 'iteration =', iter
    end if

    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

  end subroutine postprocess_channel

end module m_case_channel
