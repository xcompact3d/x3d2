module m_turbine_dummy
  !! Dummy (no-op) turbine model.
  !!
  !! Implements turbine_model_t but adds zero forcing. It is the explicit
  !! "no turbine" model selected when iturbine = 0, and also lets the
  !! wind-turbine carrier case (m_case_wind_turbine) build and run before any
  !! real turbine physics exists. The actuator disc model and
  !! actuator line model provide their own turbine_model_t
  !! implementations alongside this one.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_io_session, only: reader_session_t, writer_session_t
  use m_turbine_model, only: turbine_model_t

  implicit none

  type, extends(turbine_model_t) :: turbine_dummy_t
  contains
    procedure :: init => init_dummy
    procedure :: update => update_dummy
    procedure :: project_forces => project_forces_dummy
    procedure :: write_output => write_output_dummy
    procedure :: setup_output => setup_output_dummy
    procedure :: write_checkpoint => write_checkpoint_dummy
    procedure :: read_checkpoint => read_checkpoint_dummy
    procedure :: finalise => finalise_dummy
  end type turbine_dummy_t

contains

  subroutine init_dummy(self, backend, mesh, host_allocator, dt)
    implicit none
    class(turbine_dummy_t), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), target, intent(in) :: mesh
    type(allocator_t), target, intent(in) :: host_allocator
    real(dp), intent(in) :: dt
    ! no state to initialise
  end subroutine init_dummy

  subroutine update_dummy(self, t, dt)
    implicit none
    class(turbine_dummy_t), intent(inout) :: self
    real(dp), intent(in) :: t, dt
    ! no kinematics to advance
  end subroutine update_dummy

  subroutine project_forces_dummy(self, du, dv, dw, u, v, w)
    implicit none
    class(turbine_dummy_t), intent(inout) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    ! no momentum source added
  end subroutine project_forces_dummy

  subroutine write_output_dummy(self, iter, is_root)
    implicit none
    class(turbine_dummy_t), intent(inout) :: self
    integer, intent(in) :: iter
    logical, intent(in) :: is_root
    ! no diagnostics to write
  end subroutine write_output_dummy

  subroutine setup_output_dummy(self, is_root, append)
    class(turbine_dummy_t), intent(inout) :: self
    logical, intent(in) :: is_root, append
  end subroutine setup_output_dummy

  subroutine write_checkpoint_dummy(self, writer)
    class(turbine_dummy_t), intent(inout) :: self
    type(writer_session_t), intent(inout) :: writer
  end subroutine write_checkpoint_dummy

  subroutine read_checkpoint_dummy(self, reader)
    class(turbine_dummy_t), intent(inout) :: self
    type(reader_session_t), intent(inout) :: reader
  end subroutine read_checkpoint_dummy

  subroutine finalise_dummy(self)
    class(turbine_dummy_t), intent(inout) :: self
  end subroutine finalise_dummy

end module m_turbine_dummy
