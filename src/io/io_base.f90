module m_io_base
!! Base types and constants for session-based I/O architecture
!!
!! This module provides the fundamental building blocks for X3D2's modern I/O system:
!!
!! **Architecture Overview:**
!! - Session-based interface (`io_session_t` in `m_io_session`) is the primary user interface
!! - Factory pattern (`m_io_factory`) handles backend selection (only ADIOS2 currently)
!! - Concrete backends (ADIOS2) extend the base types defined here
!! - Base types provide polymorphic containers and enforce consistent interfaces
!!
!! **Design Philosophy:**
!! - Users interact only with `io_session_t` - no manual file handle management
!! - Backends are selected automatically based on compile-time availability
!! - Base types provide fallback implementations that error out (forcing proper override)
!! - Clean separation: session manages complexity, backends focus on I/O specifics
!!
!! **Usage (via session interface):**
!!   type(io_session_t) :: io_session
!!   call io_session%open("output.bp", MPI_COMM_WORLD, io_mode_write)
!!   call io_session%write_data("timestep", current_step)
!!   call io_session%write_data("velocity", u_field, start_dims, count_dims)
!!   call io_session%close()

  use m_common, only: dp, i8

  implicit none

  private
  public :: io_reader_t, io_writer_t, io_file_t
  public :: io_mode_read, io_mode_write

  integer, parameter :: io_mode_read = 1
  integer, parameter :: io_mode_write = 2

  !> Base file handle for I/O operations
  type :: io_file_t
  contains
    procedure :: close => base_close
    procedure :: begin_step => base_begin_step
    procedure :: end_step => base_end_step
  end type io_file_t

  !> Base I/O reader type for polymorphic usage
  type :: io_reader_t
  contains
    procedure :: init => base_reader_init
    procedure :: open => base_reader_open
    procedure :: finalise => base_reader_finalise
    ! Generic interfaces for session usage
    generic :: read_data => read_data_i8, read_data_integer, read_data_real, &
      read_data_array_3d
    procedure :: read_data_i8
    procedure :: read_data_integer
    procedure :: read_data_real
    procedure :: read_data_array_3d
  end type io_reader_t

  !> Base I/O writer type for polymorphic usage
  type :: io_writer_t
  contains
    procedure :: init => base_writer_init
    procedure :: open => base_writer_open
    procedure :: finalise => base_writer_finalise
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, &
      write_data_array_3d
    procedure :: write_data_i8
    procedure :: write_data_integer
    procedure :: write_data_real
    procedure :: write_data_array_3d
    generic :: write_attribute => write_attribute_string, &
      write_attribute_array_1d_real
    procedure :: write_attribute_string
    procedure :: write_attribute_array_1d_real
  end type io_writer_t

contains

  ! Base implementations (should be overridden by concrete types)
  subroutine base_close(self)
    class(io_file_t), intent(inout) :: self
    error stop "base_close should not be called - &
      & use concrete implementation"
  end subroutine base_close

  subroutine base_begin_step(self)
    class(io_file_t), intent(inout) :: self
    error stop "base_begin_step should not be called - &
      & use concrete implementation"
  end subroutine base_begin_step

  subroutine base_end_step(self)
    class(io_file_t), intent(inout) :: self
    error stop "base_end_step should not be called - &
      & use concrete implementation"
  end subroutine base_end_step

  subroutine base_reader_init(self, comm, name)
    class(io_reader_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    error stop "base_reader_init should not be called - &
      & use concrete implementation"
  end subroutine base_reader_init

  function base_reader_open(self, filename, mode, comm) result(file_handle)
    class(io_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), pointer :: file_handle
    error stop "base_reader_open should not be called - &
      & use concrete implementation"
  end function base_reader_open

  subroutine base_reader_finalise(self)
    class(io_reader_t), intent(inout) :: self
    error stop "base_reader_finalise should not be called - &
      & use concrete implementation"
  end subroutine base_reader_finalise

  subroutine base_writer_init(self, comm, name)
    class(io_writer_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    error stop "base_writer_init should not be called - &
      & use concrete implementation"
  end subroutine base_writer_init

  function base_writer_open(self, filename, mode, comm) result(file_handle)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), pointer :: file_handle
    error stop "base_writer_open should not be called - &
      & use concrete implementation"
  end function base_writer_open

  subroutine base_writer_finalise(self)
    class(io_writer_t), intent(inout) :: self
    error stop "base_writer_finalise should not be called - &
      & use concrete implementation"
  end subroutine base_writer_finalise

  ! Base read implementations
  subroutine read_data_i8(self, variable_name, value, file_handle)
    class(io_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "read_data_i8 should not be called - &
      & use concrete implementation"
  end subroutine read_data_i8

  subroutine read_data_integer(self, variable_name, value, file_handle)
    class(io_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "read_data_integer should not be called - &
      & use concrete implementation"
  end subroutine read_data_integer

  subroutine read_data_real(self, variable_name, value, file_handle)
    class(io_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "read_data_real should not be called - &
      & use concrete implementation"
  end subroutine read_data_real

  subroutine read_data_array_3d( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims &
    )
    class(io_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)
    error stop "read_data_array_3d should not be called - &
      & use concrete implementation"
  end subroutine read_data_array_3d

  ! Base write implementations
  subroutine write_data_i8(self, variable_name, value, file_handle)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "write_data_i8 should not be called - &
      & use concrete implementation"
  end subroutine write_data_i8

  subroutine write_data_integer(self, variable_name, value, file_handle)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "write_data_integer should not be called - &
      & use concrete implementation"
  end subroutine write_data_integer

  subroutine write_data_real(self, variable_name, value, file_handle)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "write_data_real should not be called - &
      & use concrete implementation"
  end subroutine write_data_real

  subroutine write_data_array_3d( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims &
    )
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)
    error stop "write_data_array_3d should not be called - &
      & use concrete implementation"
  end subroutine write_data_array_3d

  subroutine write_attribute_string(self, attribute_name, value, file_handle)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    error stop "write_attribute_string should not be called - &
      & use concrete implementation"
  end subroutine write_attribute_string

  subroutine write_attribute_array_1d_real( &
    self, attribute_name, values, file_handle &
    )
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    real(dp), intent(in) :: values(:)
    class(io_file_t), intent(inout) :: file_handle
    error stop "write_attribute_array_1d_real should not be called - &
      & use concrete implementation"
  end subroutine write_attribute_array_1d_real

end module m_io_base
