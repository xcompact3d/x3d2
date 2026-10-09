module m_io_base
!! @brief Provides the abstract base types and interfaces for the session-based
!! I/O architecture.
!!
!! @details This internal module defines the fundamental building blocks of
!! the I/O system. It establishes a polymorphic layer that allows the
!! high-level user session to interact with various I/O backends through a
!! consistent interface.
!!
!! The architecture is designed in distinct layers:
!! User code
!! - interacts only with the Session layer
!!
!! Session layer (`m_io_session`)
!! - manages all I/O complexity (file handles, state, etc.)
!! - instantiates the I/O backend selected at compile-time
!! - provides `reader_session_t` and `writer_session_t` for users
!!
!! Backend layer (`m_io_backend`)
!! - concrete implementation of an I/O backed (e.g., ADIOS2)
!! - extends the abstract base types defined in this module
!!
!! Base layer (`m_io_base`, this module)
!! - provides abstract `io_reader_t`, `io_writer_t` and `io_file_t` types
!! - deferred bindings make every backend implement the full interface at
!!   compile time; only optional capabilities have default implementations
!!
!! @note This is an internal module and should not be used directly by users.
!! The sole public interface for I/O is the high-level session API provided in
!! `m_io_session`.

  use m_common, only: dp, i8
  use m_field, only: field_t
  use m_base_backend, only: base_backend_t

  implicit none

  private
  public :: io_reader_t, io_writer_t, io_file_t
  public :: io_mode_read, io_mode_write

  integer, parameter :: io_mode_read = 1
  integer, parameter :: io_mode_write = 2

  !> Base file handle for I/O operations
  type, abstract :: io_file_t
  contains
    procedure(file_op), deferred :: close
    procedure(file_op), deferred :: begin_step
    procedure(file_op), deferred :: end_step
    procedure :: is_file_functional => base_is_file_functional
  end type io_file_t

  !> Base I/O reader type for polymorphic usage
  type, abstract :: io_reader_t
  contains
    procedure(reader_init), deferred :: init
    procedure(reader_open), deferred :: open
    procedure(reader_finalise), deferred :: finalise
    ! Generic interfaces for session usage
    generic :: read_data => read_data_i8, read_data_integer, read_data_real, &
      read_data_array_3d
    procedure(read_data_i8), deferred :: read_data_i8
    procedure(read_data_integer), deferred :: read_data_integer
    procedure(read_data_real), deferred :: read_data_real
    procedure(read_data_array_3d), deferred :: read_data_array_3d
  end type io_reader_t

  !> Base I/O writer type for polymorphic usage
  type, abstract :: io_writer_t
  contains
    procedure(writer_init), deferred :: init
    procedure(writer_open), deferred :: open
    procedure(writer_finalise), deferred :: finalise
    procedure :: supports_device_field_write => &
      base_supports_device_field_write
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, &
      write_data_array_3d
    procedure(write_data_i8), deferred :: write_data_i8
    procedure(write_data_integer), deferred :: write_data_integer
    procedure(write_data_real), deferred :: write_data_real
    procedure(write_data_array_3d), deferred :: write_data_array_3d
    procedure(write_field_from_solver), deferred :: write_field_from_solver
    generic :: write_attribute => write_attribute_string, &
      write_attribute_array_1d_real
    procedure(write_attribute_string), deferred :: write_attribute_string
    procedure(write_attribute_array_1d_real), deferred :: &
      write_attribute_array_1d_real
  end type io_writer_t

  abstract interface
    subroutine file_op(self)
      import :: io_file_t
      class(io_file_t), intent(inout) :: self
    end subroutine file_op

    subroutine reader_init(self, comm, name)
      import :: io_reader_t
      class(io_reader_t), intent(inout) :: self
      integer, intent(in) :: comm
      character(len=*), intent(in) :: name
    end subroutine reader_init

    function reader_open(self, filename, mode, comm) result(file_handle)
      import :: io_reader_t, io_file_t
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: filename
      integer, intent(in) :: mode
      integer, intent(in) :: comm
      class(io_file_t), allocatable :: file_handle
    end function reader_open

    subroutine reader_finalise(self)
      import :: io_reader_t
      class(io_reader_t), intent(inout) :: self
    end subroutine reader_finalise

    subroutine read_data_i8(self, variable_name, value, file_handle)
      import :: io_reader_t, io_file_t, i8
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      integer(i8), intent(out) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine read_data_i8

    subroutine read_data_integer(self, variable_name, value, file_handle)
      import :: io_reader_t, io_file_t
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      integer, intent(out) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine read_data_integer

    subroutine read_data_real(self, variable_name, value, file_handle)
      import :: io_reader_t, io_file_t, dp
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(out) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine read_data_real

    subroutine read_data_array_3d( &
      self, variable_name, array, file_handle, &
      shape_dims, start_dims, count_dims &
      )
      import :: io_reader_t, io_file_t, dp, i8
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(inout) :: array(:, :, :)
      class(io_file_t), intent(inout) :: file_handle
      integer(i8), intent(in), optional :: shape_dims(3)
      integer(i8), intent(in), optional :: start_dims(3)
      integer(i8), intent(in), optional :: count_dims(3)
    end subroutine read_data_array_3d

    subroutine writer_init(self, comm, name)
      import :: io_writer_t
      class(io_writer_t), intent(inout) :: self
      integer, intent(in) :: comm
      character(len=*), intent(in) :: name
    end subroutine writer_init

    function writer_open(self, filename, mode, comm) result(file_handle)
      import :: io_writer_t, io_file_t
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: filename
      integer, intent(in) :: mode
      integer, intent(in) :: comm
      class(io_file_t), allocatable :: file_handle
    end function writer_open

    subroutine writer_finalise(self)
      import :: io_writer_t
      class(io_writer_t), intent(inout) :: self
    end subroutine writer_finalise

    subroutine write_data_i8(self, variable_name, value, file_handle)
      import :: io_writer_t, io_file_t, i8
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      integer(i8), intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_data_i8

    subroutine write_data_integer(self, variable_name, value, file_handle)
      import :: io_writer_t, io_file_t
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      integer, intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_data_integer

    subroutine write_data_real(self, variable_name, value, file_handle, &
                               use_sp)
      import :: io_writer_t, io_file_t, dp
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
      logical, intent(in), optional :: use_sp
    end subroutine write_data_real

    subroutine write_data_array_3d( &
      self, variable_name, array, file_handle, &
      shape_dims, start_dims, count_dims, use_sp &
      )
      import :: io_writer_t, io_file_t, dp, i8
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(in) :: array(:, :, :)
      class(io_file_t), intent(inout) :: file_handle
      integer(i8), intent(in) :: shape_dims(3)
      integer(i8), intent(in) :: start_dims(3)
      integer(i8), intent(in) :: count_dims(3)
      logical, intent(in), optional :: use_sp
    end subroutine write_data_array_3d

    subroutine write_field_from_solver( &
      self, variable_name, field, file_handle, backend, &
      shape_dims, start_dims, count_dims, use_sp &
      )
      !! Write the first count_dims points of a solver field. Backends
      !! may write device-resident fields without staging through host
      !! memory; `backend` provides the layout-aware packing for that.
      import :: io_writer_t, io_file_t, field_t, base_backend_t, i8
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      class(field_t), intent(in) :: field
      class(io_file_t), intent(inout) :: file_handle
      class(base_backend_t), intent(inout) :: backend
      integer(i8), intent(in) :: shape_dims(3)
      integer(i8), intent(in) :: start_dims(3)
      integer(i8), intent(in) :: count_dims(3)
      logical, intent(in), optional :: use_sp
    end subroutine write_field_from_solver

    subroutine write_attribute_string( &
      self, attribute_name, value, file_handle &
      )
      import :: io_writer_t, io_file_t
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: attribute_name
      character(len=*), intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_attribute_string

    subroutine write_attribute_array_1d_real( &
      self, attribute_name, values, file_handle &
      )
      import :: io_writer_t, io_file_t, dp
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: attribute_name
      real(dp), intent(in) :: values(:)
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_attribute_array_1d_real
  end interface

contains

  logical function base_supports_device_field_write(self, field, backend)
    !! Default capability: fields are always staged through host memory.
    class(io_writer_t), intent(in) :: self
    class(field_t), intent(in) :: field
    class(base_backend_t), intent(in) :: backend
    base_supports_device_field_write = .false.
  end function base_supports_device_field_write

  function base_is_file_functional(self) result(is_functional)
    class(io_file_t), intent(in) :: self
    logical :: is_functional
    is_functional = .true.
  end function base_is_file_functional

end module m_io_base
