module m_io_base
!! This module provides a general-purpose I/O abstraction that can be used
!! by any module requiring file I/O operations. It defines abstract interfaces
!! for I/O operations and re-exports the factory function for creating readers.

  use m_common, only: dp, i8

  implicit none

  private
  public :: io_reader_t, io_writer_t, io_file_t
  public :: io_mode_read, io_mode_write

  integer, parameter :: io_mode_read = 1
  integer, parameter :: io_mode_write = 2

  type, abstract :: io_file_t
  contains
    procedure(file_close), deferred :: close
    procedure(file_begin_step), deferred :: begin_step
    procedure(file_end_step), deferred :: end_step
  end type io_file_t

  type, abstract :: io_reader_t
  contains
    procedure(reader_init), deferred :: init
    procedure(reader_open), deferred :: open
    procedure(reader_finalise), deferred :: finalise
    
    generic :: read_data => read_data_i8, read_data_integer, read_data_real, read_data_array_3d
    procedure(read_data_i8), deferred :: read_data_i8
    procedure(read_data_integer), deferred :: read_data_integer
    procedure(read_data_real), deferred :: read_data_real
    procedure(read_data_array_3d), deferred :: read_data_array_3d
  end type io_reader_t

  type, abstract :: io_writer_t
  contains
    procedure(writer_init), deferred :: init
    procedure(writer_open), deferred :: open
    procedure(writer_finalise), deferred :: finalise
    
    generic :: write_data => write_data_i8, write_data_integer, write_data_real, write_data_array_3d
    procedure(write_data_i8), deferred :: write_data_i8
    procedure(write_data_integer), deferred :: write_data_integer
    procedure(write_data_real), deferred :: write_data_real
    procedure(write_data_array_3d), deferred :: write_data_array_3d

    generic :: write_attribute => write_attribute_string, write_attribute_array_1d_real
    procedure(write_attribute_string), deferred :: write_attribute_string
    procedure(write_attribute_array_1d_real), deferred :: write_attribute_array_1d_real
  end type io_writer_t

  abstract interface
    subroutine file_close(self)
      import :: io_file_t
      class(io_file_t), intent(inout) :: self
    end subroutine file_close

    subroutine file_begin_step(self)
      import :: io_file_t
      class(io_file_t), intent(inout) :: self
    end subroutine file_begin_step

    subroutine file_end_step(self)
      import :: io_file_t
      class(io_file_t), intent(inout) :: self
    end subroutine file_end_step

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
      class(io_file_t), pointer :: file_handle
    end function reader_open

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

    subroutine read_data_array_3d(self, variable_name, array, file_handle, shape_dims, start_dims, count_dims)
      import :: io_reader_t, io_file_t, dp, i8
      class(io_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(inout) :: array(:, :, :)
      class(io_file_t), intent(inout) :: file_handle
      integer(i8), intent(in), optional :: shape_dims(3)
      integer(i8), intent(in), optional :: start_dims(3)
      integer(i8), intent(in), optional :: count_dims(3)
    end subroutine read_data_array_3d

    subroutine reader_finalise(self)
      import :: io_reader_t
      class(io_reader_t), intent(inout) :: self
    end subroutine reader_finalise

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
      class(io_file_t), pointer :: file_handle
    end function writer_open

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

    subroutine write_data_real(self, variable_name, value, file_handle)
      import :: io_writer_t, io_file_t, dp
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_data_real

    subroutine write_data_array_3d(self, variable_name, array, file_handle, shape_dims, start_dims, count_dims)
      import :: io_writer_t, io_file_t, dp, i8
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: variable_name
      real(dp), intent(in) :: array(:, :, :)
      class(io_file_t), intent(inout) :: file_handle
      integer(i8), intent(in), optional :: shape_dims(3)
      integer(i8), intent(in), optional :: start_dims(3)
      integer(i8), intent(in), optional :: count_dims(3)
    end subroutine write_data_array_3d

    subroutine writer_finalise(self)
      import :: io_writer_t
      class(io_writer_t), intent(inout) :: self
    end subroutine writer_finalise

    subroutine write_attribute_string(self, attribute_name, value, file_handle)
      import :: io_writer_t, io_file_t
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: attribute_name
      character(len=*), intent(in) :: value
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_attribute_string

    subroutine write_attribute_array_1d_real(self, attribute_name, values, file_handle)
      import :: io_writer_t, io_file_t, dp
      class(io_writer_t), intent(inout) :: self
      character(len=*), intent(in) :: attribute_name
      real(dp), intent(in) :: values(:)
      class(io_file_t), intent(inout) :: file_handle
    end subroutine write_attribute_array_1d_real
  end interface

end module m_io_base
