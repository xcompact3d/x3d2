module m_io_base
  !! Abstract base types and interfaces for session-based I/O architecture.
  !!
  !! This internal module defines the fundamental building blocks of the I/O
  !! system. It establishes a polymorphic layer that allows the high-level
  !! user session to interact with various I/O backends (e.g., ADIOS2, dummy)
  !! through a consistent interface.
  !!
  !! **Architecture Layers:**
  !!
  !! 1. **User Code** - interacts only with the Session layer
  !!
  !! 2. **Session Layer** (`m_io_session`)
  !!    - Manages all I/O complexity (file handles, state, etc.)
  !!    - Instantiates the I/O backend selected at compile-time
  !!    - Provides `reader_session_t` and `writer_session_t` for users
  !!
  !! 3. **Backend Layer** (`m_io_backend`)
  !!    - Concrete implementation of an I/O backend (e.g., ADIOS2)
  !!    - Extends the abstract base types defined in this module
  !!
  !! 4. **Base Layer** (`m_io_base`, this module)
  !!    - Provides abstract `io_reader_t` and `io_writer_t` types
  !!    - Enforces a consistent interface for all backends
  !!
  !! **Note:** This is an internal module and should not be used directly by
  !! users. The sole public interface for I/O is the high-level session API
  !! provided in `m_io_session`.

  use m_common, only: dp, i8

  implicit none

  private
  public :: io_reader_t, io_writer_t, io_file_t
  public :: io_mode_read, io_mode_write

  integer, parameter :: io_mode_read = 1   !! Read mode flag for opening files
  integer, parameter :: io_mode_write = 2  !! Write mode flag for opening files

  type :: io_file_t
    !! Base file handle for I/O operations.
    !!
    !! This abstract type represents an open file handle. Concrete backends
    !! extend this type to implement backend-specific file operations.
    !! Provides step-based I/O for time-series data.
  contains
    procedure :: close => base_close                     !! Close the file
    procedure :: begin_step => base_begin_step           !! Begin a new I/O step
    procedure :: end_step => base_end_step               !! End current I/O step
    procedure :: is_file_functional => base_is_file_functional !! Check if file is operational
  end type io_file_t

  type :: io_reader_t
    !! Base I/O reader type for polymorphic usage.
    !!
    !! This abstract type provides the interface for reading data from files.
    !! Concrete backends (e.g., ADIOS2) extend this type to implement
    !! backend-specific reading operations. Supports reading scalars and
    !! 3D arrays with optional hyperslab selection.
  contains
    procedure :: init => base_reader_init                !! Initialise reader
    procedure :: open => base_reader_open                !! Open file for reading
    procedure :: finalise => base_reader_finalise        !! Finalise and clean up
    ! Generic interfaces for session usage
    generic :: read_data => read_data_i8, read_data_integer, read_data_real, &
      read_data_array_3d                                 !! Read data (generic interface)
    procedure :: read_data_i8                            !! Read 64-bit integer
    procedure :: read_data_integer                       !! Read default integer
    procedure :: read_data_real                          !! Read double precision real
    procedure :: read_data_array_3d                      !! Read 3D array
  end type io_reader_t

  type :: io_writer_t
    !! Base I/O writer type for polymorphic usage.
    !!
    !! This abstract type provides the interface for writing data to files.
    !! Concrete backends (e.g., ADIOS2) extend this type to implement
    !! backend-specific writing operations. Supports writing scalars,
    !! 3D arrays, and attributes.
  contains
    procedure :: init => base_writer_init                !! Initialise writer
    procedure :: open => base_writer_open                !! Open file for writing
    procedure :: finalise => base_writer_finalise        !! Finalise and clean up
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, &
      write_data_array_3d                                !! Write data (generic interface)
    procedure :: write_data_i8                           !! Write 64-bit integer
    procedure :: write_data_integer                      !! Write default integer
    procedure :: write_data_real                         !! Write double precision real
    procedure :: write_data_array_3d                     !! Write 3D array
    generic :: write_attribute => write_attribute_string, &
      write_attribute_array_1d_real                      !! Write attribute (generic interface)
    procedure :: write_attribute_string                  !! Write string attribute
    procedure :: write_attribute_array_1d_real           !! Write 1D real array attribute
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
    class(io_file_t), allocatable :: file_handle
    type(io_file_t) :: temp_handle
    file_handle = temp_handle
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
    class(io_file_t), allocatable :: file_handle
    type(io_file_t) :: temp_handle
    file_handle = temp_handle
    error stop "base_writer_open should not be called - &
      & use concrete implementation"
  end function base_writer_open

  subroutine base_writer_finalise(self)
    class(io_writer_t), intent(inout) :: self
    error stop "base_writer_finalise should not be called - &
      & use concrete implementation"
  end subroutine base_writer_finalise

  function base_is_file_functional(self) result(is_functional)
    class(io_file_t), intent(in) :: self
    logical :: is_functional
    is_functional = .true.
  end function base_is_file_functional

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

  subroutine write_data_real(self, variable_name, value, file_handle, use_sp)
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    logical, intent(in), optional :: use_sp
    error stop "write_data_real should not be called - &
      & use concrete implementation"
  end subroutine write_data_real

  subroutine write_data_array_3d( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp
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
