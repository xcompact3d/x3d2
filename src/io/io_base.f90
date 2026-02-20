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
!! - provides abstract `reader_base_t` and `writer_base_t` types
!! - enforces a consistent interface for all backends
!!
!! @note This is an internal module and should not be used directly by users.
!! The sole public interface for I/O is the high-level session API provided in
!! `m_io_session`.

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
    procedure :: is_file_functional => base_is_file_functional
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
    procedure :: supports_device_field_write => &
      base_supports_device_field_write
    procedure :: sync_device => base_sync_device
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, &
      write_data_array_3d
    procedure :: write_data_i8
    procedure :: write_data_integer
    procedure :: write_data_real
    procedure :: write_data_array_3d
    procedure :: write_field_from_solver
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

  logical function base_supports_device_field_write(self)
    class(io_writer_t), intent(in) :: self
    base_supports_device_field_write = .false.
  end function base_supports_device_field_write

  subroutine base_sync_device(self)
    !! Ensure all device operations complete before I/O.
    !! Does nothing for non-GPU backends. CUDA-aware backends override
    !! this to call cudaDeviceSynchronize.
    class(io_writer_t), intent(inout) :: self
  end subroutine base_sync_device

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

  subroutine write_field_from_solver( &
    self, variable_name, field, file_handle, backend, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    !! Write field data directly, with backend-specific optimisations
    !! Each backend implements this to use GPU-aware I/O when available
    !! NOTE: field and backend are polymorphic - concrete types defined in implementations
    class(io_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    class(*), intent(in) :: field
    class(io_file_t), intent(inout) :: file_handle
    class(*), intent(in) :: backend
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp
    error stop "write_field_from_solver should not be called - &
      & use concrete implementation"
  end subroutine write_field_from_solver

end module m_io_base
