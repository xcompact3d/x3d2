module m_io_backend
  !! Dummy (non-functional) I/O backend for when no real backend is available.
  !!
  !! This module provides a fallback implementation of the I/O backend
  !! interface used when no real I/O backend (e.g., ADIOS2) is enabled at
  !! compile time. It allows the code to compile and link without a functional
  !! I/O library.
  !!
  !! **Purpose:**
  !!
  !! - Enables compilation without external I/O library dependencies
  !! - Provides informative error messages when I/O operations are attempted
  !! - Allows code structure to remain consistent regardless of I/O backend
  !!
  !! **Behaviour:**
  !!
  !! - Write operations are silently ignored (no-op)
  !! - Read operations terminate with error message directing user to recompile
  !! - File open/close operations are tracked but perform no actual I/O
  !!
  !! **Use Cases:**
  !!
  !! - Testing/debugging without I/O overhead
  !! - Systems where ADIOS2 is unavailable
  !! - Dry runs to validate simulation setup
  !!
  !! **Warning:** This is a non-functional stub. If you require actual file I/O,
  !! recompile with `-DWITH_ADIOS2=ON` to enable the ADIOS2 backend.
  use iso_fortran_env, only: stderr => error_unit
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, io_mode_read, &
                       io_mode_write
  use m_common, only: dp, i8

  implicit none

  private
  public :: allocate_io_reader, allocate_io_writer
  public :: get_default_backend, IO_BACKEND_DUMMY, IO_BACKEND_ADIOS2

  logical, save :: write_warning_shown = .false. !! Track if warning has been displayed

  integer, parameter :: IO_BACKEND_DUMMY = 0   !! Dummy backend identifier
  integer, parameter :: IO_BACKEND_ADIOS2 = 1  !! ADIOS2 backend identifier

  type, extends(io_file_t) :: io_dummy_file_t
    !! Dummy file handle (tracks state but performs no I/O).
    logical :: is_open = .false. !! File open state flag
  contains
    procedure :: close => file_close_dummy                  !! Close file (no-op)
    procedure :: begin_step => file_begin_step_dummy        !! Begin step (no-op)
    procedure :: end_step => file_end_step_dummy            !! End step (no-op)
    procedure :: is_file_functional => is_file_functional_dummy !! Check if functional
  end type io_dummy_file_t

  type, extends(io_reader_t) :: io_dummy_reader_t
    !! Dummy reader (errors on read attempts).
    logical :: initialised = .false. !! Initialisation state flag
  contains
    procedure :: init => reader_init_dummy                  !! Initialise reader
    procedure :: open => reader_open_dummy                  !! Open file (returns non-functional handle)
    procedure :: finalise => reader_finalise_dummy          !! Finalise (no-op)
    procedure :: read_data_i8 => read_data_i8_dummy         !! Read i8 (errors)
    procedure :: read_data_integer => read_data_integer_dummy !! Read integer (errors)
    procedure :: read_data_real => read_data_real_dummy     !! Read real (errors)
    procedure :: read_data_array_3d => read_data_array_3d_dummy !! Read 3D array (errors)
  end type io_dummy_reader_t

  type, extends(io_writer_t) :: io_dummy_writer_t
    !! Dummy writer (silently ignores write operations).
    logical :: initialised = .false. !! Initialisation state flag
  contains
    procedure :: init => writer_init_dummy                  !! Initialise writer
    procedure :: open => writer_open_dummy                  !! Open file (returns non-functional handle)
    procedure :: finalise => writer_finalise_dummy          !! Finalise (no-op)
    procedure :: write_data_i8 => write_data_i8_dummy       !! Write i8 (no-op)
    procedure :: write_data_integer => write_data_integer_dummy !! Write integer (no-op)
    procedure :: write_data_real => write_data_real_dummy   !! Write real (no-op)
    procedure :: write_data_array_3d => write_data_array_3d_dummy !! Write 3D array (no-op)
    procedure :: write_attribute_string => write_attribute_string_dummy !! Write string attribute (no-op)
    procedure :: write_attribute_array_1d_real => &
      write_attribute_array_1d_real_dummy                   !! Write 1D real array attribute (no-op)
  end type io_dummy_writer_t

contains

  subroutine allocate_io_reader(reader)
    class(io_reader_t), allocatable, intent(out) :: reader
    allocate (io_dummy_reader_t :: reader)
  end subroutine allocate_io_reader

  subroutine allocate_io_writer(writer)
    class(io_writer_t), allocatable, intent(out) :: writer
    allocate (io_dummy_writer_t :: writer)
  end subroutine allocate_io_writer

  function get_default_backend() result(backend)
    integer :: backend
    backend = IO_BACKEND_DUMMY
  end function get_default_backend

  subroutine report_read_error(variable_name)
    character(len=*), intent(in) :: variable_name
    error stop "ERROR: Cannot read '"//trim(variable_name)// &
      "' - recompile with -DWITH_ADIOS2=ON"
  end subroutine report_read_error

  subroutine file_close_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    ! silently ignore file close operations
    self%is_open = .false.
  end subroutine file_close_dummy

  subroutine file_begin_step_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    ! silently ignore begin_step operations
  end subroutine file_begin_step_dummy

  subroutine file_end_step_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    ! silently ignore end_step operations
  end subroutine file_end_step_dummy

  subroutine reader_init_dummy(self, comm, name)
    class(io_dummy_reader_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    self%initialised = .true.
  end subroutine reader_init_dummy

  function reader_open_dummy(self, filename, mode, comm) result(file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), allocatable :: file_handle

    write (stderr, '(A)') "ERROR: Cannot open file '"//trim(filename)// &
      "' for reading - ADIOS2 not available"
    write (stderr, '(A)') "File reading requires ADIOS2 support"
    write (stderr, '(A)') "Please recompile with -DWITH_ADIOS2=ON"

    allocate (io_dummy_file_t :: file_handle)
    select type (file_handle)
    type is (io_dummy_file_t)
      file_handle%is_open = .false.
    end select
  end function reader_open_dummy

  function is_file_functional_dummy(self) result(is_functional)
    class(io_dummy_file_t), intent(in) :: self
    logical :: is_functional
    is_functional = self%is_open
  end function is_file_functional_dummy

  subroutine read_data_i8_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    value = 0_i8
    call report_read_error(variable_name)
  end subroutine read_data_i8_dummy

  subroutine read_data_integer_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    value = 0
    call report_read_error(variable_name)
  end subroutine read_data_integer_dummy

  subroutine read_data_real_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    value = 0.0_dp
    call report_read_error(variable_name)
  end subroutine read_data_real_dummy

  subroutine read_data_array_3d_dummy( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims &
    )
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    array = 0.0_dp
    call report_read_error(variable_name)
  end subroutine read_data_array_3d_dummy

  subroutine reader_finalise_dummy(self)
    class(io_dummy_reader_t), intent(inout) :: self
    ! silently finalise
    self%initialised = .false.
  end subroutine reader_finalise_dummy

  subroutine writer_init_dummy(self, comm, name)
    class(io_dummy_writer_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    self%initialised = .true.
  end subroutine writer_init_dummy

  function writer_open_dummy(self, filename, mode, comm) result(file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), allocatable :: file_handle

    ! Show warning for write operations (once)
    if (.not. write_warning_shown) then
      write (stderr, '(A)') "WARNING: Cannot save file '"//trim(filename)// &
        "' - ADIOS2 not available"
      write (stderr, '(A)') "Checkpoints and snapshots will not be written"
      write_warning_shown = .true.
    end if

    ! Silently create dummy file handle
    allocate (io_dummy_file_t :: file_handle)
    select type (file_handle)
    type is (io_dummy_file_t)
      file_handle%is_open = .false.
    end select
  end function writer_open_dummy

  subroutine write_data_i8_dummy(self, variable_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    ! silently ignore write operations
  end subroutine write_data_i8_dummy

  subroutine write_data_integer_dummy(self, variable_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    ! silently ignore write operations
  end subroutine write_data_integer_dummy

  subroutine write_data_real_dummy(self, variable_name, value, file_handle, &
                                   use_sp)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    logical, intent(in), optional :: use_sp
    ! silently ignore write operations
  end subroutine write_data_real_dummy

  subroutine write_data_array_3d_dummy( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp
    ! silently ignore write operations
  end subroutine write_data_array_3d_dummy

  subroutine writer_finalise_dummy(self)
    class(io_dummy_writer_t), intent(inout) :: self
    ! silently finalise
    self%initialised = .false.
  end subroutine writer_finalise_dummy

  subroutine write_attribute_string_dummy( &
    self, attribute_name, value, file_handle &
    )
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    ! silently ignore attribute writes
  end subroutine write_attribute_string_dummy

  subroutine write_attribute_array_1d_real_dummy( &
    self, attribute_name, values, file_handle &
    )
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    real(dp), intent(in) :: values(:)
    class(io_file_t), intent(inout) :: file_handle
    ! silently ignore attribute writes
  end subroutine write_attribute_array_1d_real_dummy

end module m_io_backend
