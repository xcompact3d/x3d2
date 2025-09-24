module m_io_dummy
!! Dummy implementation of the general I/O interface for when ADIOS2 is not available

  use iso_fortran_env, only: stderr => error_unit
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, io_mode_read, io_mode_write
  use m_common, only: dp, i8

  implicit none

  private
  public :: io_dummy_reader_t, io_dummy_writer_t

  type, extends(io_file_t) :: io_dummy_file_t
    logical :: is_open = .false.
  contains
    procedure :: close => file_close_dummy
    procedure :: begin_step => file_begin_step_dummy
    procedure :: end_step => file_end_step_dummy
  end type io_dummy_file_t

  type, extends(io_reader_t) :: io_dummy_reader_t
    logical :: initialised = .false.
  contains
    procedure :: init => reader_init_dummy
    procedure :: open => reader_open_dummy
    procedure :: finalise => reader_finalise_dummy
    
    generic :: read_data => read_data_i8_dummy, read_data_integer_dummy, read_data_real_dummy, read_data_array_3d_dummy
    procedure :: read_data_i8 => read_data_i8_dummy
    procedure :: read_data_integer => read_data_integer_dummy
    procedure :: read_data_real => read_data_real_dummy
    procedure :: read_data_array_3d => read_data_array_3d_dummy
  end type io_dummy_reader_t

  type, extends(io_writer_t) :: io_dummy_writer_t
    logical :: initialised = .false.
  contains
    procedure :: init => writer_init_dummy
    procedure :: open => writer_open_dummy
    procedure :: finalise => writer_finalise_dummy
    
    generic :: write_data => write_data_i8_dummy, write_data_integer_dummy, write_data_real_dummy, write_data_array_3d_dummy
    procedure :: write_data_i8 => write_data_i8_dummy
    procedure :: write_data_integer => write_data_integer_dummy
    procedure :: write_data_real => write_data_real_dummy
    procedure :: write_data_array_3d => write_data_array_3d_dummy

    generic :: write_attribute => write_attribute_string_dummy, write_attribute_array_1d_real_dummy
    procedure :: write_attribute_string => write_attribute_string_dummy
    procedure :: write_attribute_array_1d_real => write_attribute_array_1d_real_dummy
  end type io_dummy_writer_t

contains

  subroutine file_close_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    write(stderr, '(A)') "WARNING: Dummy I/O backend - file close operation ignored"
    self%is_open = .false.
  end subroutine file_close_dummy

  subroutine file_begin_step_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    write(stderr, '(A)') "WARNING: Dummy I/O backend - begin_step operation ignored"
  end subroutine file_begin_step_dummy

  subroutine file_end_step_dummy(self)
    class(io_dummy_file_t), intent(inout) :: self
    write(stderr, '(A)') "WARNING: Dummy I/O backend - end_step operation ignored"
  end subroutine file_end_step_dummy

  subroutine reader_init_dummy(self, comm, name)
    class(io_dummy_reader_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    write(stderr, '(A)') "WARNING: Using dummy I/O backend - ADIOS2 not available"
    write(stderr, '(A)') "I/O operations will not function properly"
    self%initialised = .true.
  end subroutine reader_init_dummy

  function reader_open_dummy(self, filename, mode, comm) result(file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), pointer :: file_handle

    type(io_dummy_file_t), allocatable :: dummy_file

    write(stderr, '(A)') "ERROR: Cannot open file '" // trim(filename) // "' - ADIOS2 not available"
    write(stderr, '(A)') "Compilation was done without ADIOS2 support"
    write(stderr, '(A)') "Please recompile with -DWITH_ADIOS2=ON to enable I/O functionality"

    allocate(io_dummy_file_t :: file_handle)
    select type(file_handle)
    type is (io_dummy_file_t)
      file_handle%is_open = .false.
    end select
  end function reader_open_dummy

  subroutine read_data_i8_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot read scalar i8 '" // trim(variable_name) // "' - ADIOS2 not available"
    value = 0_i8
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine read_data_i8_dummy

  subroutine read_data_integer_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot read scalar integer '" // trim(variable_name) // "' - ADIOS2 not available"
    value = 0
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine read_data_integer_dummy

  subroutine read_data_real_dummy(self, variable_name, value, file_handle)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot read scalar real '" // trim(variable_name) // "' - ADIOS2 not available"
    value = 0.0_dp
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine read_data_real_dummy

  subroutine read_data_array_3d_dummy(self, variable_name, array, file_handle, shape_dims, start_dims, count_dims)
    class(io_dummy_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    write(stderr, '(A)') "ERROR: Cannot read 3D array '" // trim(variable_name) // "' - ADIOS2 not available"
    array = 0.0_dp
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine read_data_array_3d_dummy

  subroutine reader_finalise_dummy(self)
    class(io_dummy_reader_t), intent(inout) :: self
    write(stderr, '(A)') "WARNING: Dummy I/O backend - finalise operation ignored"
    self%initialised = .false.
  end subroutine reader_finalise_dummy

  subroutine writer_init_dummy(self, comm, name)
    class(io_dummy_writer_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name
    write(stderr, '(A)') "WARNING: Using dummy I/O backend - ADIOS2 not available"
    write(stderr, '(A)') "I/O operations will not function properly"
    self%initialised = .true.
  end subroutine writer_init_dummy

  function writer_open_dummy(self, filename, mode, comm) result(file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm
    class(io_file_t), pointer :: file_handle

    write(stderr, '(A)') "ERROR: Cannot open file '" // trim(filename) // "' for writing - ADIOS2 not available"
    write(stderr, '(A)') "Compilation was done without ADIOS2 support"
    write(stderr, '(A)') "Please recompile with -DWITH_ADIOS2=ON to enable I/O functionality"

    allocate(io_dummy_file_t :: file_handle)
    select type(file_handle)
    type is (io_dummy_file_t)
      file_handle%is_open = .false.
    end select
  end function writer_open_dummy

  subroutine write_data_i8_dummy(self, variable_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot write scalar i8 '" // trim(variable_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_data_i8_dummy

  subroutine write_data_integer_dummy(self, variable_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot write scalar integer '" // trim(variable_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_data_integer_dummy

  subroutine write_data_real_dummy(self, variable_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot write scalar real '" // trim(variable_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_data_real_dummy

  subroutine write_data_array_3d_dummy(self, variable_name, array, file_handle, shape_dims, start_dims, count_dims)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    write(stderr, '(A)') "ERROR: Cannot write 3D array '" // trim(variable_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_data_array_3d_dummy

  subroutine writer_finalise_dummy(self)
    class(io_dummy_writer_t), intent(inout) :: self
    write(stderr, '(A)') "WARNING: Dummy I/O backend - finalise operation ignored"
    self%initialised = .false.
  end subroutine writer_finalise_dummy

  subroutine write_attribute_string_dummy(self, attribute_name, value, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot write string attribute '" // trim(attribute_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_attribute_string_dummy

  subroutine write_attribute_array_1d_real_dummy(self, attribute_name, values, file_handle)
    class(io_dummy_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    real(dp), intent(in) :: values(:)
    class(io_file_t), intent(inout) :: file_handle

    write(stderr, '(A)') "ERROR: Cannot write real array attribute '" // trim(attribute_name) // "' - ADIOS2 not available"
    error stop "I/O operation failed: ADIOS2 support not compiled"
  end subroutine write_attribute_array_1d_real_dummy

end module m_io_dummy
