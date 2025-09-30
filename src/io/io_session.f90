module m_io_session
!! **PRIMARY USER INTERFACE FOR X3D2 I/O OPERATIONS**
!!
!! This module provides the ONLY interface that users need for I/O operations in X3D2.
!! All file reading and writing should be done through the session types provided here.
!!
!! **Key Features:**
!! - Type-safe specialised sessions: `reader_session_t` for reading, `writer_session_t` for writing
!! - Automatic backend selection (currently ADIOS2, extensible to other formats)
!! - Simple session-based workflow: open -> read/write -> close
!! - No manual file handle management required
!! - Compile-time prevention of mixing read/write operations
!!
!! **Usage Pattern:**
!! ```fortran
!! use m_io_session, only: writer_session_t, reader_session_t
!!
!! ! For writing data
!! type(writer_session_t) :: writer_session
!! call writer_session%open("output.bp", comm)
!! call writer_session%write_data("temperature", temp_field)
!! call writer_session%close()
!!
!! ! For reading data
!! type(reader_session_t) :: reader_session
!! call reader_session%open("input.bp", comm)
!! call reader_session%read_data("temperature", temp_field)
!! call reader_session%close()
!! ```
!!
!! **Note:** Users should NOT directly use lower-level modules like `m_io_base`,
!! `m_io_factory`, or backend-specific modules. This module abstracts away all
!! the complexity and provides everything needed for I/O operations.

  use m_common, only: dp, i8
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, &
                       io_mode_read, io_mode_write
  use m_io_factory, only: allocate_io_reader, allocate_io_writer

  implicit none

  private

  !! Public session types for user interaction
  public :: reader_session_t, writer_session_t

  !> Base type for common session functionality
  type :: io_session_base_t
    private
    class(io_file_t), pointer :: file => null()
    logical :: is_open = .false.
    logical :: is_functional = .true.  ! false for dummy I/O
  contains
    procedure :: is_session_open
    procedure :: is_session_functional
    procedure :: close => session_base_close
    procedure :: finalise => session_base_finalise
    procedure :: get_file
  end type io_session_base_t

  !> **PRIMARY TYPE FOR READING DATA** - Use this for all file reading operations
  !! This is the only interface users should use for reading data.
  !! Provides type-safe reading operations with automatic backend selection.
  !!
  !! Usage example:
  !!   type(reader_session_t) :: reader_session
  !!   call reader_session%open("checkpoint.bp", MPI_COMM_WORLD)
  !!   call reader_session%read_data("timestep", timestep)
  !!   call reader_session%read_data("velocity_u", u_field, start_dims, count_dims)
  !!   call reader_session%close()
  type, extends(io_session_base_t) :: reader_session_t
    private
    class(io_reader_t), pointer :: reader => null()
  contains
    ! Open/close operations
    procedure :: open => reader_session_open
    ! Generic read_data interface
    generic :: read_data => read_data_i8, read_data_integer, &
      read_data_real, read_data_array_3d
    procedure, private :: read_data_i8
    procedure, private :: read_data_integer
    procedure, private :: read_data_real
    procedure, private :: read_data_array_3d
  end type reader_session_t

  !> **PRIMARY TYPE FOR WRITING DATA** - Use this for all file writing operations
  !! This is the only interface users should use for writing data.
  !! Provides type-safe writing operations with automatic backend selection.
  !!
  !! Usage example:
  !!   type(writer_session_t) :: writer_session
  !!   call writer_session%open("output.bp", MPI_COMM_WORLD)
  !!   call writer_session%write_data("timestep", current_step)
  !!   call writer_session%write_data("pressure", p_field, start_dims, count_dims)
  !!   call writer_session%close()
  !!   call writer_session%write_attribute("ParaView", "vtk_xml_content")
  !!   call writer_session%close()
  type, extends(io_session_base_t) :: writer_session_t
    private
    class(io_writer_t), pointer :: writer => null()
  contains
    ! Open/close operations
    procedure :: open => writer_session_open
    ! Generic write_data interface
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, write_data_array_3d
    procedure, private :: write_data_i8
    procedure, private :: write_data_integer
    procedure, private :: write_data_real
    procedure, private :: write_data_array_3d
    ! Write attribute interface
    procedure :: write_attribute => session_write_attribute
  end type writer_session_t

contains

  ! Base session procedures
  logical function is_session_open(self)
    class(io_session_base_t), intent(in) :: self
    is_session_open = self%is_open
  end function is_session_open

  logical function is_session_functional(self)
    class(io_session_base_t), intent(in) :: self
    is_session_functional = self%is_functional
  end function is_session_functional

  function get_file(self) result(file_ptr)
    class(io_session_base_t), intent(in) :: self
    class(io_file_t), pointer :: file_ptr
    file_ptr => self%file
  end function get_file

  subroutine session_base_close(self)
    class(io_session_base_t), intent(inout) :: self
    if (.not. self%is_open) return
    call self%file%close()
    self%is_open = .false.
  end subroutine session_base_close

  subroutine session_base_finalise(self)
    class(io_session_base_t), intent(inout) :: self
    if (self%is_open) call self%close()
    if (associated(self%file)) self%file => null()
  end subroutine session_base_finalise

  ! Reader session procedures
  subroutine reader_session_open(self, filename, comm)
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    if (self%is_open) error stop "IO session already open"

    call allocate_io_reader(self%reader)
    call self%reader%init(comm, "session_reader")
    self%file => self%reader%open(filename, io_mode_read, comm)
    call self%file%begin_step()
    self%is_open = .true.
  end subroutine reader_session_open

  subroutine read_data_i8(self, variable_name, value)
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_i8

  subroutine read_data_integer(self, variable_name, value)
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_integer

  subroutine read_data_real(self, variable_name, value)
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_real

  subroutine read_data_array_3d( &
    self, variable_name, array, start_dims, count_dims &
    )
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data( &
      variable_name, array, self%file, &
      start_dims=start_dims, count_dims=count_dims &
      )
  end subroutine read_data_array_3d

  ! Writer session procedures
  subroutine writer_session_open(self, filename, comm)
    use m_io_dummy, only: io_dummy_file_t
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    if (self%is_open) error stop "IO session already open"

    call allocate_io_writer(self%writer)
    call self%writer%init(comm, "session_writer")
    self%file => self%writer%open(filename, io_mode_write, comm)
    call self%file%begin_step()

    ! check if file was actually opened (dummy I/O returns is_open = .false.)
    self%is_functional = .true.
    select type (file => self%file)
    type is (io_dummy_file_t)
      self%is_functional = file%is_open
    end select

    self%is_open = .true.  ! always mark session as open so operations don't fail
  end subroutine writer_session_open

  subroutine write_data_i8(self, variable_name, value)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_i8

  subroutine write_data_integer(self, variable_name, value)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_integer

  subroutine write_data_real(self, variable_name, value)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_real

  subroutine write_data_array_3d( &
    self, variable_name, array, start_dims, count_dims &
    )
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data( &
      variable_name, array, self%file, &
      start_dims=start_dims, count_dims=count_dims &
      )
  end subroutine write_data_array_3d

  subroutine session_write_attribute(self, attribute_name, attribute_value)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: attribute_value

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_attribute( &
      attribute_name, attribute_value, self%file &
      )
  end subroutine session_write_attribute

end module m_io_session
