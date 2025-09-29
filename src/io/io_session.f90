module m_io_session
!! High-level I/O session providing a unified session-based interface
!! This module abstracts away the complexity of the underlying I/O system,
!! providing the io_session_t type for efficient multi-variable I/O operations
!! and factory functions for creating the appropriate I/O backend.

  use mpi, only: MPI_COMM_WORLD
  use m_common, only: dp, i8
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, io_mode_read, io_mode_write
  use m_io_factory, only: allocate_io_reader, allocate_io_writer

  implicit none

  private
  
  ! Re-export factory functions through session module for proper layering
  public :: allocate_io_reader, allocate_io_writer
  
  ! Session-based API
  public :: io_session_t

  !> Type for reading and writing multiple variables from/to the same file efficiently
  !!
  !! Usage example for reading:
  !!   type(io_session_t) :: io_session
  !!   call io_session%open("checkpoint.bp", MPI_COMM_WORLD, io_mode_read)
  !!   call io_session%read_data("timestep", timestep)
  !!   call io_session%read_data("velocity_u", u_field, start_dims, count_dims)
  !!   call io_session%close()
  !!
  !! Usage example for writing:
  !!   type(io_session_t) :: io_session
  !!   call io_session%open("output.bp", MPI_COMM_WORLD, io_mode_write)
  !!   call io_session%write_data("timestep", current_step)
  !!   call io_session%write_data("pressure", p_field, start_dims, count_dims)
  !!   call io_session%write_attribute("ParaView", "vtk_xml_content")
  !!   call io_session%close()
  !!
  !! The session automatically handles backend selection (currently just ADIOS2 but can be extended to MPI-IO, etc.)
  !! and manages file handles, readers, and writers internally.
  type :: io_session_t
    private
    class(io_reader_t), pointer :: reader => null()
    class(io_writer_t), pointer :: writer => null()
    class(io_file_t), pointer :: file => null()
    logical :: is_open = .false.
    logical :: is_write_mode = .false.
  contains
    ! Open/close operations
    procedure :: open => session_open
    procedure :: close => session_close
    
    ! Generic read_data interface
    generic :: read_data => read_data_i8, read_data_integer, read_data_real, read_data_array_3d
    procedure, private :: read_data_i8
    procedure, private :: read_data_integer
    procedure, private :: read_data_real
    procedure, private :: read_data_array_3d
    
    ! Generic write_data interface
    generic :: write_data => write_data_i8, write_data_integer, write_data_real, write_data_array_3d
    procedure, private :: write_data_i8
    procedure, private :: write_data_integer
    procedure, private :: write_data_real
    procedure, private :: write_data_array_3d
    
    ! Write attribute interface
    procedure :: write_attribute => session_write_attribute
  end type io_session_t

contains

  subroutine read_data_i8(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_i8

  subroutine read_data_integer(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_integer

  subroutine read_data_real(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, value, self%file)
  end subroutine read_data_real

  subroutine read_data_array_3d(self, variable_name, array, start_dims, count_dims)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data(variable_name, array, self%file, start_dims=start_dims, count_dims=count_dims)
  end subroutine read_data_array_3d

  subroutine write_data_i8(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    if (.not. self%is_write_mode) error stop "IO session not in write mode"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_i8

  subroutine write_data_integer(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    if (.not. self%is_write_mode) error stop "IO session not in write mode"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_integer

  subroutine write_data_real(self, variable_name, value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value

    if (.not. self%is_open) error stop "IO session not open"
    if (.not. self%is_write_mode) error stop "IO session not in write mode"
    call self%writer%write_data(variable_name, value, self%file)
  end subroutine write_data_real

  subroutine write_data_array_3d(self, variable_name, array, start_dims, count_dims)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    if (.not. self%is_open) error stop "IO session not open"
    if (.not. self%is_write_mode) error stop "IO session not in write mode"
    call self%writer%write_data(variable_name, array, self%file, start_dims=start_dims, count_dims=count_dims)
  end subroutine write_data_array_3d

  subroutine session_write_attribute(self, attribute_name, attribute_value)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: attribute_value

    if (.not. self%is_open) error stop "IO session not open"
    if (.not. self%is_write_mode) error stop "IO session not in write mode"
    call self%writer%write_attribute(attribute_name, attribute_value, self%file)
  end subroutine session_write_attribute

  subroutine session_open(self, filename, comm, mode)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(in), optional :: mode

    integer :: use_mode

    if (self%is_open) error stop "IO session already open"
    
    ! Default to read mode if not specified
    use_mode = io_mode_read
    if (present(mode)) use_mode = mode
    
    if (use_mode == io_mode_read) then
      call allocate_io_reader(self%reader)
      call self%reader%init(comm, "session_reader")
      self%file => self%reader%open(filename, io_mode_read, comm)
      self%is_write_mode = .false.
    else if (use_mode == io_mode_write) then
      call allocate_io_writer(self%writer)
      call self%writer%init(comm, "session_writer")
      self%file => self%writer%open(filename, io_mode_write, comm)
      self%is_write_mode = .true.
    else
      error stop "Invalid I/O mode"
    end if
    
    call self%file%begin_step()
    self%is_open = .true.
  end subroutine session_open

  subroutine session_close(self)
    class(io_session_t), intent(inout) :: self

    if (.not. self%is_open) return
    call self%file%close()
    self%is_open = .false.
    self%is_write_mode = .false.
  end subroutine session_close

end module m_io_session
