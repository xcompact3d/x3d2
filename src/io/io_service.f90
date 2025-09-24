module m_io_service
!! High-level I/O service providing simple, user-friendly functions
!! This module abstracts away the complexity of the underlying I/O system,
!! providing convenient functions for common I/O operations and factory functions
!! for creating the appropriate I/O backend.

  use mpi, only: MPI_COMM_WORLD
  use m_common, only: dp, i8
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, io_mode_read, io_mode_write
  use m_io_factory, only: make_reader, make_writer, make_reader_ptr, make_writer_ptr, &
                          allocate_io_reader_ptr, allocate_io_writer_ptr

  implicit none

  private
  
  ! Factory functions for creating I/O objects
  public :: create_io_reader, create_io_writer, allocate_io_writer_ptr, allocate_io_reader_ptr
  
  ! High-level convenience functions
  public :: read_scalar, read_array, read_multiple

  ! Session-based API
  public :: io_session_t

  !> Type for reading and writing multiple variables from/to the same file efficiently
  type :: io_session_t
    private
    class(io_reader_t), pointer :: reader => null()
    class(io_writer_t), pointer :: writer => null()
    class(io_file_t), pointer :: file => null()
    logical :: is_open = .false.
    logical :: is_write_mode = .false.
  contains
    ! Open/close operations
    generic :: open => session_open
    procedure, private :: session_open
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
  end type io_session_t

contains

  function create_io_reader() result(reader)
    class(io_reader_t), allocatable :: reader
    reader = make_reader()
  end function create_io_reader

  function create_io_writer() result(writer)
    class(io_writer_t), allocatable :: writer
    writer = make_writer()
  end function create_io_writer

  !> Allocate an I/O writer pointer directly without copying (avoids ADIOS2 C++ object corruption)
  subroutine allocate_io_writer_ptr(writer)
    class(io_writer_t), pointer, intent(out) :: writer
    call make_writer_ptr(writer)
  end subroutine allocate_io_writer_ptr

  subroutine allocate_io_reader_ptr(reader)
    class(io_reader_t), pointer, intent(out) :: reader
    call make_reader_ptr(reader)
  end subroutine allocate_io_reader_ptr

  !> Read a scalar integer(i8) value from a file
  subroutine read_scalar(filename, variable_name, value, comm)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value
    integer, intent(in), optional :: comm

    class(io_reader_t), allocatable :: reader
    class(io_file_t), allocatable :: file
    integer :: use_comm

    use_comm = MPI_COMM_WORLD
    if (present(comm)) use_comm = comm

    reader = create_io_reader()
    call reader%init(use_comm, "io_service_reader")
    
    file = reader%open(filename, io_mode_read, use_comm)
    call file%begin_step()
    call reader%read_data(variable_name, value, file)
    
    call file%close()
    call reader%finalise()
  end subroutine read_scalar

  !> Read a 3D array from a file
  subroutine read_array(filename, variable_name, array, start_dims, count_dims, comm)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)
    integer, intent(in), optional :: comm

    class(io_reader_t), allocatable :: reader
    class(io_file_t), allocatable :: file
    integer :: use_comm

    use_comm = MPI_COMM_WORLD
    if (present(comm)) use_comm = comm

    reader = create_io_reader()
    call reader%init(use_comm, "io_service_reader")
    
    file = reader%open(filename, io_mode_read, use_comm)
    call file%begin_step()
    call reader%read_data(variable_name, array, file, start_dims=start_dims, count_dims=count_dims)
    
    call file%close()
    call reader%finalise()
  end subroutine read_array

  !> Read multiple variables from the same file efficiently
  subroutine read_multiple(filename, comm, read_operations)
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: comm
    
    interface
      subroutine read_operations(io_session)
        import :: io_session_t
        type(io_session_t), intent(inout) :: io_session
      end subroutine read_operations
    end interface

    type(io_session_t) :: session
    integer :: use_comm

    use_comm = MPI_COMM_WORLD
    if (present(comm)) use_comm = comm

    call session%open(filename, use_comm)
    call read_operations(session)
    call session%close()
  end subroutine read_multiple

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
      call allocate_io_reader_ptr(self%reader)
      call self%reader%init(comm, "session_reader")
      self%file => self%reader%open(filename, io_mode_read, comm)
      self%is_write_mode = .false.
    else if (use_mode == io_mode_write) then
      call allocate_io_writer_ptr(self%writer)
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

end module m_io_service
