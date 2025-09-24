module m_io_service
!! High-level I/O service providing simple, user-friendly functions
!! This module abstracts away the complexity of the underlying I/O system,
!! providing convenient functions for common I/O operations and factory functions
!! for creating the appropriate I/O backend.

  use mpi, only: MPI_COMM_WORLD
  use m_common, only: dp, i8
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, io_mode_read
#ifdef WITH_ADIOS2
  use m_io_adios2, only: io_adios2_reader_t, io_adios2_writer_t
#endif
  use m_io_dummy, only: io_dummy_reader_t, io_dummy_writer_t

  implicit none

  private
  
  ! Factory functions for creating I/O objects
  public :: create_io_reader, create_io_writer
  
  ! High-level convenience functions
  public :: read_scalar, read_array, read_multiple

  ! Session-based API
  public :: io_session_t

#ifdef WITH_ADIOS2
  logical, parameter :: use_adios2_backend = .true.
#else
  logical, parameter :: use_adios2_backend = .false.
#endif

  !> Type for reading multiple variables from the same file efficiently
  type :: io_session_t
    private
    class(io_reader_t), allocatable :: reader
    class(io_file_t), allocatable :: file
    logical :: is_open = .false.
  contains
    procedure :: open => session_open
    procedure :: close => session_close
    
    ! Generic read_data interface
    generic :: read_data => read_data_i8, read_data_real, read_data_array_3d
    procedure, private :: read_data_i8
    procedure, private :: read_data_real
    procedure, private :: read_data_array_3d
  end type io_session_t

contains

  function create_io_reader() result(reader)
    class(io_reader_t), allocatable :: reader

    if (use_adios2_backend) then
      allocate(io_adios2_reader_t :: reader)
    else
      allocate(io_dummy_reader_t :: reader)
    end if
  end function create_io_reader

  function create_io_writer() result(writer)
    class(io_writer_t), allocatable :: writer

    if (use_adios2_backend) then
      allocate(io_adios2_writer_t :: writer)
    else
      allocate(io_dummy_writer_t :: writer)
    end if
  end function create_io_writer

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

  subroutine session_open(self, filename, comm)
    class(io_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    if (self%is_open) error stop "IO session already open"
    
    allocate(self%reader, source=create_io_reader())
    call self%reader%init(comm, "session_reader")
    allocate(self%file, source=self%reader%open(filename, io_mode_read, comm))
    call self%file%begin_step()
    self%is_open = .true.
  end subroutine session_open

  subroutine session_close(self)
    class(io_session_t), intent(inout) :: self

    if (.not. self%is_open) return
    call self%file%close()
    self%is_open = .false.
  end subroutine session_close

end module m_io_service
