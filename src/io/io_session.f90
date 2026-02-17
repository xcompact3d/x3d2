module m_io_session
!! @brief Provides high-level, session-based user interface for all I/O
!! operations
!!
!! @details This module is the sole entry point for file reading and writing.
!! It abstracts away all backend details and provides a type-safe interface
!! for all I/O tasks.
!!
!! Key features:
!! - Type-safe sessions: specialised `reader_session_t` and `writer_session_t`
!!   types for reading and writing operations, respectively.
!! - Automatic backend selection: based on compile-time options
!! - Resource cleanup: memory is automatically freed when sessions
!!   go out of scope (using final subroutines).
!! - Simplified workflow - user only needs to manage a simple
!! `open -> read/write -> close` workflow, with no need for manual file handle
!!   management or explicit cleanup calls.
!!
!! @example
!! A typical usage pattern for reading data and writing data:
!!
!! @code{.f90}
!! use m_io_session, only: writer_session_t, reader_session_t
!!
!! implicit none
!!
!! real, dimension(:,:,:), allocatable :: temp_field
!! type(writer_session_t)           :: writer
!! type(reader_session_t)           :: reader
!!
!! ! For writing data
!! call writer%open("output.bp")
!! call writer%write_data("temperature", temp_field)
!! call writer%close()
!! ! Note: writer is automatically cleaned up when it goes out of scope
!!
!! ! For reading data
!! call reader%open("input.bp")
!! call reader%read_data("temperature", temp_field)
!! call reader%close()
!! ! Note: reader is automatically cleaned up when it goes out of scope
!! @endcode
!!
!! @note Users should only use the types provided by this module. The lower-level
!! modules like `m_io_base` and `m_io_backend` are internal components and should
!! never be used directly in user code.
  use m_common, only: dp, i8
  use m_field, only: field_t
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, &
                       io_mode_read, io_mode_write
  use m_io_backend, only: allocate_io_reader, allocate_io_writer

  implicit none

  private

  !! Public session types for user interaction
  public :: reader_session_t, writer_session_t

  !> Base type for common session functionality
  type :: io_session_base_t
    private
    class(io_file_t), allocatable :: file
    logical :: is_open = .false.
    logical :: is_functional = .true.  ! false for dummy I/O
  contains
    procedure :: is_session_open
    procedure :: is_session_functional
    procedure :: close => session_base_close
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
    class(io_reader_t), allocatable :: reader
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
    final :: reader_session_finaliser
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
    class(io_writer_t), allocatable :: writer
  contains
    ! Open/close operations
    procedure :: open => writer_session_open
    procedure :: begin_step => writer_session_begin_step
    procedure :: end_step => writer_session_end_step
    ! Generic write_data interface
    generic :: write_data => write_data_i8, write_data_integer, &
      write_data_real, write_data_array_3d
    procedure, private :: write_data_i8
    procedure, private :: write_data_integer
    procedure, private :: write_data_real
    procedure, private :: write_data_array_3d
    ! Field-from-solver interface (GPU-aware when available)
    procedure :: write_field_from_solver => session_write_field_from_solver
    procedure :: supports_device_field_write => &
      session_supports_device_field_write
    procedure :: sync_device => session_sync_device
    ! Write attribute interface
    procedure :: write_attribute => session_write_attribute
    final :: writer_session_finaliser
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

  subroutine session_base_close(self)
    class(io_session_base_t), intent(inout) :: self
    if (.not. self%is_open) return
    call self%file%close()
    self%is_open = .false.
  end subroutine session_base_close

  ! Reader session procedures
  subroutine reader_session_open(self, filename, comm)
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    if (self%is_open) error stop "IO session already open"

    call allocate_io_reader(self%reader)
    call self%reader%init(comm, "session_reader")
    self%file = self%reader%open(filename, io_mode_read, comm)
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
    self, variable_name, array, start_dims, count_dims, shape_dims &
    )
    class(reader_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)
    integer(i8), intent(in), optional :: shape_dims(3)

    if (.not. self%is_open) error stop "IO session not open"
    call self%reader%read_data( &
      variable_name, array, self%file, &
      start_dims=start_dims, count_dims=count_dims, shape_dims=shape_dims &
      )
  end subroutine read_data_array_3d

  ! Writer session procedures
  subroutine writer_session_open(self, filename, comm)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    if (self%is_open) error stop "IO session already open"

    call allocate_io_writer(self%writer)
    call self%writer%init(comm, "session_writer")
    self%file = self%writer%open(filename, io_mode_write, comm)
    call self%file%begin_step()

    ! check if backend is functional
    self%is_functional = self%file%is_file_functional()

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

  subroutine write_data_real(self, variable_name, value, use_sp)
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    logical, intent(in), optional :: use_sp

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data(variable_name, value, self%file, use_sp)
  end subroutine write_data_real

  subroutine write_data_array_3d( &
    self, variable_name, array, shape_dims, start_dims, count_dims, use_sp &
    )
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_data( &
      variable_name, array, self%file, &
      shape_dims, start_dims, count_dims, use_sp &
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

  subroutine session_write_field_from_solver( &
    self, variable_name, field, backend, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    !! Write field data with backend-specific optimisations
    class(writer_session_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    class(*), intent(in) :: field
    class(*), intent(in) :: backend
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp

    if (.not. self%is_open) error stop "IO session not open"
    call self%writer%write_field_from_solver( &
      variable_name, field, self%file, backend, &
      shape_dims, start_dims, count_dims, use_sp &
      )
  end subroutine session_write_field_from_solver

  logical function session_supports_device_field_write(self)
    class(writer_session_t), intent(in) :: self
    session_supports_device_field_write = &
      self%writer%supports_device_field_write()
  end function session_supports_device_field_write

  subroutine session_sync_device(self)
    !! Synchronise device before a batch of field writes.
    !! Does nothing for non-GPU backends.
    class(writer_session_t), intent(inout) :: self
    call self%writer%sync_device()
  end subroutine session_sync_device

  subroutine writer_session_begin_step(self)
    !! Begin a new timestep for writing (used for time-series in single file)
    class(writer_session_t), intent(inout) :: self

    if (.not. self%is_open) error stop "IO session not open"
    call self%file%begin_step()
  end subroutine writer_session_begin_step

  subroutine writer_session_end_step(self)
    !! End the current timestep for writing
    class(writer_session_t), intent(inout) :: self

    if (.not. self%is_open) error stop "IO session not open"
    call self%file%end_step()
  end subroutine writer_session_end_step

  !> Finalisation for reader_session_t
  !! Called automatically when a reader_session_t goes out of scope
  !! Ensures proper cleanup even if user forgets to call close
  subroutine reader_session_finaliser(self)
    type(reader_session_t) :: self
    if (self%is_open) call self%close()

    if (allocated(self%file)) then
      deallocate (self%file)
    end if

    if (allocated(self%reader)) then
      call self%reader%finalise()
      deallocate (self%reader)
    end if
  end subroutine reader_session_finaliser

  !> Finalisation for writer_session_t
  !! Called automatically when a writer_session_t goes out of scope
  !! Ensures proper cleanup even if user forgets to call close
  subroutine writer_session_finaliser(self)
    type(writer_session_t) :: self

    if (self%is_open) call self%close()

    if (allocated(self%file)) then
      deallocate (self%file)
    end if

    if (allocated(self%writer)) then
      call self%writer%finalise()
      deallocate (self%writer)
    end if
  end subroutine writer_session_finaliser

end module m_io_session
