module m_adios2_io
!! This module contains ADIOS2 (ADaptable Input Output System version 2)
!! operations for reading and writing data. ADIOS2 transports data as
!! groups of self-describing variables and attributes across different
!! media types (e.g., files, memory, network).
!! ADIOS2 APIs are based on:
!! - MPI (although non-MPI serial code is also supported)
!! - Deferred/prefetch/grouped variables transport mode by default
!! - Engine abstraction for reusing the APIs for different transport modes
  use adios2, only: adios2_adios, adios2_io, adios2_engine, &
                    adios2_variable, adios2_attribute, &
                    adios2_mode_sync, adios2_mode_write, &
                    adios2_mode_deferred, adios2_mode_read, &
                    adios2_step_mode_append, adios2_step_mode_read, &
                    adios2_init, adios2_finalize, &
                    adios2_declare_io, adios2_set_engine, &
                    adios2_open, adios2_close, &
                    adios2_begin_step, adios2_end_step, &
                    adios2_define_variable, adios2_inquire_variable, &
                    adios2_define_attribute, &
                    adios2_set_selection, adios2_put, &
                    adios2_get, adios2_remove_all_variables, &
                    adios2_found, adios2_constant_dims, &
                    adios2_type_dp, adios2_type_integer4
  use mpi, only: MPI_COMM_NULL, MPI_Initialized, MPI_Comm_rank
  use m_common, only: dp, i8
  implicit none

  private
  public :: adios2_writer_t, adios2_reader_t, adios2_file_t, &
            adios2_mode_write, adios2_mode_read

  !> ADIOS2 base type
  !> Abstract base module for ADIOS2 operations
  !> Defines the abstract base type `base_adios2_t` for common ADIOS2 components
  type, abstract :: base_adios2_t
    private
    type(adios2_adios) :: adios              !! ADIOS2 global handler
    type(adios2_io) :: io                    !! ADIOS2 IO object for managing I/O
    type(adios2_engine) :: engine            !! ADIOS2 engine for data reading/writing
    logical :: is_step_active = .false.      !! Flag to track if a step is active
    integer :: comm = MPI_COMM_NULL          !! MPI communicator

  contains
    procedure, public :: init                !! Initialises ADIOS2 handler
    procedure, public :: open                !! Opens an ADIOS2 engine
    procedure, public :: close               !! Closes the ADIOS2 session
    procedure, public :: end_step            !! Ends a step in the ADIOS2 engine
    procedure, public :: handle_error        !! Error handling for ADIOS2 operations
    procedure, public :: finalise            !! Finalises ADIOS2 handler

    procedure(begin_step), deferred, public :: begin_step !! Begins a step in the ADIOS2 engine
  end type base_adios2_t

  !> ADIOS2 writer type
  type, extends(base_adios2_t) :: adios2_writer_t
  contains
    procedure, public :: begin_step => begin_step_writer

    generic, public :: write_data => write_scalar_int, &
      write_scalar_real, &
      write_array_1d_real, &
      write_array_2d_real, &
      write_array_3d_real, &
      write_array_1d_int, &
      write_array_4d_real
    generic, public :: write_attribute => write_attribute_string, &
      write_attribute_array_1d_real

    procedure, private :: write_scalar_int
    procedure, private :: write_scalar_real
    procedure, private :: write_array_1d_int
    procedure, private :: write_array_1d_real
    procedure, private :: write_array_2d_real
    procedure, private :: write_array_3d_real
    procedure, private :: write_array_4d_real
    procedure, private :: write_attribute_string
    procedure, private :: write_attribute_array_1d_real
  end type adios2_writer_t

  !> ADIOS2 reader type
  type, extends(base_adios2_t) :: adios2_reader_t
  contains
    procedure, public :: begin_step => begin_step_reader

    generic, public :: read_data => read_scalar_integer, &
      read_scalar_i8, &
      read_scalar_real, &
      read_array_2d_real, &
      read_array_3d_real

    procedure, private :: read_scalar_integer
    procedure, private :: read_scalar_i8
    procedure, private :: read_scalar_real
    procedure, private :: read_array_2d_real
    procedure, private :: read_array_3d_real
  end type adios2_reader_t

  !> ADIOS2 file type
  type :: adios2_file_t
    type(adios2_engine) :: engine
  end type adios2_file_t

  abstract interface
    !> Begins a step in the ADIOS2 engine
    subroutine begin_step(self, file)
      import :: base_adios2_t, adios2_file_t
      class(base_adios2_t), intent(inout) :: self
      type(adios2_file_t), intent(inout) :: file
    end subroutine begin_step
  end interface
contains

  !> Initialises ADIOS2
  !> self: Instance of `base_adios2_t`
  !> comm: MPI communicator (use `MPI_COMM_WORLD` for parallel runs)
  !> io_name: unique name associated with IO component inside ADIOS2
  subroutine init(self, comm, io_name)
    class(base_adios2_t), intent(inout) :: self
    integer, intent(in) :: comm
    !> io that spawns an engine based on its configuration
    character(len=*), intent(in) :: io_name

    logical :: is_mpi_initialised
    integer :: ierr, comm_rank

    if (comm == MPI_COMM_NULL) &
      call self%handle_error(1, "Invalid MPI communicator")
    call MPI_Initialized(is_mpi_initialised, ierr)
    if (.not. is_mpi_initialised) &
       call self%handle_error(1, "MPI must be initialised &
                              & before calling ADIOS2 init")

    self%comm = comm
    call MPI_Comm_rank(self%comm, comm_rank, ierr)
    call self%handle_error(ierr, "Failed to get MPI rank")

    ! create adios handler passing communicator
    call adios2_init(self%adios, comm, ierr)
    call self%handle_error(ierr, "Failed to initialise ADIOS2")

    ! declare IO process configuration inside adios
    call adios2_declare_io(self%io, self%adios, io_name, ierr)
    call self%handle_error(ierr, "Failed to declare ADIOS2 IO object")

    ! hardcode engine type to BP5
    call adios2_set_engine(self%io, "BP5", ierr)

    if (.not. self%io%valid) &
      call self%handle_error(1, "Failed to create ADIOS2 IO object")
  end subroutine init

  !> Opens an ADIOS2 engine
  !> filename: Unique engine identifier within io
  !> mode: Opening mode (write, append, read)
  function open (self, filename, mode, comm) result(file)
    class(base_adios2_t), intent(inout) :: self
    character(len=*), intent(in) :: filename   !! Unique engine identifier within io
    integer, intent(in) :: mode                !! Opening mode (write, append, read)
    integer, intent(in), optional :: comm      !! MPI communicator (optional)
    integer :: ierr, use_comm
    type(adios2_file_t) :: file   !! ADIOS2 file object

    ! if opening in write mode, we are starting a new independent dataset
    ! remove all old variables from the IO object
    if (mode == adios2_mode_write) then
      call adios2_remove_all_variables(self%io, ierr)
      call self%handle_error(ierr, "Failed to remove old ADIOS2 variables &
                             & before open")
    end if

    use_comm = self%comm
    if (present(comm)) use_comm = comm
    if (.not. self%io%valid) &
      call self%handle_error(1, "ADIOS2 IO object is not valid")

    call adios2_open(file%engine, self%io, filename, mode, use_comm, ierr)
    call self%handle_error(ierr, "Failed to open ADIOS2 engine")
  end function open

  !> Closes ADIOS2 session
  subroutine close (self, file)
    class(base_adios2_t), intent(inout) :: self
    type(adios2_file_t), intent(inout) :: file
    integer :: ierr

    if (self%is_step_active) call self%end_step(file)

    if (file%engine%valid) then
      call adios2_close(file%engine, ierr)
      call self%handle_error(ierr, "Failed to close ADIOS2 engine")
    end if
  end subroutine close

  !> Ends a step in the ADIOS2 engine
  subroutine end_step(self, file)
    class(base_adios2_t), intent(inout) :: self
    type(adios2_file_t), intent(inout) :: file
    integer :: ierr

    if (.not. self%is_step_active) return
    call adios2_end_step(file%engine, ierr)
    call self%handle_error(ierr, "Failed to end ADIOS2 step")
    self%is_step_active = .false.
  end subroutine end_step

  !> Finalises ADIOS2 handler
  subroutine finalise(self)
    class(base_adios2_t), intent(inout) :: self
    integer :: ierr

    if (self%adios%valid) then
      call adios2_finalize(self%adios, ierr)
      call self%handle_error(ierr, "Failed to finalise ADIOS2")
    end if
  end subroutine finalise

  !> Handles ADIOS2 errors
  subroutine handle_error(self, ierr, message)
    class(base_adios2_t), intent(inout) :: self
    integer, intent(in) :: ierr              !! Error code from ADIOS2 operations
    character(len=*), intent(in) :: message  !! Error message to display

    if (ierr /= 0) then
      print *, "ADIOS2 Error: ", message
      print *, "Error code: ", ierr
      error stop
    end if
  end subroutine handle_error

  !> Begin a step for ADIOS2 writer type
  subroutine begin_step_writer(self, file)
    class(adios2_writer_t), intent(inout) :: self
    type(adios2_file_t), intent(inout) :: file
    integer :: ierr

    if (self%is_step_active) return

    call adios2_begin_step(file%engine, adios2_step_mode_append, ierr)
    call self%handle_error(ierr, "Error beginning  ADIOS2 step")
    self%is_step_active = .true.
  end subroutine begin_step_writer

  !> Write scalar integer data
  subroutine write_scalar_int(self, name, data, file)
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    integer, intent(in) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, adios2_type_integer4, &
                                  ierr)
      call self%handle_error(ierr, &
                             "Error defining ADIOS2 scalar integer variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 scalar integer data")
  end subroutine write_scalar_int

  !> Write scalar real data
  subroutine write_scalar_real(self, name, data, file)
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, &
                                  adios2_type_dp, ierr)
      call self%handle_error(ierr, "Error defining ADIOS2 scalar &
                                   & double precision real variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 scalar &
                                 & double precision real data")
  end subroutine write_scalar_real

  !> Write 1d array integer data
  subroutine write_array_1d_int( &
    self, name, data, file, shape_dims, start_dims, count_dims &
    )
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name                              !! unique variable identifier within io
    integer, dimension(:), intent(in) :: data                         !! scalar real data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(:), intent(in), optional :: shape_dims, &  !! Global shape
                                                       start_dims, &  !! Local offset
                                                       count_dims     !! Local size
    type(adios2_variable) :: var                                      !! handler to newly defined variable
    integer :: ierr
    integer(i8), dimension(1) :: local_shape, local_start, local_count

    if (present(shape_dims)) then
      local_shape = shape_dims
    else
      local_shape = int(size(data), i8)
    end if

    if (present(start_dims)) then
      local_start = start_dims
    else
      local_start = 0_i8
    end if

    if (present(count_dims)) then
      local_count = count_dims
    else
      local_count = int(size(data), i8)
    end if

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      ! define adios2 variable to be written in given file format
      call adios2_define_variable(var, self%io, name, adios2_type_integer4, &
                                  1, local_shape, local_start, &
                                  local_count, adios2_constant_dims, ierr)
      call self%handle_error(ierr, &
                             "Error defining ADIOS2 1D array integer variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 1D array integer data")
  end subroutine write_array_1d_int

  !> Write 1d array real data
  subroutine write_array_1d_real( &
    self, name, data, file, shape_dims, start_dims, count_dims &
    )
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(1), intent(in), optional :: shape_dims, &
                                                       start_dims, &
                                                       count_dims
    type(adios2_variable) :: var
    integer :: ierr
    integer(i8), dimension(1) :: local_shape, local_start, local_count

    if (present(shape_dims)) then
      local_shape = shape_dims
    else
      local_shape = int(size(data), i8)
    end if

    if (present(start_dims)) then
      local_start = start_dims
    else
      local_start = 0_i8
    end if

    if (present(count_dims)) then
      local_count = count_dims
    else
      local_count = int(size(data), i8)
    end if

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, adios2_type_dp, &
                                  1, local_shape, local_start, &
                                  local_count, adios2_constant_dims, ierr)
      call self%handle_error(ierr, "Error defining ADIOS2 1D array &
                                   & double precision real variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 1D array &
                                & double precision real data")
  end subroutine write_array_1d_real

  !> Write 2d array real data
  subroutine write_array_2d_real( &
    self, name, data, file, shape_dims, start_dims, count_dims &
    )
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:, :), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(2), intent(in) :: shape_dims, &
                                             start_dims, &
                                             count_dims
    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, adios2_type_dp, &
                                  2, shape_dims, start_dims, &
                                  count_dims, adios2_constant_dims, ierr)
      call self%handle_error(ierr, "Error defining ADIOS2 2D array &
                                   & double precision real variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 2D array &
                                 & double precision real data")
  end subroutine write_array_2d_real

  !> Write 3d array real data
  subroutine write_array_3d_real( &
    self, name, data, file, shape_dims, start_dims, count_dims &
    )
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:, :, :), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(3), intent(in) :: shape_dims, &
                                             start_dims, &
                                             count_dims
    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, adios2_type_dp, &
                                  3, shape_dims, start_dims, &
                                  count_dims, adios2_constant_dims, ierr)
      call self%handle_error(ierr, "Error defining ADIOS2 3D array &
                                   & double precision real variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 3D array &
                                 & double precision real data")
  end subroutine write_array_3d_real

  !> Write 4d array real data
  subroutine write_array_4d_real( &
    self, name, data, file, shape_dims, start_dims, count_dims &
    )
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:, :, :, :), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(4), intent(in) :: shape_dims, &
                                             start_dims, &
                                             count_dims
    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)

    if (ierr /= adios2_found) then
      call adios2_define_variable(var, self%io, name, adios2_type_dp, &
                                  4, shape_dims, start_dims, &
                                  count_dims, adios2_constant_dims, ierr)
      call self%handle_error(ierr, "Error defining ADIOS2 4D array &
                                   & double precision real variable")
    end if

    call adios2_put(file%engine, var, data, adios2_mode_deferred, ierr)
    call self%handle_error(ierr, "Error writing ADIOS2 4D array &
                                 & double precision real data")
  end subroutine write_array_4d_real

  !> Write string attribute for Paraview
  subroutine write_attribute_string(self, name, value, file)
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name, value
    type(adios2_file_t), intent(inout) :: file

    type(adios2_attribute) :: attr
    integer :: ierr

    call adios2_define_attribute(attr, self%io, name, value, ierr)
    call self%handle_error(ierr, "Error defining ADIOS2 attribute")
  end subroutine write_attribute_string

  !> Write 1D array attribute for Paraview
  subroutine write_attribute_array_1d_real(self, name, data, file)
    class(adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:), intent(in) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_attribute) :: attr
    integer :: ierr
    integer :: num_elements

    num_elements = size(data)

    call adios2_define_attribute(attr, self%io, name, data, num_elements, ierr)
    call self%handle_error(ierr, &
                           "Error defining ADIOS2 1D real array attribute")
  end subroutine write_attribute_array_1d_real

  !> Begin a step for ADIOS2 reader type
  subroutine begin_step_reader(self, file)
    class(adios2_reader_t), intent(inout) :: self
    type(adios2_file_t), intent(inout) :: file

    integer :: ierr

    if (self%is_step_active) return

    call adios2_begin_step(file%engine, adios2_step_mode_read, ierr)
    call self%handle_error(ierr, "Error beginning  ADIOS2 step")
    self%is_step_active = .true.  ! set flag after successful step begin
  end subroutine begin_step_reader

  !> Read scalar integer data
  subroutine read_scalar_integer(self, name, data, file)
    class(adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    integer, intent(out) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)
    call self%handle_error(ierr, "Failed to inquire variable"//name)

    if (ierr == adios2_found) then
      call adios2_get(file%engine, var, data, adios2_mode_sync, ierr)
      call self%handle_error(ierr, "Failed to read variable"//name)
    end if
  end subroutine read_scalar_integer

  !> Read scalar long integer data
  subroutine read_scalar_i8(self, name, data, file)
    class(adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    integer(i8), intent(out) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)
    call self%handle_error(ierr, "Failed to inquire variable"//name)

    if (ierr == adios2_found) then
      call adios2_get(file%engine, var, data, adios2_mode_sync, ierr)
      call self%handle_error(ierr, "Failed to read variable"//name)
    end if
  end subroutine read_scalar_i8

  !> Read scalar real data
  subroutine read_scalar_real(self, name, data, file)
    class(adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: data
    type(adios2_file_t), intent(inout) :: file

    type(adios2_variable) :: var
    integer :: ierr

    ! retrieve a variable hanlder within current io handler
    call adios2_inquire_variable(var, self%io, name, ierr)
    call self%handle_error(ierr, "Failed to inquire variable"//name)

    if (ierr == adios2_found) then
      call adios2_get(file%engine, var, data, adios2_mode_sync, ierr)
      call self%handle_error(ierr, "Failed to read variable"//name)
    end if
  end subroutine read_scalar_real

  !> Read 2d array real data
  subroutine read_array_2d_real(self, name, data, file, start_dims, count_dims)
    class(adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:, :), allocatable, intent(out) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(2), intent(in) :: start_dims, count_dims

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)
    call self%handle_error(ierr, "Failed to inquire variable"//name)

    if (ierr == adios2_found) then
      if (.not. allocated(data)) allocate (data(count_dims(1), count_dims(2)))

      call adios2_set_selection(var, 2, start_dims, count_dims, ierr)
      call self%handle_error(ierr, &
                             "Failed to set selection for variable"//name)

      call adios2_get(file%engine, var, data, adios2_mode_sync, ierr)
      call self%handle_error(ierr, "Failed to read variable"//name)
    end if
  end subroutine read_array_2d_real

  !> Read 3d array real data
  subroutine read_array_3d_real(self, name, data, file, start_dims, count_dims)
    class(adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), dimension(:, :, :), allocatable, intent(out) :: data
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(3), intent(in) :: start_dims, count_dims

    type(adios2_variable) :: var
    integer :: ierr

    call adios2_inquire_variable(var, self%io, name, ierr)
    call self%handle_error(ierr, "Failed to inquire variable"//name)

    if (ierr == adios2_found) then
      if (.not. allocated(data)) &
        allocate (data(count_dims(1), count_dims(2), count_dims(3)))

      call adios2_set_selection(var, 3, start_dims, count_dims, ierr)
      call self%handle_error(ierr, &
                             "Failed to set selection for variable"//name)

      call adios2_get(file%engine, var, data, adios2_mode_sync, ierr)
      call self%handle_error(ierr, "Failed to read variable"//name)
    end if
  end subroutine read_array_3d_real

end module m_adios2_io
