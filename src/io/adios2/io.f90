module m_io_backend
!! @brief Provides ADIOS2-specific implementation of the I/O backend interface
!!
!! @details This module contains the concrete backend implementation for ADIOS2
!! (ADaptive Input Output System v2) library. It acts as a translation layer
!! converting generic I/O calls from the session interface into specific calls
!! to the ADIOS2 API.
!!
!! The `adios2_reader_t` and `adios2_writer_t` types defined here extend the
!! abstract base types from `m_io_base` and implement required procedures
!!
!! This backend leverages several key features of the underlying ADIOS2 library
!! - engine abstraction - the same API can be used for different transport
!! methods (e.g. BP4, BP5, HDF5)
!! - Asynchronous I/O - by default ADIOS2 uses a deferred transport mode
!! which can improve performance by overlapping computation and I/O
!! - MPI integration - it is designed for large-scale paralle I/O and
!! integrates with MPI, though serial operation is also supported
!!
!! @note This is an internal backend module and should never be used directly.
!! All user interaction must go through `m_io_session`.
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
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, &
                       io_mode_read, io_mode_write

  implicit none

  private
  public :: allocate_io_reader, allocate_io_writer
  public :: get_default_backend, IO_BACKEND_DUMMY, IO_BACKEND_ADIOS2

  integer, parameter :: IO_BACKEND_DUMMY = 0
  integer, parameter :: IO_BACKEND_ADIOS2 = 1

  type, extends(io_reader_t) :: io_adios2_reader_t
    private
    type(adios2_adios) :: adios              !! ADIOS2 global handler
    type(adios2_io) :: io_handle             !! ADIOS2 IO object for managing I/O
    logical :: is_step_active = .false.      !! Flag to track if a step is active
    integer :: comm = MPI_COMM_NULL          !! MPI communicator
  contains
    procedure :: init => reader_init_adios2
    procedure :: open => reader_open_adios2
    procedure :: read_data_i8 => read_data_i8_adios2
    procedure :: read_data_integer => read_data_integer_adios2
    procedure :: read_data_real => read_data_real_adios2
    procedure :: read_data_array_3d => read_data_array_3d_adios2
    procedure :: finalise => finalise_reader_adios2
    procedure, private :: handle_error => handle_error_reader
  end type io_adios2_reader_t

  type, extends(io_writer_t) :: io_adios2_writer_t
    private
    type(adios2_adios) :: adios              !! ADIOS2 global handler
    type(adios2_io) :: io_handle             !! ADIOS2 IO object for managing I/O
    logical :: is_step_active = .false.      !! Flag to track if a step is active
    integer :: comm = MPI_COMM_NULL          !! MPI communicator
  contains
    procedure :: init => writer_init_adios2
    procedure :: open => writer_open_adios2
    procedure :: write_data_i8 => write_data_i8_adios2
    procedure :: write_data_integer => write_data_integer_adios2
    procedure :: write_data_real => write_data_real_adios2
    procedure :: write_data_array_3d => write_data_array_3d_adios2
    procedure :: write_attribute_string => write_attribute_string_adios2
    procedure :: write_attribute_array_1d_real => &
      write_attribute_array_1d_real_adios2
    procedure :: finalise => finalise_writer_adios2
    procedure, private :: handle_error => handle_error_writer
  end type io_adios2_writer_t

  type, extends(io_file_t) :: io_adios2_file_t
    private
    type(adios2_engine) :: engine            !! ADIOS2 engine for data reading/writing
    logical :: is_step_active = .false.      !! Flag to track if a step is active
    logical :: is_writer = .false.           !! Flag to track if this is for writing
  contains
    procedure :: close => file_close_adios2
    procedure :: begin_step => file_begin_step_adios2
    procedure :: end_step => file_end_step_adios2
    procedure, private :: handle_error => handle_error_file
  end type io_adios2_file_t

contains

  subroutine allocate_io_reader(reader)
    class(io_reader_t), allocatable, intent(out) :: reader
    allocate (io_adios2_reader_t :: reader)
  end subroutine allocate_io_reader

  subroutine allocate_io_writer(writer)
    class(io_writer_t), allocatable, intent(out) :: writer
    allocate (io_adios2_writer_t :: writer)
  end subroutine allocate_io_writer

  function get_default_backend() result(backend)
    integer :: backend
    backend = IO_BACKEND_ADIOS2
  end function get_default_backend

  subroutine reader_init_adios2(self, comm, name)
    class(io_adios2_reader_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name

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
    call adios2_declare_io(self%io_handle, self%adios, name, ierr)
    call self%handle_error(ierr, "Failed to declare ADIOS2 IO object")

    ! hardcode engine type to BP5
    call adios2_set_engine(self%io_handle, "BP5", ierr)

    if (.not. self%io_handle%valid) &
      call self%handle_error(1, "Failed to create ADIOS2 IO object")
  end subroutine reader_init_adios2

  function reader_open_adios2(self, filename, mode, comm) result(file_handle)
    class(io_adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm

    class(io_file_t), allocatable :: file_handle
    type(io_adios2_file_t) :: temp_handle
    integer :: ierr, use_comm

    use_comm = comm
    if (.not. self%io_handle%valid) &
      call self%handle_error(1, "ADIOS2 IO object is not valid")

    call adios2_open( &
      temp_handle%engine, self%io_handle, filename, &
      adios2_mode_read, use_comm, ierr)
    call self%handle_error(ierr, "Failed to open ADIOS2 engine for reading")
    temp_handle%is_writer = .false.

    file_handle = temp_handle
  end function reader_open_adios2

  subroutine read_data_i8_adios2(self, variable_name, value, file_handle)
    class(io_adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr == adios2_found) then
        call adios2_get(file_handle%engine, var, value, adios2_mode_sync, ierr)
        call self%handle_error(ierr, "Failed to read variable "//variable_name)
      else
        call self%handle_error(1, "Variable " &
                               //trim(variable_name)//" not found in file")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine read_data_i8_adios2

  subroutine read_data_integer_adios2(self, variable_name, value, file_handle)
    class(io_adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr == adios2_found) then
        call adios2_get(file_handle%engine, var, value, adios2_mode_sync, ierr)
        call self%handle_error(ierr, "Failed to read variable "//variable_name)
      else
        call self%handle_error(1, "Variable " &
                               //trim(variable_name)//" not found in file")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine read_data_integer_adios2

  subroutine read_data_real_adios2(self, variable_name, value, file_handle)
    class(io_adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(out) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      ! retrieve a variable handler within current io handler
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr == adios2_found) then
        call adios2_get(file_handle%engine, var, value, adios2_mode_sync, ierr)
        call self%handle_error(ierr, "Failed to read variable "//variable_name)
      else
        call self%handle_error(1, "Variable " &
                               //trim(variable_name)//" not found in file")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine read_data_real_adios2

  subroutine read_data_array_3d_adios2( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims &
    )
    class(io_adios2_reader_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(inout) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    type(adios2_variable) :: var
    integer :: ierr
    integer(i8) :: local_start(3), local_count(3)

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr == adios2_found) then
        if (present(start_dims)) then
          local_start = start_dims
        else
          local_start = 0_i8
        end if
        if (present(count_dims)) then
          local_count = count_dims
        else
          local_count = int(shape(array), i8)
        end if

        ! apply selection only when explicitly requested (partial array reads)
        if (present(start_dims) .or. present(count_dims)) then
          call adios2_set_selection(var, 3, local_start, local_count, ierr)
          call self%handle_error(ierr, &
                                 "Failed to set selection for variable " &
                                 //trim(variable_name))
        end if

        call adios2_get(file_handle%engine, var, array, adios2_mode_sync, ierr)
        call self%handle_error(ierr, &
                               "Failed to read variable "//trim(variable_name))
      else
        call self%handle_error(1, "Variable " &
                               //trim(variable_name)//" not found in file")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine read_data_array_3d_adios2

  subroutine finalise_reader_adios2(self)
    class(io_adios2_reader_t), intent(inout) :: self
    integer :: ierr

    if (self%adios%valid) then
      call adios2_finalize(self%adios, ierr)
      call self%handle_error(ierr, "Failed to finalise ADIOS2")
    end if
  end subroutine finalise_reader_adios2

  subroutine writer_init_adios2(self, comm, name)
    class(io_adios2_writer_t), intent(inout) :: self
    integer, intent(in) :: comm
    character(len=*), intent(in) :: name

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
    call adios2_declare_io(self%io_handle, self%adios, name, ierr)
    call self%handle_error(ierr, "Failed to declare ADIOS2 IO object")

    ! hardcode engine type to BP5
    call adios2_set_engine(self%io_handle, "BP5", ierr)

    if (.not. self%io_handle%valid) &
      call self%handle_error(1, "Failed to create ADIOS2 IO object")
  end subroutine writer_init_adios2

  function writer_open_adios2(self, filename, mode, comm) result(file_handle)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(in) :: comm

    class(io_file_t), allocatable :: file_handle
    type(io_adios2_file_t) :: temp_handle
    integer :: ierr, use_comm

    use_comm = comm
    if (.not. self%io_handle%valid) &
      call self%handle_error(1, "ADIOS2 IO object is not valid")

    ! if opening in write mode, we are starting a new independent dataset
    ! remove all old variables from the IO object
    if (mode == io_mode_write) then
      call adios2_remove_all_variables(self%io_handle, ierr)
      call self%handle_error(ierr, "Failed to remove old ADIOS2 variables &
                             & before open")
    end if

    call adios2_open( &
      temp_handle%engine, self%io_handle, filename, &
      adios2_mode_write, use_comm, ierr)
    call self%handle_error(ierr, "Failed to open ADIOS2 engine for writing")
    temp_handle%is_writer = .true.

    file_handle = temp_handle
  end function writer_open_adios2

  subroutine write_data_i8_adios2(self, variable_name, value, file_handle)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        ! Use integer4 type for i8 compatibility - ADIOS2 will handle the conversion
        call adios2_define_variable( &
          var, self%io_handle, variable_name, adios2_type_integer4, ierr)
        call self%handle_error(ierr, &
                               "Error defining ADIOS2 scalar i8 variable")
      end if

      call adios2_put( &
        file_handle%engine, var, value, adios2_mode_deferred, ierr)
      call self%handle_error(ierr, "Error writing ADIOS2 scalar i8 data")
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_i8_adios2

  subroutine write_data_integer_adios2(self, variable_name, value, file_handle)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer, intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable( &
          var, self%io_handle, variable_name, adios2_type_integer4, ierr)
        call self%handle_error(ierr, &
                               "Error defining ADIOS2 scalar integer variable")
      end if

      call adios2_put( &
        file_handle%engine, var, value, adios2_mode_deferred, ierr)
      call self%handle_error(ierr, "Error writing ADIOS2 scalar integer data")
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_integer_adios2

  subroutine write_data_real_adios2(self, variable_name, value, file_handle)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable(var, self%io_handle, variable_name, &
                                    adios2_type_dp, ierr)
        call self%handle_error(ierr, "Error defining ADIOS2 scalar &
                                     & double precision real variable")
      end if

      call adios2_put( &
        file_handle%engine, var, value, adios2_mode_deferred, ierr)
      call self%handle_error(ierr, "Error writing ADIOS2 scalar &
                                   & double precision real data")
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_real_adios2

  subroutine write_data_array_3d_adios2( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims &
    )
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in), optional :: shape_dims(3)
    integer(i8), intent(in), optional :: start_dims(3)
    integer(i8), intent(in), optional :: count_dims(3)

    type(adios2_variable) :: var
    integer :: ierr
    integer(i8) :: local_shape(3), local_start(3), local_count(3)

    select type (file_handle)
    type is (io_adios2_file_t)
      if (present(shape_dims)) then
        local_shape = shape_dims
      else
        local_shape = int(shape(array), i8)
      end if

      if (present(start_dims)) then
        local_start = start_dims
      else
        local_start = 0_i8
      end if

      if (present(count_dims)) then
        local_count = count_dims
      else
        local_count = int(shape(array), i8)
      end if

      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable(var, self%io_handle, variable_name, &
                                    adios2_type_dp, 3, local_shape, &
                                    local_start, local_count, &
                                    adios2_constant_dims, ierr)
        call self%handle_error(ierr, "Error defining ADIOS2 3D array &
                                     & double precision real variable")
      end if

      call adios2_put( &
        file_handle%engine, var, array, adios2_mode_deferred, ierr)
      call self%handle_error(ierr, "Error writing ADIOS2 3D array &
                                   & double precision real data")
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_array_3d_adios2

  subroutine write_attribute_string_adios2( &
    self, attribute_name, value, file_handle &
    )
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    character(len=*), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_attribute) :: attr
    integer :: ierr

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_define_attribute( &
        attr, self%io_handle, attribute_name, value, ierr)
      call self%handle_error(ierr, &
                             "Error defining ADIOS2 string attribute " &
                             //trim(attribute_name))
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_attribute_string_adios2

  subroutine write_attribute_array_1d_real_adios2( &
    self, attribute_name, values, file_handle &
    )
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: attribute_name
    real(dp), intent(in) :: values(:)
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_attribute) :: attr
    integer :: ierr
    integer :: num_elements

    select type (file_handle)
    type is (io_adios2_file_t)
      num_elements = size(values)
      call adios2_define_attribute( &
        attr, self%io_handle, attribute_name, values, num_elements, ierr)
      call self%handle_error( &
                             ierr, "Error defining ADIOS2 real &
                             & array attribute " &
                             //trim(attribute_name))
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_attribute_array_1d_real_adios2

  subroutine finalise_writer_adios2(self)
    class(io_adios2_writer_t), intent(inout) :: self
    integer :: ierr

    if (self%adios%valid) then
      call adios2_finalize(self%adios, ierr)
      call self%handle_error(ierr, "Failed to finalise ADIOS2")
    end if
  end subroutine finalise_writer_adios2

  subroutine file_close_adios2(self)
    class(io_adios2_file_t), intent(inout) :: self

    integer :: ierr

    if (self%is_step_active) call self%end_step()

    if (self%engine%valid) then
      call adios2_close(self%engine, ierr)
      call self%handle_error(ierr, "Failed to close ADIOS2 engine")
    end if
  end subroutine file_close_adios2

  subroutine file_begin_step_adios2(self)
    class(io_adios2_file_t), intent(inout) :: self

    integer :: ierr

    if (self%is_step_active) return

    if (self%is_writer) then
      call adios2_begin_step(self%engine, adios2_step_mode_append, ierr)
      call self%handle_error(ierr, "Error beginning ADIOS2 step for writing")
    else
      call adios2_begin_step(self%engine, adios2_step_mode_read, ierr)
      call self%handle_error(ierr, "Error beginning ADIOS2 step for reading")
    end if

    self%is_step_active = .true.
  end subroutine file_begin_step_adios2

  subroutine file_end_step_adios2(self)
    class(io_adios2_file_t), intent(inout) :: self

    integer :: ierr

    if (.not. self%is_step_active) return
    call adios2_end_step(self%engine, ierr)
    call self%handle_error(ierr, "Failed to end ADIOS2 step")
    self%is_step_active = .false.
  end subroutine file_end_step_adios2

  subroutine handle_error_reader(self, ierr, message)
    class(io_adios2_reader_t), intent(inout) :: self
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message

    if (ierr /= 0) then
      print *, "ADIOS2 Reader Error: ", message
      print *, "Error code: ", ierr
      error stop
    end if
  end subroutine handle_error_reader

  subroutine handle_error_writer(self, ierr, message)
    class(io_adios2_writer_t), intent(inout) :: self
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message

    if (ierr /= 0) then
      print *, "ADIOS2 Writer Error: ", message
      print *, "Error code: ", ierr
      error stop
    end if
  end subroutine handle_error_writer

  subroutine handle_error_file(self, ierr, message)
    class(io_adios2_file_t), intent(inout) :: self
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message

    if (ierr /= 0) then
      print *, "ADIOS2 File Error: ", message
      print *, "Error code: ", ierr
      error stop
    end if
  end subroutine handle_error_file

end module m_io_backend
