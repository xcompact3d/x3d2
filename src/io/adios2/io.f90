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
                    adios2_begin_step, adios2_end_step, adios2_perform_puts, &
                    adios2_define_variable, adios2_inquire_variable, &
                    adios2_variable_type, adios2_define_attribute, &
                    adios2_set_selection, adios2_put, &
                    adios2_get, adios2_remove_all_variables, &
                    adios2_found, adios2_constant_dims, &
                    adios2_type_dp, adios2_type_integer4, adios2_type_real
#ifdef X3D2_ADIOS2_CUDA
  use adios2, only: adios2_set_memory_space, adios2_memory_space_gpu
  use cudafor, only: cudaDeviceSynchronize, c_devloc, c_devptr
#endif
  use mpi, only: MPI_COMM_NULL, MPI_Initialized, MPI_Comm_rank, &
                 MPI_Allreduce, MPI_SUM, MPI_MAX, &
                 MPI_DOUBLE_PRECISION, MPI_Wtime
  use m_common, only: dp, i8, sp, is_sp
  use m_io_base, only: io_reader_t, io_writer_t, io_file_t, &
                       io_mode_read, io_mode_write
  use iso_c_binding, only: c_double, c_float, c_int, c_char, c_null_char, &
                           c_ptr, c_loc

  implicit none

  private
  public :: allocate_io_reader, allocate_io_writer
  public :: get_default_backend, IO_BACKEND_DUMMY, IO_BACKEND_ADIOS2

  integer, parameter :: IO_BACKEND_DUMMY = 0
  integer, parameter :: IO_BACKEND_ADIOS2 = 1

  integer, parameter :: gpu_write_mode_auto = 0
  integer, parameter :: gpu_write_mode_force_device = 1
  integer, parameter :: gpu_write_mode_force_host = 2
  real(c_double), parameter :: bytes_real_dp = &
    real(storage_size(0.0_dp)/8, c_double)
  real(c_double), parameter :: bytes_real_sp = &
    real(storage_size(0.0_sp)/8, c_double)
  real(c_double), parameter :: bytes_integer_default = &
    real(storage_size(0)/8, c_double)
  real(c_double), parameter :: bytes_integer_i8 = &
    real(storage_size(0_i8)/8, c_double)

  logical, save :: runtime_options_initialised = .false.
  logical, save :: runtime_options_reported = .false.
  logical, save :: runtime_bench_enabled = .false.
  logical, save :: runtime_bench_verbose = .true.
  logical, save :: runtime_nvtx_enabled = .true.
  integer, save :: runtime_bench_warmup_steps = 2
  integer, save :: runtime_gpu_write_mode = gpu_write_mode_auto
  character(len=16), save :: runtime_gpu_write_mode_name = "auto"

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
    procedure :: write_field_from_solver => write_field_from_solver_adios2
    procedure :: supports_device_field_write => &
      supports_device_field_write_adios2
#ifdef X3D2_ADIOS2_CUDA
    procedure :: write_data_array_3d_device => &
      write_data_array_3d_device_adios2
    procedure :: sync_device => sync_device_adios2
#endif
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
    integer :: comm = MPI_COMM_NULL
    character(len=256) :: bench_file_name = ""
    logical :: bench_enabled = .false.
    logical :: bench_verbose = .true.
    integer :: bench_warmup_steps = 0
    integer :: bench_step_counter = 0
    integer :: bench_measured_steps = 0
    real(c_double) :: bench_step_put_time_local = 0.0_c_double
    real(c_double) :: bench_step_put_bytes_local = 0.0_c_double
    real(c_double) :: bench_sum_put_time = 0.0_c_double
    real(c_double) :: bench_sum_put_time2 = 0.0_c_double
    real(c_double) :: bench_sum_end_step_time = 0.0_c_double
    real(c_double) :: bench_sum_end_step_time2 = 0.0_c_double
    real(c_double) :: bench_sum_total_time = 0.0_c_double
    real(c_double) :: bench_sum_total_time2 = 0.0_c_double
    real(c_double) :: bench_sum_throughput_gib_s = 0.0_c_double
    real(c_double) :: bench_sum_throughput_gib_s2 = 0.0_c_double
    real(c_double) :: bench_sum_bytes = 0.0_c_double
    real(c_double) :: bench_close_time = 0.0_c_double
    logical :: nvtx_step_range_active = .false.
  contains
    procedure :: close => file_close_adios2
    procedure :: begin_step => file_begin_step_adios2
    procedure :: end_step => file_end_step_adios2
    procedure, private :: handle_error => handle_error_file
  end type io_adios2_file_t

#ifdef X3D2_ADIOS2_CUDA
  ! C wrapper for GPU-aware ADIOS2 put.
  ! The engine and variable f2c handles are integer(8) in ADIOS2's
  ! Fortran bindings. The C wrapper casts them to the opaque C types.
  ! Memory space must be set to GPU via the native Fortran
  ! adios2_set_memory_space before calling this.
  interface
    subroutine adios2_put_gpu(engine_f2c, variable_f2c, data, mode, ierr) &
      bind(C, name='adios2_put_gpu')
      use iso_c_binding, only: c_ptr, c_int, c_int64_t
      integer(c_int64_t), intent(in) :: engine_f2c
      integer(c_int64_t), intent(in) :: variable_f2c
      type(c_ptr), value :: data
      integer(c_int), intent(in) :: mode
      integer(c_int), intent(out) :: ierr
    end subroutine adios2_put_gpu

    function nvtx_range_push_a(name) bind(C, name='nvtxRangePushA') &
      result(status)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value :: name
      integer(c_int) :: status
    end function nvtx_range_push_a

    function nvtx_range_pop() bind(C, name='nvtxRangePop') result(status)
      use iso_c_binding, only: c_int
      integer(c_int) :: status
    end function nvtx_range_pop
  end interface
#endif

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

  function get_adios2_vartype(use_sp) result(vartype)
    logical, intent(in) :: use_sp    !! flag for single precision output
    integer :: vartype

    if (use_sp) then
      vartype = adios2_type_real
    else if (is_sp) then
      vartype = adios2_type_real
    else
      vartype = adios2_type_dp
    end if
  end function get_adios2_vartype

  pure function to_lower_ascii(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i, code

    lowered = text
    do i = 1, len(text)
      code = iachar(lowered(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lowered(i:i) = achar(code + 32)
      end if
    end do
  end function to_lower_ascii

  logical function env_to_logical(name, default_value) result(value)
    character(len=*), intent(in) :: name
    logical, intent(in) :: default_value

    character(len=64) :: raw_value
    character(len=64) :: parsed_value
    integer :: status, value_length

    value = default_value

    call get_environment_variable(name, raw_value, length=value_length, &
                                  status=status)
    if (status /= 0 .or. value_length <= 0) return

    parsed_value = to_lower_ascii(adjustl(raw_value(1:value_length)))
    select case (trim(parsed_value))
    case ("1", "true", "yes", "on")
      value = .true.
    case ("0", "false", "no", "off")
      value = .false.
    end select
  end function env_to_logical

  integer function env_to_integer(name, default_value, min_value) result(value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    integer, intent(in) :: min_value

    character(len=64) :: raw_value
    integer :: status, value_length, ios, parsed_value

    value = default_value
    call get_environment_variable(name, raw_value, length=value_length, &
                                  status=status)
    if (status /= 0 .or. value_length <= 0) return

    read (raw_value(1:value_length), *, iostat=ios) parsed_value
    if (ios /= 0) return

    value = max(parsed_value, min_value)
  end function env_to_integer

  subroutine init_runtime_options(comm)
    integer, intent(in) :: comm

    character(len=64) :: raw_value
    character(len=64) :: mode_value
    integer :: status, value_length
    integer :: ierr, comm_rank

    if (runtime_options_initialised) return
    runtime_options_initialised = .true.

    runtime_bench_enabled = env_to_logical("X3D2_ADIOS2_IO_BENCH", .false.)
    runtime_bench_verbose = env_to_logical("X3D2_ADIOS2_IO_BENCH_VERBOSE", &
                                           .true.)
    runtime_nvtx_enabled = env_to_logical("X3D2_ADIOS2_NVTX", .true.)
    runtime_bench_warmup_steps = env_to_integer( &
                                 "X3D2_ADIOS2_IO_BENCH_WARMUP", 2, 0)

    runtime_gpu_write_mode = gpu_write_mode_auto
    runtime_gpu_write_mode_name = "auto"
    call get_environment_variable("X3D2_ADIOS2_GPU_WRITE_MODE", raw_value, &
                                  length=value_length, status=status)
    if (status == 0 .and. value_length > 0) then
      mode_value = to_lower_ascii(adjustl(raw_value(1:value_length)))
      select case (trim(mode_value))
      case ("auto")
        runtime_gpu_write_mode = gpu_write_mode_auto
        runtime_gpu_write_mode_name = "auto"
      case ("gpu", "device", "direct")
        runtime_gpu_write_mode = gpu_write_mode_force_device
        runtime_gpu_write_mode_name = "gpu"
      case ("host", "d2h", "staged")
        runtime_gpu_write_mode = gpu_write_mode_force_host
        runtime_gpu_write_mode_name = "host"
      case default
        runtime_gpu_write_mode = gpu_write_mode_auto
        runtime_gpu_write_mode_name = "auto"
      end select
    end if

#ifndef X3D2_ADIOS2_CUDA
    if (runtime_gpu_write_mode == gpu_write_mode_force_device) then
      runtime_gpu_write_mode = gpu_write_mode_auto
      runtime_gpu_write_mode_name = "auto"
    end if
#endif

    call MPI_Comm_rank(comm, comm_rank, ierr)
    if (ierr /= 0) return

    if (comm_rank == 0 .and. .not. runtime_options_reported) then
      runtime_options_reported = .true.
#ifdef X3D2_ADIOS2_CUDA
      print '(A,A)', "ADIOS2 GPU write mode: ", &
        trim(runtime_gpu_write_mode_name)
#else
      print '(A)', "ADIOS2 GPU write mode: host (CUDA path unavailable)"
#endif
      if (runtime_bench_enabled) then
        print '(A,I0)', "ADIOS2 I/O benchmark enabled; warm-up steps: ", &
          runtime_bench_warmup_steps
      end if
      print '(A,L1)', "ADIOS2 NVTX markers enabled: ", runtime_nvtx_enabled
    end if
  end subroutine init_runtime_options

  subroutine nvtx_push_if_enabled(range_name)
    character(len=*), intent(in) :: range_name
#ifdef X3D2_ADIOS2_CUDA
    integer :: nvtx_status
    character(kind=c_char), allocatable, target :: c_name(:)
    integer :: i, n
    if (runtime_nvtx_enabled) then
      n = len_trim(range_name)
      allocate (c_name(n + 1))
      do i = 1, n
        c_name(i) = achar(iachar(range_name(i:i)), kind=c_char)
      end do
      c_name(n + 1) = c_null_char
      nvtx_status = nvtx_range_push_a(c_loc(c_name(1)))
    end if
#endif
  end subroutine nvtx_push_if_enabled

  subroutine nvtx_pop_if_enabled()
#ifdef X3D2_ADIOS2_CUDA
    integer :: nvtx_status
    if (runtime_nvtx_enabled) nvtx_status = nvtx_range_pop()
#endif
  end subroutine nvtx_pop_if_enabled

  subroutine bench_reset_step(self)
    class(io_adios2_file_t), intent(inout) :: self

    self%bench_step_put_time_local = 0.0_c_double
    self%bench_step_put_bytes_local = 0.0_c_double
  end subroutine bench_reset_step

  subroutine bench_record_put(self, put_time, put_bytes)
    class(io_adios2_file_t), intent(inout) :: self
    real(c_double), intent(in) :: put_time
    real(c_double), intent(in) :: put_bytes

    if (.not. self%bench_enabled) return

    self%bench_step_put_time_local = self%bench_step_put_time_local + &
                                     max(put_time, 0.0_c_double)
    self%bench_step_put_bytes_local = self%bench_step_put_bytes_local + &
                                      max(put_bytes, 0.0_c_double)
  end subroutine bench_record_put

  pure real(c_double) function safe_std(sum_value, sum_sq_value, count) &
    result(std_value)
    real(c_double), intent(in) :: sum_value
    real(c_double), intent(in) :: sum_sq_value
    integer, intent(in) :: count

    real(c_double) :: mean_value
    real(c_double) :: variance_value

    if (count <= 1) then
      std_value = 0.0_c_double
      return
    end if

    mean_value = sum_value/real(count, c_double)
    variance_value = sum_sq_value/real(count, c_double) - &
                     mean_value*mean_value
    std_value = sqrt(max(variance_value, 0.0_c_double))
  end function safe_std

  subroutine bench_record_step(self, end_step_time_local)
    class(io_adios2_file_t), intent(inout) :: self
    real(c_double), intent(in) :: end_step_time_local

    real(c_double) :: send_values(3)
    real(c_double) :: sum_values(3)
    real(c_double) :: max_values(3)
    real(c_double) :: step_put_time
    real(c_double) :: step_end_step_time
    real(c_double) :: step_total_time
    real(c_double) :: step_bytes
    real(c_double) :: step_throughput_gib_s
    logical :: is_warmup
    integer :: ierr, comm_rank

    if (.not. self%bench_enabled) return

    send_values(1) = self%bench_step_put_time_local
    send_values(2) = max(end_step_time_local, 0.0_c_double)
    send_values(3) = self%bench_step_put_bytes_local

    call MPI_Allreduce(send_values, sum_values, 3, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, self%comm, ierr)
    call self%handle_error(ierr, "Failed to reduce ADIOS2 bench sums")

    call MPI_Allreduce(send_values, max_values, 3, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, self%comm, ierr)
    call self%handle_error(ierr, "Failed to reduce ADIOS2 bench maxima")

    self%bench_step_counter = self%bench_step_counter + 1

    step_put_time = max_values(1)
    step_end_step_time = max_values(2)
    step_total_time = step_put_time + step_end_step_time
    step_bytes = sum_values(3)

    if (step_end_step_time > 0.0_c_double) then
      step_throughput_gib_s = step_bytes/(1024.0_c_double**3) / &
                              step_end_step_time
    else
      step_throughput_gib_s = 0.0_c_double
    end if

    is_warmup = self%bench_step_counter <= self%bench_warmup_steps
    if (.not. is_warmup) then
      self%bench_measured_steps = self%bench_measured_steps + 1
      self%bench_sum_put_time = self%bench_sum_put_time + step_put_time
      self%bench_sum_put_time2 = self%bench_sum_put_time2 + &
                                 step_put_time*step_put_time
      self%bench_sum_end_step_time = self%bench_sum_end_step_time + &
                                     step_end_step_time
      self%bench_sum_end_step_time2 = self%bench_sum_end_step_time2 + &
                                      step_end_step_time*step_end_step_time
      self%bench_sum_total_time = self%bench_sum_total_time + step_total_time
      self%bench_sum_total_time2 = self%bench_sum_total_time2 + &
                                   step_total_time*step_total_time
      self%bench_sum_throughput_gib_s = self%bench_sum_throughput_gib_s + &
                                        step_throughput_gib_s
      self%bench_sum_throughput_gib_s2 = self%bench_sum_throughput_gib_s2 + &
                                         step_throughput_gib_s* &
                                         step_throughput_gib_s
      self%bench_sum_bytes = self%bench_sum_bytes + step_bytes
    end if

    if (self%bench_verbose) then
      call MPI_Comm_rank(self%comm, comm_rank, ierr)
      call self%handle_error(ierr, "Failed to get rank for ADIOS2 bench")

      if (comm_rank == 0) then
        if (is_warmup) then
          print '(A,I0,A,ES11.3,A,ES11.3,A)', &
            "ADIOS2 bench step ", self%bench_step_counter, &
            " (warm-up): put_max=", step_put_time, &
            " s end_step_max=", step_end_step_time, " s"
        else
          print '(A,I0,A,ES11.3,A,ES11.3,A,F10.3,A)', &
            "ADIOS2 bench step ", self%bench_step_counter, &
            ": put_max=", step_put_time, &
            " s end_step_max=", step_end_step_time, &
            " s throughput=", step_throughput_gib_s, " GiB/s"
        end if
      end if
    end if

    call bench_reset_step(self)
  end subroutine bench_record_step

  subroutine bench_print_summary(self)
    class(io_adios2_file_t), intent(in) :: self

    integer :: ierr, comm_rank, n
    real(c_double) :: avg_put
    real(c_double) :: std_put
    real(c_double) :: avg_end_step
    real(c_double) :: std_end_step
    real(c_double) :: avg_total
    real(c_double) :: std_total
    real(c_double) :: avg_throughput
    real(c_double) :: std_throughput
    real(c_double) :: total_gib

    if (.not. self%bench_enabled .or. .not. self%is_writer) return

    call MPI_Comm_rank(self%comm, comm_rank, ierr)
    if (ierr /= 0 .or. comm_rank /= 0) return

    print '(A)', "ADIOS2 I/O benchmark summary"
    print '(A,A)', "  File: ", trim(self%bench_file_name)
    print '(A,A)', "  GPU write mode: ", trim(runtime_gpu_write_mode_name)
    print '(A,I0)', "  Steps seen: ", self%bench_step_counter
    print '(A,I0)', "  Warm-up steps discarded: ", self%bench_warmup_steps

    n = self%bench_measured_steps
    if (n <= 0) then
      print '(A)', "  No measured steps after warm-up discard."
      if (self%bench_close_time > 0.0_c_double) then
        print '(A,ES11.3,A)', "  adios2_close time: ", &
          self%bench_close_time, " s"
      end if
      return
    end if

    avg_put = self%bench_sum_put_time/real(n, c_double)
    std_put = safe_std(self%bench_sum_put_time, self%bench_sum_put_time2, n)

    avg_end_step = self%bench_sum_end_step_time/real(n, c_double)
    std_end_step = safe_std(self%bench_sum_end_step_time, &
                            self%bench_sum_end_step_time2, n)

    avg_total = self%bench_sum_total_time/real(n, c_double)
    std_total = safe_std(self%bench_sum_total_time, &
                         self%bench_sum_total_time2, n)

    avg_throughput = self%bench_sum_throughput_gib_s/real(n, c_double)
    std_throughput = safe_std(self%bench_sum_throughput_gib_s, &
                              self%bench_sum_throughput_gib_s2, n)

    total_gib = self%bench_sum_bytes/(1024.0_c_double**3)

    print '(A,I0)', "  Measured steps: ", n
    print '(A,ES11.3,A,ES11.3,A)', "  Put time (max rank): mean=", &
      avg_put, " s std=", std_put, " s"
    print '(A,ES11.3,A,ES11.3,A)', "  EndStep time (max rank): mean=", &
      avg_end_step, " s std=", std_end_step, " s"
    print '(A,ES11.3,A,ES11.3,A)', "  Put+EndStep (max rank): mean=", &
      avg_total, " s std=", std_total, " s"
    print '(A,F10.3,A,F10.3)', "  Throughput (GiB/s): mean=", &
      avg_throughput, " std=", std_throughput
    print '(A,F12.3,A)', "  Total data written (all ranks): ", total_gib, &
      " GiB"

    if (self%bench_close_time > 0.0_c_double) then
      print '(A,ES11.3,A)', "  adios2_close time: ", &
        self%bench_close_time, " s"
    end if
  end subroutine bench_print_summary

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
                              &before calling ADIOS2 init")

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
    temp_handle%comm = use_comm

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
    integer :: vtype
    real(c_double) :: val_dp_temp
    real(c_float) :: val_sp_temp

    select type (file_handle)
    type is (io_adios2_file_t)
      ! retrieve a variable handler within current io handler
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr == adios2_found) then
        call adios2_variable_type(vtype, var, ierr)

        if (vtype == adios2_type_dp) then
          ! file is double precision
          call adios2_get(file_handle%engine, var, val_dp_temp, &
                          adios2_mode_sync, ierr)
          value = real(val_dp_temp, dp)
        else if (vtype == adios2_type_real) then
          ! file is single precision
          call adios2_get(file_handle%engine, var, val_sp_temp, &
                          adios2_mode_sync, ierr)
          value = real(val_sp_temp, dp)
        else
          call adios2_get(file_handle%engine, var, value, &
                          adios2_mode_sync, ierr)
        end if
        call self%handle_error(ierr, "Failed to read variable " &
                               //trim(variable_name))
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
    integer :: vtype
    real(c_double), allocatable :: arr_dp_temp(:, :, :)
    real(c_float), allocatable :: arr_sp_temp(:, :, :)

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

        call adios2_variable_type(vtype, var, ierr)

        if (vtype == adios2_type_dp) then
          ! file is double precision
          allocate (arr_dp_temp(local_count(1), local_count(2), &
                                local_count(3)))
          call adios2_get(file_handle%engine, var, arr_dp_temp, &
                          adios2_mode_sync, ierr)
          array = real(arr_dp_temp, dp)
          deallocate (arr_dp_temp)
        else if (vtype == adios2_type_real) then
          ! file is single precision
          allocate (arr_sp_temp(local_count(1), local_count(2), &
                                local_count(3)))
          call adios2_get(file_handle%engine, var, arr_sp_temp, &
                          adios2_mode_sync, ierr)
          array = real(arr_sp_temp, dp)
          deallocate (arr_sp_temp)
        else
          call adios2_get(file_handle%engine, var, array, &
                          adios2_mode_sync, ierr)
        end if

        call self%handle_error(ierr, "Failed to read variable " &
                               //trim(variable_name))
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
                              &before calling ADIOS2 init")

    self%comm = comm
    call init_runtime_options(comm)
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
                             &before open")
    end if

    call adios2_open( &
      temp_handle%engine, self%io_handle, filename, &
      adios2_mode_write, use_comm, ierr)
    call self%handle_error(ierr, "Failed to open ADIOS2 engine for writing")
    temp_handle%is_writer = .true.
    temp_handle%comm = use_comm
    temp_handle%bench_file_name = trim(filename)
    temp_handle%bench_enabled = runtime_bench_enabled
    temp_handle%bench_verbose = runtime_bench_verbose
    temp_handle%bench_warmup_steps = runtime_bench_warmup_steps

    file_handle = temp_handle
  end function writer_open_adios2

  subroutine write_data_i8_adios2(self, variable_name, value, file_handle)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    integer(i8), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle

    type(adios2_variable) :: var
    integer :: ierr
    real(c_double) :: t0_put
    real(c_double) :: put_bytes

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

      if (file_handle%bench_enabled) t0_put = MPI_Wtime()
      call nvtx_push_if_enabled("ADIOS2_Put")
      call adios2_put(file_handle%engine, var, value, adios2_mode_deferred, &
                      ierr)
      call nvtx_pop_if_enabled()
      if (file_handle%bench_enabled) then
        put_bytes = bytes_integer_i8
        call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
      end if
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
    real(c_double) :: t0_put
    real(c_double) :: put_bytes

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable( &
          var, self%io_handle, variable_name, adios2_type_integer4, ierr)
        call self%handle_error(ierr, &
                               "Error defining ADIOS2 scalar integer variable")
      end if

      if (file_handle%bench_enabled) t0_put = MPI_Wtime()
      call nvtx_push_if_enabled("ADIOS2_Put")
      call adios2_put(file_handle%engine, var, value, adios2_mode_deferred, &
                      ierr)
      call nvtx_pop_if_enabled()
      if (file_handle%bench_enabled) then
        put_bytes = bytes_integer_default
        call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
      end if
      call self%handle_error(ierr, "Error writing ADIOS2 scalar integer data")
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_integer_adios2

  subroutine write_data_real_adios2(self, variable_name, value, file_handle, &
                                    use_sp)
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: value
    class(io_file_t), intent(inout) :: file_handle
    logical, intent(in), optional :: use_sp

    type(adios2_variable) :: var
    integer :: ierr
    integer :: vartype
    real(sp) :: value_sp
    logical :: convert_to_sp
    real(c_double) :: t0_put
    real(c_double) :: put_bytes

    ! Determine if we should convert to single precision
    convert_to_sp = .false.
    if (present(use_sp)) convert_to_sp = use_sp

    ! Get the appropriate ADIOS2 variable type
    vartype = get_adios2_vartype(convert_to_sp)

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable(var, self%io_handle, variable_name, &
                                    vartype, ierr)
        call self%handle_error(ierr, "Error defining ADIOS2 &
                                     &scalar real variable")
      end if

      ! Write data - convert to single precision if needed
      if (convert_to_sp .and. .not. is_sp) then
        value_sp = real(value, sp)
        ! Use sync mode to ensure data is copied before value_sp goes out of scope
        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put(file_handle%engine, var, value_sp, adios2_mode_sync, &
                        ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = bytes_real_sp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        call self%handle_error(ierr, "Error writing ADIOS2 scalar &
                                     &single precision real data")
      else
        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put(file_handle%engine, var, value, &
                        adios2_mode_deferred, ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = bytes_real_dp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        call self%handle_error(ierr, "Error writing ADIOS2 scalar real data")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_real_adios2

  subroutine write_data_array_3d_adios2( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp

    type(adios2_variable) :: var
    integer :: ierr
    integer :: vartype
    real(sp), allocatable :: array_sp(:, :, :)
    logical :: convert_to_sp
    real(c_double) :: t0_put
    real(c_double) :: put_bytes

    ! Determine if we should convert to single precision
    convert_to_sp = .false.
    if (present(use_sp)) convert_to_sp = use_sp

    ! Get the appropriate ADIOS2 variable type
    vartype = get_adios2_vartype(convert_to_sp)

    select type (file_handle)
    type is (io_adios2_file_t)
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)

      if (ierr /= adios2_found) then
        call adios2_define_variable(var, self%io_handle, variable_name, &
                                    vartype, 3, shape_dims, &
                                    start_dims, count_dims, &
                                    adios2_constant_dims, ierr)
        call self%handle_error(ierr, "Error defining ADIOS2 &
                                     &3D array real variable")
      end if

      ! Write data - convert to single precision if needed
      if (convert_to_sp .and. .not. is_sp) then
        ! Allocate temporary single precision buffer
        allocate (array_sp(size(array, 1), size(array, 2), size(array, 3)))
        array_sp = real(array, sp)
        ! Use sync mode to ensure data is copied before buffer is deallocated
        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put(file_handle%engine, var, array_sp, adios2_mode_sync, &
                        ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = real(size(array_sp, kind=i8), c_double) * &
                      bytes_real_sp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        deallocate (array_sp)
        call self%handle_error(ierr, "Error writing ADIOS2 3D array &
                                     &single precision real data")
      else
        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put(file_handle%engine, var, array, &
                        adios2_mode_deferred, ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = real(size(array, kind=i8), c_double) * &
                      bytes_real_dp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        call self%handle_error(ierr, "Error writing ADIOS2 &
                                     &3D array real data")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_array_3d_adios2

#ifdef X3D2_ADIOS2_CUDA
  subroutine write_data_array_3d_device_adios2( &
    self, variable_name, array, file_handle, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    !! GPU-aware I/O: passes device pointer directly to ADIOS2.
    !! Uses the native Fortran adios2_set_memory_space to tell ADIOS2
    !! the buffer is on GPU, then calls a C wrapper for adios2_put
    !! (the Fortran generic adios2_put does not accept device arrays).
    !! Requires ADIOS2 built with -DADIOS2_USE_CUDA=ON.
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    real(dp), device, target, intent(in) :: array(:, :, :)
    class(io_file_t), intent(inout) :: file_handle
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp

    type(adios2_variable) :: var
    type(c_devptr) :: devptr
    type(c_ptr) :: device_ptr
    integer :: ierr, vartype
    real(sp), allocatable, device, target :: array_sp(:, :, :)
    logical :: convert_to_sp
    real(c_double) :: t0_put
    real(c_double) :: put_bytes

    convert_to_sp = .false.
    if (present(use_sp)) convert_to_sp = use_sp
    vartype = get_adios2_vartype(convert_to_sp)

    select type (file_handle)
    type is (io_adios2_file_t)
      ! Define or inquire variable
      call adios2_inquire_variable(var, self%io_handle, variable_name, ierr)
      if (ierr /= adios2_found) then
        call adios2_define_variable(var, self%io_handle, variable_name, &
                                    vartype, 3, shape_dims, &
                                    start_dims, count_dims, &
                                    adios2_constant_dims, ierr)
        call self%handle_error(ierr, "Error defining ADIOS2 variable (GPU)")
      end if

      ! Tell ADIOS2 this variable holds GPU memory
      call adios2_set_memory_space(var, adios2_memory_space_gpu, ierr)
      call self%handle_error(ierr, "Error setting GPU memory space")

      ! Get device pointer and pass via C wrapper
      ! Use sync mode only when converting dp->sp (temporary buffer must
      ! be flushed before deallocation); deferred mode otherwise.
      if (convert_to_sp .and. .not. is_sp) then
        allocate (array_sp(size(array, 1), size(array, 2), size(array, 3)))
        array_sp = real(array, sp)
        devptr = c_devloc(array_sp)
        device_ptr = transfer(devptr, device_ptr)

        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put_gpu(file_handle%engine%f2c, var%f2c, &
                            device_ptr, adios2_mode_sync, ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = real(size(array_sp, kind=i8), c_double) * &
                      bytes_real_sp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        deallocate (array_sp)
        call self%handle_error(ierr, "Error in GPU-aware ADIOS2 put (sync)")
      else
        devptr = c_devloc(array)
        device_ptr = transfer(devptr, device_ptr)

        if (file_handle%bench_enabled) t0_put = MPI_Wtime()
        call nvtx_push_if_enabled("ADIOS2_Put")
        call adios2_put_gpu(file_handle%engine%f2c, var%f2c, &
                            device_ptr, adios2_mode_deferred, ierr)
        call nvtx_pop_if_enabled()
        if (file_handle%bench_enabled) then
          put_bytes = real(size(array, kind=i8), c_double) * &
                      bytes_real_dp
          call bench_record_put(file_handle, MPI_Wtime() - t0_put, put_bytes)
        end if
        call self%handle_error(ierr, "Error in GPU-aware ADIOS2 &
                              &put (deferred)")
      end if
    class default
      call self%handle_error(1, "Invalid file handle type for ADIOS2")
    end select
  end subroutine write_data_array_3d_device_adios2
#endif

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
                             &array attribute " &
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
    real(c_double) :: t0_close

    if (self%is_step_active) call self%end_step()

    if (self%engine%valid) then
      if (self%is_writer .and. self%bench_enabled) t0_close = MPI_Wtime()
      call adios2_close(self%engine, ierr)
      if (self%is_writer .and. self%bench_enabled) then
        self%bench_close_time = self%bench_close_time + &
                                (MPI_Wtime() - t0_close)
      end if
      call self%handle_error(ierr, "Failed to close ADIOS2 engine")
    end if

    if (self%is_writer .and. self%bench_enabled) then
      call bench_print_summary(self)
    end if

    if (self%nvtx_step_range_active) then
      call nvtx_pop_if_enabled()
      self%nvtx_step_range_active = .false.
    end if
  end subroutine file_close_adios2

  subroutine file_begin_step_adios2(self)
    class(io_adios2_file_t), intent(inout) :: self

    integer :: ierr

    if (self%is_step_active) return

    if (self%is_writer) then
      call adios2_begin_step(self%engine, adios2_step_mode_append, ierr)
      call self%handle_error(ierr, "Error beginning ADIOS2 step for writing")
      if (self%bench_enabled) call bench_reset_step(self)
      if (.not. self%nvtx_step_range_active) then
        call nvtx_push_if_enabled("ADIOS2_I/O_Step")
        self%nvtx_step_range_active = .true.
      end if
    else
      call adios2_begin_step(self%engine, adios2_step_mode_read, ierr)
      call self%handle_error(ierr, "Error beginning ADIOS2 step for reading")
    end if

    self%is_step_active = .true.
  end subroutine file_begin_step_adios2

  subroutine file_end_step_adios2(self)
    class(io_adios2_file_t), intent(inout) :: self

    integer :: ierr
    real(c_double) :: t0_end_step
    real(c_double) :: end_step_time_local

    if (.not. self%is_step_active) return

    if (self%is_writer .and. self%bench_enabled) t0_end_step = MPI_Wtime()
    if (self%is_writer) call nvtx_push_if_enabled("ADIOS2_EndStep")
    call adios2_end_step(self%engine, ierr)
    if (self%is_writer) call nvtx_pop_if_enabled()
    call self%handle_error(ierr, "Failed to end ADIOS2 step")

    if (self%is_writer .and. self%bench_enabled) then
      end_step_time_local = MPI_Wtime() - t0_end_step
      call bench_record_step(self, end_step_time_local)
    end if

    if (self%nvtx_step_range_active) then
      call nvtx_pop_if_enabled()
      self%nvtx_step_range_active = .false.
    end if

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

  logical function supports_device_field_write_adios2(self)
    class(io_adios2_writer_t), intent(in) :: self
    call init_runtime_options(self%comm)
#ifdef X3D2_ADIOS2_CUDA
    supports_device_field_write_adios2 = &
      runtime_gpu_write_mode /= gpu_write_mode_force_host
#else
    supports_device_field_write_adios2 = .false.
#endif
  end function supports_device_field_write_adios2

#ifdef X3D2_ADIOS2_CUDA
  subroutine sync_device_adios2(self)
    !! Synchronise the GPU device before a batch of I/O operations.
    !! Call once before writing multiple fields to avoid per-field syncs.
    class(io_adios2_writer_t), intent(inout) :: self
    integer :: ierr
    ierr = cudaDeviceSynchronize()
    if (ierr /= 0) error stop "cudaDeviceSynchronize failed before I/O"
  end subroutine sync_device_adios2
#endif

  subroutine write_field_from_solver_adios2( &
    self, variable_name, field, file_handle, backend, &
    shape_dims, start_dims, count_dims, use_sp &
    )
    !! Write field with automatic GPU-aware optimisation when available
    use m_field, only: field_t
#ifdef X3D2_ADIOS2_CUDA
    use m_cuda_allocator, only: cuda_field_t
#endif
    class(io_adios2_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: variable_name
    class(*), intent(in) :: field
    class(io_file_t), intent(inout) :: file_handle
    class(*), intent(in) :: backend
    integer(i8), intent(in) :: shape_dims(3)
    integer(i8), intent(in) :: start_dims(3)
    integer(i8), intent(in) :: count_dims(3)
    logical, intent(in), optional :: use_sp

    call init_runtime_options(self%comm)

#ifdef X3D2_ADIOS2_CUDA
    if (runtime_gpu_write_mode /= gpu_write_mode_force_host) then
      ! GPU-aware ADIOS2: write directly from device memory
      select type (field_typed => field)
      type is (cuda_field_t)
        call self%write_data_array_3d_device( &
          variable_name, field_typed%data_d, &
          file_handle, shape_dims, start_dims, count_dims, use_sp &
          )
        return
      end select
    end if
#endif

    ! Non-CUDA backend: standard host path
    select type (field_typed => field)
    type is (field_t)
      call self%write_data_array_3d( &
        variable_name, field_typed%data, file_handle, &
        shape_dims, start_dims, count_dims, use_sp &
        )
    class default
      error stop "write_field_from_solver: Unsupported field type"
    end select
  end subroutine write_field_from_solver_adios2

end module m_io_backend
