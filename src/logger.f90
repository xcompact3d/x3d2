module m_logger
  !! Single point for x3d2 to report progress, warnings and errors, and to
  !! stop the run when something goes wrong.
  !!
  !! Messages carry a level: debug, info, warning or error. Those below the
  !! threshold set with `set_log_level` (info by default) are dropped. Debug
  !! and info lines go to standard output, warnings and errors to standard
  !! error, and every line is flushed so output from different ranks and
  !! from a crash is not lost in a buffer.
  !!
  !! The logging calls are made on every rank, but by default only rank 0
  !! writes, so a message appears once however many ranks the run has. Pass
  !! `all_ranks=.true.` for something that only some ranks know about; each
  !! rank that makes the call then writes, and on a multi-rank run the line
  !! starts with the rank that wrote it.
  !!
  !! `abort_run` is how x3d2 stops on an error. It writes the message and
  !! takes down the whole job: a plain `error stop` on one rank can leave the
  !! others waiting forever in MPI. `check_mpi` turns an MPI return code into
  !! such an abort. Backends add their own checks on top of `abort_run` (see
  !! m_cuda_errors for CUDA).
  !!
  !! `str` turns numbers into text, so a message can be assembled with `//`:
  !!
  !!   call log_info('Checkpoint frequency: '//str(freq))
  use iso_fortran_env, only: int32, int64, real32, real64, &
                             stdout => output_unit, stderr => error_unit

  use m_mpi, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_MAX_ERROR_STRING, &
                   MPI_Initialized, MPI_Finalized, MPI_Comm_rank, &
                   MPI_Comm_size, MPI_Barrier, MPI_Abort, MPI_Error_string

  implicit none

  private
  public :: LOG_LEVEL_DEBUG, LOG_LEVEL_INFO, LOG_LEVEL_WARN, LOG_LEVEL_ERROR
  public :: set_log_level, get_log_level
  public :: log_debug, log_info, log_warn, log_error
  public :: abort_run, check_mpi
  public :: str

  integer, parameter :: LOG_LEVEL_DEBUG = 10
  integer, parameter :: LOG_LEVEL_INFO = 20
  integer, parameter :: LOG_LEVEL_WARN = 30
  integer, parameter :: LOG_LEVEL_ERROR = 40

  integer :: log_level = LOG_LEVEL_INFO

  !> Exit code of a run stopped by abort_run.
  integer, parameter :: abort_code = 1

  !> Default edit descriptors `str` uses for reals: enough significant digits
  !> (9 and 17) that no value loses precision when printed.
  character(len=*), parameter :: real32_fmt = '(es15.8)'
  character(len=*), parameter :: real64_fmt = '(es24.16e3)'

  interface str
    module procedure str_int32
    module procedure str_int64
    module procedure str_real32
    module procedure str_real64
    module procedure str_logical
  end interface str

contains

  subroutine set_log_level(level)
    !! Drop messages below `level` (one of the LOG_LEVEL_* parameters).
    integer, intent(in) :: level

    log_level = level
  end subroutine set_log_level

  integer function get_log_level()
    get_log_level = log_level
  end function get_log_level

  subroutine log_debug(message, all_ranks)
    !! Detail for developers, hidden unless the level is LOG_LEVEL_DEBUG.
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: all_ranks

    call write_message(LOG_LEVEL_DEBUG, 'DEBUG: ', message, all_ranks)
  end subroutine log_debug

  subroutine log_info(message, all_ranks)
    !! Progress and setup information.
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: all_ranks

    call write_message(LOG_LEVEL_INFO, '', message, all_ranks)
  end subroutine log_info

  subroutine log_warn(message, all_ranks)
    !! Something is suspicious, but the run can carry on.
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: all_ranks

    call write_message(LOG_LEVEL_WARN, 'WARNING: ', message, all_ranks)
  end subroutine log_warn

  subroutine log_error(message, all_ranks)
    !! An error the caller handles itself. Use abort_run to stop the run.
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: all_ranks

    call write_message(LOG_LEVEL_ERROR, 'ERROR: ', message, all_ranks)
  end subroutine log_error

  subroutine abort_run(message, collective)
    !! Report a fatal error and stop every rank of the run.
    !!
    !! By default the rank that calls this writes the message, prefixed with
    !! its rank on a multi-rank run, and aborts the job. That is safe however
    !! many ranks get here, so it is the right choice for a failure only some
    !! ranks see: a CUDA error, a missing file, a bad MPI return code.
    !!
    !! Pass `collective=.true.` when every rank is certain to reach this call
    !! with the same message, typically while checking the input. Rank 0 then
    !! writes the message once and aborts, while the other ranks wait for it.
    !! Never use it for an error some ranks may not see: those ranks would
    !! carry on while the rest wait, and the run hangs.
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: collective

    logical :: is_collective, mpi_running
    integer :: rank, nproc, ierr

    is_collective = .false.
    if (present(collective)) is_collective = collective

    call query_mpi(mpi_running, rank, nproc)

    if (is_collective .and. rank /= 0) then
      ! Rank 0 is reporting this same error and will abort the job. Waiting
      ! here, rather than aborting too, keeps this rank from killing the job
      ! before rank 0 has written the message. The barrier never completes
      ! because rank 0 does not join it.
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end if

    call emit(stderr, line_prefix(.not. is_collective, rank, nproc) &
              //'ERROR: '//message)

    if (mpi_running .and. nproc > 1) then
      call MPI_Abort(MPI_COMM_WORLD, abort_code, ierr)
    end if
    ! A single rank, or MPI not running: stop this process.
    error stop abort_code
  end subroutine abort_run

  subroutine check_mpi(ierr, operation)
    !! Abort the run if an MPI call did not succeed.
    !!
    !! `operation` says what was being done, for example 'MPI_Waitall in
    !! sendrecv_fields', and appears in the message with MPI's own
    !! description of the error.
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation

    character(len=MPI_MAX_ERROR_STRING) :: description
    integer :: description_len, ierr_string

    if (ierr == MPI_SUCCESS) return

    call MPI_Error_string(ierr, description, description_len, ierr_string)
    if (ierr_string /= MPI_SUCCESS) then
      description = 'unknown error'
      description_len = len_trim(description)
    end if
    call abort_run(operation//' failed: '//description(1:description_len) &
                   //' (MPI error code '//str(ierr)//')')
  end subroutine check_mpi

  subroutine write_message(level, tag, message, all_ranks)
    integer, intent(in) :: level
    character(len=*), intent(in) :: tag, message
    logical, intent(in), optional :: all_ranks

    logical :: every_rank, mpi_running
    integer :: rank, nproc, unit

    if (level < log_level) return

    every_rank = .false.
    if (present(all_ranks)) every_rank = all_ranks

    call query_mpi(mpi_running, rank, nproc)
    if (.not. every_rank .and. rank /= 0) return

    unit = stdout
    if (level >= LOG_LEVEL_WARN) unit = stderr
    call emit(unit, line_prefix(every_rank, rank, nproc)//tag//message)
  end subroutine write_message

  subroutine query_mpi(mpi_running, rank, nproc)
    !! Rank and size of MPI_COMM_WORLD, which can only be asked between
    !! MPI_Init and MPI_Finalize. Outside that window this process counts as
    !! a single rank 0.
    logical, intent(out) :: mpi_running
    integer, intent(out) :: rank, nproc

    logical :: initialised, finalised
    integer :: ierr

    rank = 0
    nproc = 1
    call MPI_Initialized(initialised, ierr)
    call MPI_Finalized(finalised, ierr)
    mpi_running = initialised .and. .not. finalised
    if (mpi_running) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    end if
  end subroutine query_mpi

  function line_prefix(show_rank, rank, nproc) result(prefix)
    !! Names the writing rank, when there is more than one to tell apart.
    logical, intent(in) :: show_rank
    integer, intent(in) :: rank, nproc
    character(len=:), allocatable :: prefix

    if (show_rank .and. nproc > 1) then
      prefix = '[rank '//str(rank)//'] '
    else
      prefix = ''
    end if
  end function line_prefix

  subroutine emit(unit, line)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: line

    write (unit, '(a)') line
    flush (unit)
  end subroutine emit

  function str_int32(value) result(text)
    integer(int32), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=11) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)
  end function str_int32

  function str_int64(value) result(text)
    integer(int64), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=20) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)
  end function str_int64

  function str_real32(value, fmt) result(text)
    !! `fmt` is a format string such as '(f8.3)'; the default is scientific
    !! notation with nine significant digits.
    real(real32), intent(in) :: value
    character(len=*), intent(in), optional :: fmt
    character(len=:), allocatable :: text

    character(len=64) :: buffer

    if (present(fmt)) then
      write (buffer, fmt) value
    else
      write (buffer, real32_fmt) value
    end if
    text = trim(adjustl(buffer))
  end function str_real32

  function str_real64(value, fmt) result(text)
    !! `fmt` is a format string such as '(f8.3)'; the default is scientific
    !! notation with seventeen significant digits.
    real(real64), intent(in) :: value
    character(len=*), intent(in), optional :: fmt
    character(len=:), allocatable :: text

    character(len=64) :: buffer

    if (present(fmt)) then
      write (buffer, fmt) value
    else
      write (buffer, real64_fmt) value
    end if
    text = trim(adjustl(buffer))
  end function str_real64

  function str_logical(value) result(text)
    logical, intent(in) :: value
    character(len=:), allocatable :: text

    text = merge('T', 'F', value)
  end function str_logical

end module m_logger
