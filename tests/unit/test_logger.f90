program test_logger
  !! Checks the text `str` produces, the level threshold, and that logging
  !! from every rank, from rank 0 only and below the threshold all return
  !! normally on any number of ranks.
  use iso_fortran_env, only: int64, real32, real64

  use m_logger, only: LOG_LEVEL_DEBUG, LOG_LEVEL_INFO, LOG_LEVEL_WARN, &
                      set_log_level, get_log_level, &
                      log_debug, log_info, log_warn, log_error, &
                      check_mpi, str
  use m_mpi, only: MPI_SUCCESS
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  integer :: nrank, nproc
  logical :: allpass

  call initialise_mpi(nrank, nproc)
  allpass = .true.

  call check_text('int32', str(42), '42')
  call check_text('negative int32', str(-7), '-7')
  call check_text('int64', str(9000000000_int64), '9000000000')
  call check_text('real64 default', str(0.5_real64), '5.000000E-01')
  call check_text('real32 default', str(-1.25_real32), '-1.250000E+00')
  call check_text('real64 with format', str(3.14159_real64, '(f6.3)'), &
                  '3.142')
  call check_text('logical', str(.true.)//str(.false.), 'TF')

  call check(get_log_level() == LOG_LEVEL_INFO, 'default level is info')
  call set_log_level(LOG_LEVEL_DEBUG)
  call check(get_log_level() == LOG_LEVEL_DEBUG, 'level can be lowered')
  call log_debug('test_logger: debug line, shown at debug level')

  call set_log_level(LOG_LEVEL_WARN)
  call log_info('test_logger: info line, must not appear')
  call set_log_level(LOG_LEVEL_INFO)

  call log_info('test_logger: info line from rank 0')
  call log_info('test_logger: info line from rank '//str(nrank), &
                all_ranks=.true.)
  call log_warn('test_logger: warning line from rank 0')
  call log_error('test_logger: non-fatal error line from rank 0')

  ! Success returns without stopping the run.
  call check_mpi(MPI_SUCCESS, 'test_logger')

  call finalise_test(allpass, nrank)

contains

  subroutine check_text(label, got, expected)
    character(len=*), intent(in) :: label, got, expected

    call check(got == expected, label//": got '"//got//"', expected '" &
               //expected//"'")
  end subroutine check_text

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      print *, 'FAIL: ', label
      allpass = .false.
    end if
  end subroutine check

end program test_logger
