program test_logger_abort
  !! Calls abort_run and must not return from it.
  !!
  !! ctest passes this test when the deliberate error message appears in the
  !! output, so it fails if the run hangs (timeout), crashes before reporting,
  !! or carries on past the abort. Run on several ranks it shows the abort
  !! takes down every rank:
  !!
  !!   - no argument: only the last rank fails, the way a CUDA or I/O error
  !!     strikes one rank while the others carry on into MPI;
  !!   - `collective`: every rank reaches abort_run with the same message,
  !!     the way an input check does.
  use m_logger, only: abort_run, str
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier
  use m_test_utils, only: initialise_mpi

  implicit none

  integer :: nrank, nproc, ierr
  character(len=32) :: mode

  call initialise_mpi(nrank, nproc)
  mode = ''
  if (command_argument_count() >= 1) call get_command_argument(1, mode)

  if (trim(mode) == 'collective') then
    call abort_run('deliberate collective abort from test_logger_abort', &
                   collective=.true.)
  else
    if (nrank == nproc - 1) then
      call abort_run('deliberate abort from rank '//str(nrank) &
                     //' in test_logger_abort')
    end if
    ! The other ranks wait for a partner that never arrives, as they would
    ! in a halo exchange with a failed rank.
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end if

  print *, 'abort_run returned: this line must never be reached'

end program test_logger_abort
