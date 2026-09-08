program test_adm_checkpoint
  !! Tests the ADIOS2 round-trip of restart-critical ADM scalar state.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_adm, only: adm_t
  use m_common, only: dp
  use m_io_session, only: reader_session_t, writer_session_t

  implicit none

  character(len=*), parameter :: checkpoint_file = 'test_adm_state.bp'
  type(adm_t) :: adm
  type(reader_session_t) :: reader
  type(writer_session_t) :: writer
  integer :: ierr, nrank
  logical :: allpass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)

  adm%n_turb = 1
  adm%sample_count = 7
  allocate (adm%disc(1))
  adm%disc(1)%U_disc = 9._dp
  adm%disc(1)%U_disc_filt = 8._dp
  adm%disc(1)%power = 700._dp
  adm%disc(1)%thrust = 80._dp
  adm%disc(1)%U_disc_mean = 7._dp
  adm%disc(1)%power_mean = 600._dp
  adm%disc(1)%thrust_mean = 70._dp

  call writer%open(checkpoint_file, MPI_COMM_WORLD)
  call adm%write_checkpoint(writer)
  call writer%close()

  adm%sample_count = 0
  adm%disc(1)%U_disc = 0._dp
  adm%disc(1)%U_disc_filt = 0._dp
  adm%disc(1)%power = 0._dp
  adm%disc(1)%thrust = 0._dp
  adm%disc(1)%U_disc_mean = 0._dp
  adm%disc(1)%power_mean = 0._dp
  adm%disc(1)%thrust_mean = 0._dp

  call reader%open(checkpoint_file, MPI_COMM_WORLD)
  call adm%read_checkpoint(reader)
  call reader%close()

  allpass = adm%sample_count == 7 &
            .and. adm%filter_initialized &
            .and. adm%disc(1)%U_disc == 9._dp &
            .and. adm%disc(1)%U_disc_filt == 8._dp &
            .and. adm%disc(1)%power == 700._dp &
            .and. adm%disc(1)%thrust == 80._dp &
            .and. adm%disc(1)%U_disc_mean == 7._dp &
            .and. adm%disc(1)%power_mean == 600._dp &
            .and. adm%disc(1)%thrust_mean == 70._dp

  if (nrank == 0) call execute_command_line('rm -rf '//checkpoint_file)
  deallocate (adm%disc)
  call MPI_Finalize(ierr)

  if (nrank == 0) then
    if (allpass) then
      write (stderr, '(A)') 'ADM checkpoint test passed.'
    else
      error stop 'ADM checkpoint test failed.'
    end if
  end if

end program test_adm_checkpoint
