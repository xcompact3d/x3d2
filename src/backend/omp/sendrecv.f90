module m_omp_sendrecv
  !! MPI halo exchange utilities for OMP backend.
  !!
  !! Provides non-blocking point-to-point communication for exchanging
  !! boundary halos between neighbouring MPI processes. Used in distributed
  !! compact finite difference schemes that require off-process data.
  !!
  !! **Communication pattern:** Bidirectional simultaneous send/recv with
  !! neighbours in one decomposition direction.
  !!
  !! **Single-process optimisation:** Direct copy when no MPI communication needed.
  use mpi

  use m_common, only: dp, MPI_X3D2_DP

  implicit none

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    !! Exchange boundary halos with neighbouring MPI processes.
    !!
    !! Performs bidirectional halo exchange using non-blocking MPI
    !! communication (MPI_Isend/MPI_Irecv). Sends data to both neighbours
    !! simultaneously and receives from both, then waits for all operations
    !! to complete.
    !!
    !! **Special case:** Single-process (nproc=1) uses direct memory copy
    !! for periodic boundaries without MPI overhead.
    !!
    !! **Communication pattern:**
    !! - Send start halo to previous process
    !! - Receive end halo from next process
    !! - Send end halo to next process
    !! - Receive start halo from previous process
    !!
    !! **Non-blocking:** All 4 operations initiated before waiting for completion.
    implicit none

    real(dp), dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e  !! Receive buffers (start/end halos)
    real(dp), dimension(:, :, :), intent(in) :: f_send_s, f_send_e   !! Send buffers (start/end halos)
    integer, intent(in) :: n_data    !! Number of data elements to transfer
    integer, intent(in) :: nproc     !! Number of processes in this direction
    integer, intent(in) :: prev      !! Rank of previous neighbour
    integer, intent(in) :: next      !! Rank of next neighbour

    integer :: req(4), err(4), ierr, tag = 1234

    if (nproc == 1) then
      f_recv_s = f_send_e
      f_recv_e = f_send_s
    else
      call MPI_Isend(f_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)
    end if

  end subroutine sendrecv_fields

end module m_omp_sendrecv
