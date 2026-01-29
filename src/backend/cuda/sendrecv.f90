module m_cuda_sendrecv
  !! MPI communication for CUDA backend using device pointers.
  !!
  !! Passes device pointers directly to MPI calls. With GPU-aware MPI
  !! implementations (e.g., OpenMPI with CUDA support, MVAPICH2-GDR),
  !! data transfers directly between GPU memories without staging through
  !! host, reducing latency and increasing bandwidth.
  !!
  !! Without GPU-aware MPI, the implementation may stage through host
  !! memory automatically, still functional but with additional overhead.
  !!
  !! - `sendrecv_fields`: Single field halo exchange
  !! - `sendrecv_3fields`: Batch exchange for three fields (velocity components
  !!   or derivatives). Batching amortises MPI overhead and enables better
  !!   network utilisation.
  use cudafor
  use mpi

  use m_common, only: dp, MPI_X3D2_DP

  implicit none

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    !! Exchange boundary halos using MPI with device pointers.
    !!
    !! MPI_Isend/Irecv allows all four communications (send to prev/next,
    !! receive from prev/next) to proceed concurrently, enabling network
    !! pipelining. MPI_Waitall synchronises only when results needed.
    !!
    !! When nproc=1, data copied directly on device without MPI.
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e  !! Device receive buffers
    real(dp), device, dimension(:, :, :), intent(in) :: f_send_s, f_send_e   !! Device send buffers
    integer, intent(in) :: n_data    !! Number of data elements
    integer, intent(in) :: nproc     !! Number of processes in direction
    integer, intent(in) :: prev      !! Previous neighbour rank
    integer, intent(in) :: next      !! Next neighbour rank

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

  subroutine sendrecv_3fields( &
    f1_recv_s, f1_recv_e, f2_recv_s, f2_recv_e, f3_recv_s, f3_recv_e, &
    f1_send_s, f1_send_e, f2_send_s, f2_send_e, f3_send_s, f3_send_e, &
    n_data, nproc, prev, next &
    )
    !! Exchange three fields simultaneously using batched MPI communication.
    !!
    !! Used for: (1) velocity component halos (u, v, w) before computing transport
    !! equation, (2) derivative field halos (du, dud, d2u) in distributed compact
    !! schemes. Batching all three fields amortises MPI setup overhead. Single
    !! MPI_Waitall for all 12 operations reduces synchronisation points.
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: &
      f1_recv_s, f1_recv_e, f2_recv_s, f2_recv_e, f3_recv_s, f3_recv_e  !! Device receive buffers
    real(dp), device, dimension(:, :, :), intent(in) :: &
      f1_send_s, f1_send_e, f2_send_s, f2_send_e, f3_send_s, f3_send_e  !! Device send buffers
    integer, intent(in) :: n_data    !! Number of data elements per field
    integer, intent(in) :: nproc     !! Number of processes
    integer, intent(in) :: prev      !! Previous neighbour rank
    integer, intent(in) :: next      !! Next neighbour rank

    integer :: req(12), err(12), ierr, tag = 1234

    if (nproc == 1) then
      f1_recv_s = f1_send_e
      f1_recv_e = f1_send_s
      f2_recv_s = f2_send_e
      f2_recv_e = f2_send_s
      f3_recv_s = f3_send_e
      f3_recv_e = f3_send_s
    else
      call MPI_Isend(f1_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f1_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f1_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f1_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      call MPI_Isend(f2_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(5), err(5))
      call MPI_Irecv(f2_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(6), err(6))
      call MPI_Isend(f2_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(7), err(7))
      call MPI_Irecv(f2_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(8), err(8))

      call MPI_Isend(f3_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(9), err(9))
      call MPI_Irecv(f3_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(10), err(10))
      call MPI_Isend(f3_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(11), err(11))
      call MPI_Irecv(f3_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(12), err(12))

      call MPI_Waitall(12, req, MPI_STATUSES_IGNORE, ierr)
    end if

  end subroutine sendrecv_3fields

end module m_cuda_sendrecv
