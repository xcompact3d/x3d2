module m_omp_sendrecv
  use mpi

  use m_common, only: dp, MPI_X3D2_DP

  implicit none

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    implicit none

    real(dp), dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e
    real(dp), dimension(:, :, :), intent(in) :: f_send_s, f_send_e
    integer, intent(in) :: n_data, nproc, prev, next

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
