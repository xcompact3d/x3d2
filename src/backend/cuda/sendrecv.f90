module m_cuda_sendrecv
  use cudafor
  use mpi

  use m_common, only: dp

  implicit none

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e
    real(dp), device, dimension(:, :, :), intent(in) :: f_send_s, f_send_e
    integer, intent(in) :: n_data, nproc, prev, next

    integer :: req(4), err(4), ierr, tag = 1234

    if (nproc == 1) then
      f_recv_s = f_send_e
      f_recv_e = f_send_s
    else
      call MPI_Isend(f_send_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f_recv_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f_send_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f_recv_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)
    end if

  end subroutine sendrecv_fields

  subroutine sendrecv_3fields( &
    f1_recv_s, f1_recv_e, f2_recv_s, f2_recv_e, f3_recv_s, f3_recv_e, &
    f1_send_s, f1_send_e, f2_send_s, f2_send_e, f3_send_s, f3_send_e, &
    n_data, nproc, prev, next &
    )
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: &
      f1_recv_s, f1_recv_e, f2_recv_s, f2_recv_e, f3_recv_s, f3_recv_e
    real(dp), device, dimension(:, :, :), intent(in) :: &
      f1_send_s, f1_send_e, f2_send_s, f2_send_e, f3_send_s, f3_send_e
    integer, intent(in) :: n_data, nproc, prev, next

    integer :: req(12), err(12), ierr, tag = 1234

    if (nproc == 1) then
      f1_recv_s = f1_send_e
      f1_recv_e = f1_send_s
      f2_recv_s = f2_send_e
      f2_recv_e = f2_send_s
      f3_recv_s = f3_send_e
      f3_recv_e = f3_send_s
    else
      call MPI_Isend(f1_send_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f1_recv_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f1_send_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f1_recv_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      call MPI_Isend(f2_send_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(5), err(5))
      call MPI_Irecv(f2_recv_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(6), err(6))
      call MPI_Isend(f2_send_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(7), err(7))
      call MPI_Irecv(f2_recv_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(8), err(8))

      call MPI_Isend(f3_send_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(9), err(9))
      call MPI_Irecv(f3_recv_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(10), err(10))
      call MPI_Isend(f3_send_e, n_data, MPI_DOUBLE_PRECISION, &
                     next, tag, MPI_COMM_WORLD, req(11), err(11))
      call MPI_Irecv(f3_recv_s, n_data, MPI_DOUBLE_PRECISION, &
                     prev, tag, MPI_COMM_WORLD, req(12), err(12))

      call MPI_Waitall(12, req, MPI_STATUSES_IGNORE, ierr)
    end if

  end subroutine sendrecv_3fields

end module m_cuda_sendrecv
