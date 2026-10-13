module m_cuda_sendrecv
  use cudafor
  use m_mpi, only: MPI_COMM_WORLD, MPI_STATUSES_IGNORE, MPI_Irecv, MPI_Isend, &
                   MPI_Waitall

  use m_common, only: dp, MPI_X3D2_DP

  implicit none

  integer(kind=cuda_stream_kind), parameter :: default_stream = 0

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e
    real(dp), device, dimension(:, :, :), intent(in) :: f_send_s, f_send_e
    integer, intent(in) :: n_data, nproc, prev, next

    integer :: req(4), err(4), ierr, tag = 1234

    if (nproc == 1) then
      call copy_halo(f_recv_s, f_send_e, n_data)
      call copy_halo(f_recv_e, f_send_s, n_data)
    else
      ! Ensure the send buffers are ready before MPI touches device memory.
      ierr = cudaStreamSynchronize(default_stream)

#ifdef MPI
      call MPI_Isend(f_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)
#endif
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
      call copy_halo(f1_recv_s, f1_send_e, n_data)
      call copy_halo(f1_recv_e, f1_send_s, n_data)
      call copy_halo(f2_recv_s, f2_send_e, n_data)
      call copy_halo(f2_recv_e, f2_send_s, n_data)
      call copy_halo(f3_recv_s, f3_send_e, n_data)
      call copy_halo(f3_recv_e, f3_send_s, n_data)
    else
      ! Ensure the send buffers are ready before MPI touches device memory.
      ierr = cudaStreamSynchronize(default_stream)

#ifdef MPI
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
#endif
    end if

  end subroutine sendrecv_3fields

  subroutine copy_halo(dst, src, n_data)
    !! Periodic wrap-around on a single rank: copies the first n_data
    !! elements of a halo buffer, the groups of the direction being
    !! exchanged. The buffers are sized for the direction with the most
    !! groups, so a whole-array assignment would copy that maximum on every
    !! exchange (16 times the data in use on an elongated mesh).
    implicit none

    real(dp), device, contiguous, dimension(:, :, :), intent(inout) :: dst
    real(dp), device, contiguous, dimension(:, :, :), intent(in) :: src
    integer, intent(in) :: n_data

    integer :: ierr

    if (n_data > size(src) .or. n_data > size(dst)) then
      error stop 'copy_halo: n_data exceeds the halo buffer.'
    end if

    ! Default stream, so it stays ordered with the kernels around it
    ierr = cudaMemcpyAsync(dst, src, n_data, default_stream)
    if (ierr /= cudaSuccess) then
      print *, 'CUDA error: ', trim(cudaGetErrorString(ierr))
      error stop 'copy_halo: device to device copy failed.'
    end if

  end subroutine copy_halo

end module m_cuda_sendrecv
