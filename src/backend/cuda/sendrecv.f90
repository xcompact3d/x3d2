module m_cuda_sendrecv
  use cudafor
  use m_mpi, only: MPI_COMM_WORLD, MPI_STATUSES_IGNORE, MPI_Irecv, MPI_Isend, &
                   MPI_Waitall

  use m_common, only: dp, MPI_X3D2_DP
  use m_logger, only: check_mpi, abort_run
  use m_cuda_errors, only: check_cuda

  implicit none

  integer(kind=cuda_stream_kind), parameter :: default_stream = 0

contains

  subroutine sendrecv_fields(f_recv_s, f_recv_e, f_send_s, f_send_e, &
                             n_data, nproc, prev, next)
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: f_recv_s, f_recv_e
    real(dp), device, dimension(:, :, :), intent(in) :: f_send_s, f_send_e
    integer, intent(in) :: n_data, nproc, prev, next

    integer :: req(4), err(4), ierr, i, tag = 1234

    if (nproc == 1) then
      call copy_halo(f_recv_s, f_send_e, n_data)
      call copy_halo(f_recv_e, f_send_s, n_data)
    else
      ! Ensure the send buffers are ready before MPI touches device memory.
      ierr = cudaStreamSynchronize(default_stream)
      call check_cuda(ierr, 'cudaStreamSynchronize in sendrecv_fields')

#ifdef MPI
      call MPI_Isend(f_send_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(1), err(1))
      call MPI_Irecv(f_recv_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(2), err(2))
      call MPI_Isend(f_send_e, n_data, MPI_X3D2_DP, &
                     next, tag, MPI_COMM_WORLD, req(3), err(3))
      call MPI_Irecv(f_recv_s, n_data, MPI_X3D2_DP, &
                     prev, tag, MPI_COMM_WORLD, req(4), err(4))

      do i = 1, 4
        call check_mpi(err(i), 'MPI_Isend/MPI_Irecv in sendrecv_fields')
      end do

      call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)
      call check_mpi(ierr, 'MPI_Waitall in sendrecv_fields')
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

    integer :: req(12), err(12), ierr, i, tag = 1234

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
      call check_cuda(ierr, 'cudaStreamSynchronize in sendrecv_3fields')

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

      do i = 1, 12
        call check_mpi(err(i), 'MPI_Isend/MPI_Irecv in sendrecv_3fields')
      end do

      call MPI_Waitall(12, req, MPI_STATUSES_IGNORE, ierr)
      call check_mpi(ierr, 'MPI_Waitall in sendrecv_3fields')
#endif
    end if

  end subroutine sendrecv_3fields

  subroutine copy_halo(dst, src, n_data)
    !! Periodic wrap-around on a single rank: copies the first n_data
    !! elements of a halo buffer, the groups of the direction being
    !! exchanged. The buffers are sized for the direction with the most
    !! groups, so a whole-array assignment would copy that maximum on every
    !! exchange (16 times the data in use on an elongated mesh).
    !!
    !! The copy is a device array section assignment, which CUDA Fortran runs
    !! on the device in the default stream like the kernels around it. The
    !! dummies are plain assumed shape: declaring them contiguous made the
    !! caller copy its arguments through a host temporary in a debug build,
    !! touching device memory from the host.
    implicit none

    real(dp), device, dimension(:, :, :), intent(inout) :: dst
    real(dp), device, dimension(:, :, :), intent(in) :: src
    integer, intent(in) :: n_data

    integer :: n_groups

    ! n_data is SZ*n_halo (or SZ for the 2x2 systems) per group
    n_groups = n_data/(size(src, 1)*size(src, 2))
    if (n_groups*size(src, 1)*size(src, 2) /= n_data .or. &
        n_groups > size(src, 3) .or. n_groups > size(dst, 3)) then
      call abort_run('copy_halo: n_data does not match the halo buffer.')
    end if

    dst(:, :, 1:n_groups) = src(:, :, 1:n_groups)

  end subroutine copy_halo

end module m_cuda_sendrecv
