module m_cuda_exec_dist
   use cudafor
   use mpi

   use m_common, only: dp
   use m_cuda_common, only: SZ
   use m_cuda_kernels_dist, only: der_univ_dist, der_univ_subs
   use m_cuda_tdsops, only: cuda_tdsops_t

   implicit none

contains

   subroutine exec_dist_tds_compact( &
      du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
      tdsops, nproc, pprev, pnext, blocks, threads &
   )
      implicit none

      ! du = d(u)
      real(dp), device, dimension(:, :, :), intent(out) :: du
      real(dp), device, dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e

      ! The ones below are intent(out) just so that we can write data in them,
      ! not because we actually need the data they store later where this
      ! subroutine is called. We absolutely don't care the data they pass back
      real(dp), device, dimension(:, :, :), intent(out) :: &
         du_send_s, du_send_e, du_recv_s, du_recv_e

      type(cuda_tdsops_t), intent(in) :: tdsops
      integer, intent(in) :: nproc, pprev, pnext
      type(dim3), intent(in) :: blocks, threads

      integer :: n_block
      integer :: ierr, tag = 1234
      integer :: srerr(4), mpireq(4)

      n_block = blocks%x

      call der_univ_dist<<<blocks, threads>>>( &
         du, du_send_s, du_send_e, u, u_recv_s, u_recv_e, &
         tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
         tdsops%n, tdsops%dist_fw_dev, tdsops%dist_bw_dev, tdsops%dist_af_dev &
      )

      ! halo exchange for 2x2 systems
      if (nproc == 1) then
         du_recv_s = du_send_e
         du_recv_e = du_send_s
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(du_send_s, SZ*1*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(du_recv_e, SZ*1*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(du_send_e, SZ*1*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(du_recv_s, SZ*1*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if

      call der_univ_subs<<<blocks, threads>>>( &
         du, du_recv_s, du_recv_e, &
         tdsops%n, tdsops%dist_sa_dev, tdsops%dist_sc_dev &
      )

   end subroutine exec_dist_tds_compact

end module m_cuda_exec_dist
