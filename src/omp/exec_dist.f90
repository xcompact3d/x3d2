module m_omp_exec_dist
   use mpi

   use m_common, only: dp
   use m_omp_common, only: SZ
   use m_omp_kernels_dist, only: der_univ_dist, der_univ_subs
   use m_tdsops, only: tdsops_t
   use m_omp_sendrecv, only: sendrecv_fields

   implicit none

contains

   subroutine exec_dist_tds_compact( &
      du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
      tdsops, nproc, pprev, pnext, n_block &
      )
      implicit none

      ! du = d(u)
      real(dp), dimension(:, :, :), intent(out) :: du
      real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e

      ! The ones below are intent(out) just so that we can write data in them,
      ! not because we actually need the data they store later where this
      ! subroutine is called. We absolutely don't care about the data they pass back
      real(dp), dimension(:, :, :), intent(out) :: &
         du_send_s, du_send_e, du_recv_s, du_recv_e

      type(tdsops_t), intent(in) :: tdsops
      integer, intent(in) :: nproc, pprev, pnext
      integer, intent(in) :: n_block

      integer :: n_data
      integer :: k

      n_data = SZ*n_block

      !$omp parallel do
      do k = 1, n_block
         call der_univ_dist( &
            du(:, :, k), du_send_s(:, :, k), du_send_e(:, :, k), u(:, :, k), &
            u_recv_s(:, :, k), u_recv_e(:, :, k), &
            tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, tdsops%n, &
            tdsops%dist_fw, tdsops%dist_bw, tdsops%dist_af &
            )
      end do
      !$omp end parallel do

      ! halo exchange for 2x2 systems
      call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                           n_data, nproc, pprev, pnext)

      !$omp parallel do
      do k = 1, n_block
         call der_univ_subs(du(:, :, k), &
                            du_recv_s(:, :, k), du_recv_e(:, :, k), &
                            tdsops%n, tdsops%dist_sa, tdsops%dist_sc)
      end do
      !$omp end parallel do

   end subroutine exec_dist_tds_compact

end module m_omp_exec_dist

