module m_cuda_exec_dist
  use cudafor
  use mpi

  use m_common, only: dp
  use m_cuda_common, only: SZ
  use m_cuda_kernels_dist, only: der_univ_dist, der_univ_subs, &
                                 transeq_3fused_dist, transeq_3fused_subs
  use m_cuda_sendrecv, only: sendrecv_fields, sendrecv_3fields
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

    integer :: n_data

    n_data = SZ*1*blocks%x

    call der_univ_dist<<<blocks, threads>>>( & !&
      du, du_send_s, du_send_e, u, u_recv_s, u_recv_e, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%tds_n, tdsops%dist_fw_dev, tdsops%dist_bw_dev, &
      tdsops%dist_af_dev &
      )

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, nproc, pprev, pnext)

    call der_univ_subs<<<blocks, threads>>>( & !&
      du, du_recv_s, du_recv_e, &
      tdsops%tds_n, tdsops%dist_sa_dev, tdsops%dist_sc_dev &
      )

  end subroutine exec_dist_tds_compact

  subroutine exec_dist_transeq_3fused( &
    r_du, u, u_recv_s, u_recv_e, v, v_recv_s, v_recv_e, &
    dud, d2u, &
    du_send_s, du_send_e, du_recv_s, du_recv_e, &
    dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
    d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
    der1st, der2nd, nu, nproc, pprev, pnext, blocks, threads &
    )
    implicit none

    ! r_du = -1/2*(v*d1(u) + d1(u*v)) + nu*d2(u)
    !> The result array, it is also used as temporary storage
    real(dp), device, dimension(:, :, :), intent(out) :: r_du
    real(dp), device, dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e
    real(dp), device, dimension(:, :, :), intent(in) :: v, v_recv_s, v_recv_e

    ! The ones below are intent(out) just so that we can write data in them,
    ! not because we actually need the data they store later where this
    ! subroutine is called. We absolutely don't care the data they pass back
    real(dp), device, dimension(:, :, :), intent(out) :: dud, d2u
    real(dp), device, dimension(:, :, :), intent(out) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e, &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e

    type(cuda_tdsops_t), intent(in) :: der1st, der2nd
    real(dp), intent(in) :: nu
    integer, intent(in) :: nproc, pprev, pnext
    type(dim3), intent(in) :: blocks, threads

    integer :: n_data

    n_data = SZ*1*blocks%x

    call transeq_3fused_dist<<<blocks, threads>>>( & !&
      r_du, dud, d2u, &
      du_send_s, du_send_e, &
      dud_send_s, dud_send_e, &
      d2u_send_s, d2u_send_e, &
      u, u_recv_s, u_recv_e, &
      v, v_recv_s, v_recv_e, der1st%tds_n, &
      der1st%coeffs_s_dev, der1st%coeffs_e_dev, der1st%coeffs_dev, &
      der1st%dist_fw_dev, der1st%dist_bw_dev, der1st%dist_af_dev, &
      der2nd%coeffs_s_dev, der2nd%coeffs_e_dev, der2nd%coeffs_dev, &
      der2nd%dist_fw_dev, der2nd%dist_bw_dev, der2nd%dist_af_dev &
      )

    ! halo exchange for 2x2 systems
    call sendrecv_3fields( &
      du_recv_s, du_recv_e, dud_recv_s, dud_recv_e, &
      d2u_recv_s, d2u_recv_e, &
      du_send_s, du_send_e, dud_send_s, dud_send_e, &
      d2u_send_s, d2u_send_e, &
      n_data, nproc, pprev, pnext &
      )

    call transeq_3fused_subs<<<blocks, threads>>>( & !&
      r_du, v, dud, d2u, &
      du_recv_s, du_recv_e, &
      dud_recv_s, dud_recv_e, &
      d2u_recv_s, d2u_recv_e, &
      der1st%dist_sa_dev, der1st%dist_sc_dev, &
      der2nd%dist_sa_dev, der2nd%dist_sc_dev, &
      der1st%tds_n, nu &
      )

  end subroutine exec_dist_transeq_3fused

end module m_cuda_exec_dist
