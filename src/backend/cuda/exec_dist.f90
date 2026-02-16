module m_cuda_exec_dist
  !! Distributed compact scheme execution on GPU.
  !!
  !! Orchestrates CUDA kernel launches and MPI halo exchange for distributed
  !! compact finite difference schemes. Handles both generic derivative operations
  !! and fused transport equation computation.
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
    !! Execute distributed compact scheme derivative $du = d(u)$ on GPU.
    !!
    !! Calls distributed kernel, exchanges halo data for $2 \times 2$ boundary
    !! systems, then applies substitution kernel.
    implicit none

    ! du = d(u)

    ! The ones below are intent(out) just so that we can write data in them,
    ! not because we actually need the data they store later where this
    ! subroutine is called. We absolutely don't care the data they pass back
    real(dp), device, dimension(:, :, :), intent(out) :: du  !! Output: derivative
    real(dp), device, dimension(:, :, :), intent(in) :: u  !! Input: field with local data
    real(dp), device, dimension(:, :, :), intent(in) :: u_recv_s, u_recv_e  !! Halo data from neighbours

    ! Temporary buffers for halo exchange (overwritten during computation)
    real(dp), device, dimension(:, :, :), intent(out) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e

    type(cuda_tdsops_t), intent(in) :: tdsops  !! Tridiagonal operators
    integer, intent(in) :: nproc, pprev, pnext  !! MPI ranks (total, previous, next)
    type(dim3), intent(in) :: blocks, threads  !! CUDA kernel configuration

    integer :: n_data

    n_data = SZ*1*blocks%x

    call der_univ_dist<<<blocks, threads>>>( & !&
      du, du_send_s, du_send_e, u, u_recv_s, u_recv_e, &
      tdsops%n_tds, tdsops%n_rhs, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%dist_fw_dev, tdsops%dist_bw_dev, tdsops%dist_af_dev &
      )

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, nproc, pprev, pnext)

    call der_univ_subs<<<blocks, threads>>>( & !&
      du, du_recv_s, du_recv_e, tdsops%n_tds, &
      tdsops%dist_sa_dev, tdsops%dist_sc_dev, tdsops%stretch_dev &
      )

  end subroutine exec_dist_tds_compact

  subroutine exec_dist_transeq_3fused( &
    r_du, u, u_recv_s, u_recv_e, v, v_recv_s, v_recv_e, &
    dud, d2u, &
    du_send_s, du_send_e, du_recv_s, du_recv_e, &
    dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
    d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
    tdsops_du, tdsops_dud, tdsops_d2u, nu, nproc, pprev, pnext, &
    blocks, threads &
    )
    !! Execute fused transport equation computation on GPU with distributed compact scheme.
    !!
    !! Computes $r\_du = -\frac{1}{2}(v \frac{\partial u}{\partial x} + \frac{\partial (uv)}{\partial x}) + \nu \frac{\partial^2 u}{\partial x^2}$
    !! Launches distributed kernel for three operators (du, dud, d2u), exchanges halo data for all
    !! boundary systems in one batch, then applies substitution kernel.
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: r_du  !! Output: transport equation RHS
    real(dp), device, dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e  !! Field u with halos
    real(dp), device, dimension(:, :, :), intent(in) :: v, v_recv_s, v_recv_e  !! Field v with halos

    ! Temporary storage for derivatives and halo exchange buffers
    real(dp), device, dimension(:, :, :), intent(out) :: dud, d2u
    real(dp), device, dimension(:, :, :), intent(out) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e, &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e

    type(cuda_tdsops_t), intent(in) :: tdsops_du, tdsops_dud, tdsops_d2u  !! Operators for each derivative
    real(dp), intent(in) :: nu  !! Kinematic viscosity
    integer, intent(in) :: nproc, pprev, pnext  !! MPI ranks
    type(dim3), intent(in) :: blocks, threads  !! CUDA kernel configuration

    integer :: n_data

    n_data = SZ*1*blocks%x

    call transeq_3fused_dist<<<blocks, threads>>>( & !&
      r_du, dud, d2u, &
      du_send_s, du_send_e, &
      dud_send_s, dud_send_e, &
      d2u_send_s, d2u_send_e, &
      u, u_recv_s, u_recv_e, &
      v, v_recv_s, v_recv_e, &
      tdsops_du%n_tds, tdsops_du%n_rhs, &
      tdsops_du%coeffs_s_dev, tdsops_du%coeffs_e_dev, tdsops_du%coeffs_dev, &
      tdsops_du%dist_fw_dev, tdsops_du%dist_bw_dev, tdsops_du%dist_af_dev, &
      tdsops_dud%coeffs_s_dev, tdsops_dud%coeffs_e_dev, &
      tdsops_dud%coeffs_dev, tdsops_dud%dist_fw_dev, tdsops_dud%dist_bw_dev, &
      tdsops_dud%dist_af_dev, &
      tdsops_d2u%coeffs_s_dev, tdsops_d2u%coeffs_e_dev, &
      tdsops_d2u%coeffs_dev, tdsops_d2u%dist_fw_dev, tdsops_d2u%dist_bw_dev, &
      tdsops_d2u%dist_af_dev &
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
      du_recv_s, du_recv_e, dud_recv_s, dud_recv_e, d2u_recv_s, d2u_recv_e, &
      tdsops_du%n_tds, nu, &
      tdsops_du%dist_sa_dev, tdsops_du%dist_sc_dev, tdsops_du%stretch_dev, &
      tdsops_dud%dist_sa_dev, tdsops_dud%dist_sc_dev, tdsops_dud%stretch_dev, &
      tdsops_d2u%dist_sa_dev, tdsops_d2u%dist_sc_dev, tdsops_d2u%stretch_dev, &
      tdsops_d2u%stretch_correct_dev &
      )

  end subroutine exec_dist_transeq_3fused

end module m_cuda_exec_dist
