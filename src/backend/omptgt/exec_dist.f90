module m_omptgt_exec_dist
  !! Drivers of the distributed compact solves for the OpenMP target backend,
  !! following cuda/exec_dist.f90. Fields are (SZ, n, n_groups) device arrays
  !! and only the first n_groups groups are processed.
  use iso_c_binding, only: c_loc
  use m_logger, only: abort_run
  use m_common, only: dp
  use m_dist_reduced, only: solve_reduced_exact
  use m_mesh_content, only: par_t
  use m_omptgt_common, only: SZ, device_sync, memcpy_d2h, memcpy_h2d
  use m_omptgt_kernels_dist, only: der_univ_dist, der_univ_subs, &
                                   der_univ_subs_exact, &
                                   der_penta_full, der_penta_periodic, &
                                   transeq_3fused_dist, transeq_3fused_subs
  use m_omptgt_sendrecv, only: sendrecv_fields, sendrecv_3fields
  use m_omptgt_tdsops, only: omptgt_tdsops_t

  implicit none

contains

  subroutine exec_dist_tds_compact( &
    du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
    tdsops, nproc, pprev, pnext, n_groups &
    )
    implicit none

    ! du = d(u)
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du, u, u_recv_s, u_recv_e
    ! Scratch for the 2x2 boundary systems
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e
    type(omptgt_tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: nproc, pprev, pnext, n_groups

    integer :: n_data

    n_data = SZ*1*n_groups

    call der_univ_dist( &
      du, du_send_s, du_send_e, u, u_recv_s, u_recv_e, &
      tdsops%n_tds, tdsops%n_rhs, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%dist_fw_dev, tdsops%dist_bw_dev, tdsops%dist_af_dev, &
      size(u, 2), n_groups &
      )

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, nproc, pprev, pnext)

    call der_univ_subs( &
      du, du_recv_s, du_recv_e, tdsops%n_tds, &
      tdsops%dist_sa_dev, tdsops%dist_sc_dev, tdsops%stretch_dev, &
      size(u, 2), n_groups &
      )

  end subroutine exec_dist_tds_compact

  subroutine exec_dist_tds_compact_exact( &
    du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
    tdsops, par, dir, n_groups &
    )
    !! As exec_dist_tds_compact, but with the reduced system that couples the
    !! ranks solved exactly (see m_dist_reduced and tdsops%dist_exact). Its
    !! right-hand sides, two planes per field, go through the host.
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du, u, u_recv_s, u_recv_e
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e
    type(omptgt_tdsops_t), intent(in) :: tdsops
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir, n_groups

    real(dp), allocatable, target :: d_s(:), d_e(:), x_s(:), x_e(:)
    integer :: n, n_data

    n = tdsops%n_tds
    n_data = SZ*n_groups

    call der_univ_dist( &
      du, du_send_s, du_send_e, u, u_recv_s, u_recv_e, &
      tdsops%n_tds, tdsops%n_rhs, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%dist_fw_dev, tdsops%dist_bw_dev, tdsops%dist_af_dev, &
      size(u, 2), n_groups &
      )

    ! The kernel must be complete before its results are copied back
    call device_sync()
    allocate (d_s(n_data), d_e(n_data), x_s(n_data), x_e(n_data))
    call memcpy_d2h(c_loc(d_s), c_loc(du_send_s), n_data)
    call memcpy_d2h(c_loc(d_e), c_loc(du_send_e), n_data)
    call solve_reduced_exact( &
      x_s, x_e, d_s, d_e, tdsops%dist_sa(1), tdsops%dist_sc(1), &
      tdsops%dist_sa(n), tdsops%dist_sc(n), par, dir)
    call memcpy_h2d(c_loc(du_recv_s), c_loc(x_s), n_data)
    call memcpy_h2d(c_loc(du_recv_e), c_loc(x_e), n_data)

    call der_univ_subs_exact( &
      du, du_recv_s, du_recv_e, tdsops%n_tds, &
      tdsops%dist_sa_dev, tdsops%dist_sc_dev, tdsops%stretch_dev, &
      size(u, 2), n_groups &
      )

  end subroutine exec_dist_tds_compact_exact

  subroutine exec_dist_transeq_3fused( &
    r_du, u, u_recv_s, u_recv_e, v, v_recv_s, v_recv_e, &
    dud, d2u, &
    du_send_s, du_send_e, du_recv_s, du_recv_e, &
    dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
    d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
    tdsops_du, tdsops_dud, tdsops_d2u, nu, nproc, pprev, pnext, n_groups &
    )
    implicit none

    ! r_du = -1/2*(v*d1(u) + d1(u*v)) + nu*d2(u)
    !> The result array, it is also used as temporary storage
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: r_du
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      u, u_recv_s, u_recv_e, v, v_recv_s, v_recv_e
    ! Scratch fields and buffers
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      dud, d2u, &
      du_send_s, du_send_e, du_recv_s, du_recv_e, &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e
    type(omptgt_tdsops_t), intent(in) :: tdsops_du, tdsops_dud, tdsops_d2u
    real(dp), intent(in) :: nu
    integer, intent(in) :: nproc, pprev, pnext, n_groups

    integer :: n_data

    n_data = SZ*1*n_groups

    call transeq_3fused_dist( &
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
      tdsops_d2u%dist_af_dev, size(u, 2), n_groups &
      )

    ! halo exchange for 2x2 systems
    call sendrecv_3fields( &
      du_recv_s, du_recv_e, dud_recv_s, dud_recv_e, &
      d2u_recv_s, d2u_recv_e, &
      du_send_s, du_send_e, dud_send_s, dud_send_e, &
      d2u_send_s, d2u_send_e, &
      n_data, nproc, pprev, pnext &
      )

    call transeq_3fused_subs( &
      r_du, v, dud, d2u, &
      du_recv_s, du_recv_e, dud_recv_s, dud_recv_e, d2u_recv_s, d2u_recv_e, &
      tdsops_du%n_tds, nu, &
      tdsops_du%dist_sa_dev, tdsops_du%dist_sc_dev, tdsops_du%stretch_dev, &
      tdsops_dud%dist_sa_dev, tdsops_dud%dist_sc_dev, tdsops_dud%stretch_dev, &
      tdsops_d2u%dist_sa_dev, tdsops_d2u%dist_sc_dev, tdsops_d2u%stretch_dev, &
      tdsops_d2u%stretch_correct_dev, size(u, 2), n_groups &
      )

  end subroutine exec_dist_transeq_3fused

  subroutine exec_dist_penta_compact(du, u, u_recv_s, u_recv_e, tdsops, &
                                     n_groups)
    !! Single-subdomain non-periodic pentadiagonal compact solve.
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du, u, u_recv_s, u_recv_e
    type(omptgt_tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n_groups

    call der_penta_full( &
      du, u, u_recv_s, u_recv_e, &
      tdsops%n_tds, tdsops%n_rhs, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%dist_fw_dev, tdsops%dist_af_dev, tdsops%dist_sa_dev, &
      tdsops%dist_bw_dev, tdsops%beta, tdsops%beta_lhs_s, &
      size(u, 2), n_groups &
      )

  end subroutine exec_dist_penta_compact

  subroutine exec_dist_penta_periodic(du, u, u_recv_s, u_recv_e, tdsops, &
                                      n_groups)
    !! Cyclic pentadiagonal compact solve via Sherman-Morrison-Woodbury.
    !! Halos u_recv_s/e must hold the periodic extension.
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      du, u, u_recv_s, u_recv_e
    type(omptgt_tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n_groups

    if (.not. associated(tdsops%penta_z_dev)) then
      call abort_run('exec_dist_penta_periodic: not a periodic pentadiagonal &
                  &operator.')
    end if

    call der_penta_periodic( &
      du, u, u_recv_s, u_recv_e, &
      tdsops%n_tds, tdsops%n_rhs, &
      tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
      tdsops%dist_fw_dev, tdsops%dist_af_dev, tdsops%dist_sa_dev, &
      tdsops%dist_bw_dev, tdsops%beta, tdsops%beta_lhs_s, &
      tdsops%alpha, tdsops%beta, tdsops%penta_z_dev, size(u, 2), n_groups &
      )

  end subroutine exec_dist_penta_periodic

end module m_omptgt_exec_dist
