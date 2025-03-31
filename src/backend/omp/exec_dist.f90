module m_omp_exec_dist
  use mpi

  use m_common, only: dp
  use m_omp_common, only: SZ
  use m_omp_kernels_dist, only: der_univ_dist, der_univ_subs, &
                                der_univ_fused_subs
  use m_tdsops, only: tdsops_t
  use m_omp_sendrecv, only: sendrecv_fields

  implicit none

contains

  subroutine exec_dist_tds_compact( &
    du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
    tdsops, nproc, pprev, pnext, n_groups)
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
    integer, intent(in) :: n_groups

    integer :: n_data
    integer :: k

    n_data = SZ*n_groups

    !$omp parallel do
    do k = 1, n_groups
      call der_univ_dist( &
        du(:, :, k), du_send_s(:, :, k), du_send_e(:, :, k), &
        u(:, :, k), u_recv_s(:, :, k), u_recv_e(:, :, k), &
        tdsops%n_tds, tdsops%n_rhs, &
        tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, &
        tdsops%dist_fw, tdsops%dist_bw, tdsops%dist_af &
        )
    end do
    !$omp end parallel do

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, nproc, pprev, pnext)

    !$omp parallel do
    do k = 1, n_groups
      call der_univ_subs( &
        du(:, :, k), du_recv_s(:, :, k), du_recv_e(:, :, k), &
        tdsops%n_tds, tdsops%dist_sa, tdsops%dist_sc, tdsops%stretch &
        )
    end do
    !$omp end parallel do

  end subroutine exec_dist_tds_compact

  subroutine exec_dist_transeq_compact( &
    rhs_du, dud, d2u, &
    du_send_s, du_send_e, du_recv_s, du_recv_e, &
    dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
    d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
    u, u_recv_s, u_recv_e, &
    v, v_recv_s, v_recv_e, &
    tdsops_du, tdsops_dud, tdsops_d2u, nu, nproc, pprev, pnext, n_groups)

    implicit none

    !> The result array, it is also used as temporary storage
    real(dp), dimension(:, :, :), intent(out) :: rhs_du
    !> Temporary storage arrays
    real(dp), dimension(:, :, :), intent(out) :: dud, d2u

    ! The ones below are intent(out) just so that we can write data in them,
    ! not because we actually need the data they store later where this
    ! subroutine is called. We absolutely don't care about the data they pass back
    real(dp), dimension(:, :, :), intent(out) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e
    real(dp), dimension(:, :, :), intent(out) :: &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e
    real(dp), dimension(:, :, :), intent(out) :: &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e

    real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e
    real(dp), dimension(:, :, :), intent(in) :: v, v_recv_s, v_recv_e

    type(tdsops_t), intent(in) :: tdsops_du, tdsops_dud, tdsops_d2u
    real(dp), intent(in) :: nu
    integer, intent(in) :: nproc, pprev, pnext
    integer, intent(in) :: n_groups

    real(dp), dimension(:, :), allocatable :: ud, ud_recv_s, ud_recv_e

    integer :: n_data, n_halo
    integer :: i, j, k

    ! TODO: don't hardcode n_halo
    n_halo = 4
    n_data = SZ*n_groups

    allocate (ud(SZ, tdsops_dud%n_tds))
    allocate (ud_recv_e(SZ, n_halo))
    allocate (ud_recv_s(SZ, n_halo))

    !$omp parallel do private(ud, ud_recv_e, ud_recv_s)
    do k = 1, n_groups
      call der_univ_dist( &
        rhs_du(:, :, k), du_send_s(:, :, k), du_send_e(:, :, k), &
        u(:, :, k), u_recv_s(:, :, k), u_recv_e(:, :, k), &
        tdsops_du%n_tds, tdsops_du%n_rhs, &
        tdsops_du%coeffs_s, tdsops_du%coeffs_e, tdsops_du%coeffs, &
        tdsops_du%dist_fw, tdsops_du%dist_bw, tdsops_du%dist_af &
        )

      call der_univ_dist( &
        d2u(:, :, k), d2u_send_s(:, :, k), d2u_send_e(:, :, k), &
        u(:, :, k), u_recv_s(:, :, k), u_recv_e(:, :, k), &
        tdsops_d2u%n_tds, tdsops_d2u%n_rhs, &
        tdsops_d2u%coeffs_s, tdsops_d2u%coeffs_e, tdsops_d2u%coeffs, &
        tdsops_d2u%dist_fw, tdsops_d2u%dist_bw, tdsops_d2u%dist_af &
        )

      ! Handle dud by locally generating u*v
      do j = 1, tdsops_dud%n_tds
        !$omp simd
        do i = 1, SZ
          ud(i, j) = u(i, j, k)*v(i, j, k)
        end do
        !$omp end simd
      end do

      do j = 1, n_halo
        !$omp simd
        do i = 1, SZ
          ud_recv_s(i, j) = u_recv_s(i, j, k)*v_recv_s(i, j, k)
          ud_recv_e(i, j) = u_recv_e(i, j, k)*v_recv_e(i, j, k)
        end do
        !$omp end simd
      end do

      call der_univ_dist( &
        dud(:, :, k), dud_send_s(:, :, k), dud_send_e(:, :, k), &
        ud(:, :), ud_recv_s(:, :), ud_recv_e(:, :), &
        tdsops_dud%n_tds, tdsops_dud%n_rhs, &
        tdsops_dud%coeffs_s, tdsops_dud%coeffs_e, tdsops_dud%coeffs, &
        tdsops_dud%dist_fw, tdsops_dud%dist_bw, tdsops_dud%dist_af &
        )

    end do
    !$omp end parallel do

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, nproc, pprev, pnext)
    call sendrecv_fields(dud_recv_s, dud_recv_e, dud_send_s, dud_send_e, &
                         n_data, nproc, pprev, pnext)
    call sendrecv_fields(d2u_recv_s, d2u_recv_e, d2u_send_s, d2u_send_e, &
                         n_data, nproc, pprev, pnext)

    !$omp parallel do
    do k = 1, n_groups
      call der_univ_fused_subs( &
        rhs_du(:, :, k), dud(:, :, k), d2u(:, :, k), v(:, :, k), &
        du_recv_s(:, :, k), du_recv_e(:, :, k), &
        dud_recv_s(:, :, k), dud_recv_e(:, :, k), &
        d2u_recv_s(:, :, k), d2u_recv_e(:, :, k), &
        nu, tdsops_du%n_tds, &
        tdsops_du%dist_sa, tdsops_du%dist_sc, tdsops_du%stretch, &
        tdsops_dud%dist_sa, tdsops_dud%dist_sc, tdsops_dud%stretch, &
        tdsops_d2u%dist_sa, tdsops_d2u%dist_sc, tdsops_d2u%stretch, &
        tdsops_d2u%stretch_correct &
        )
    end do
    !$omp end parallel do

  end subroutine exec_dist_transeq_compact

end module m_omp_exec_dist

