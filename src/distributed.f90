module m_cuda_kernels_dist
  use cudafor

  use m_common, only: dp

  implicit none

contains

  attributes(global) subroutine der_univ_dist( &
    du, send_u_s, send_u_e, u, u_s, u_e, &
    n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, ffr, fbc, faf &
    )
    implicit none

    ! Arguments
    real(dp), device, intent(out), dimension(:, :, :) :: du, send_u_s, &
                                                         send_u_e
    real(dp), device, intent(in), dimension(:, :, :) :: u, u_s, u_e
    integer, value, intent(in) :: n_tds, n_rhs
    real(dp), device, intent(in), dimension(:, :) :: coeffs_s, coeffs_e
    real(dp), device, intent(in), dimension(:) :: coeffs
    real(dp), device, intent(in), dimension(:) :: ffr, fbc, faf

    ! Local variables
    integer :: i, j, b, k, lj
    integer :: jm2, jm1, jp1, jp2

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, &
                temp_du, alpha, last_r

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)
    last_r = ffr(1)

    du(i, 1, b) = coeffs_s(1, 1)*u_s(i, 1, b) &
                  + coeffs_s(2, 1)*u_s(i, 2, b) &
                  + coeffs_s(3, 1)*u_s(i, 3, b) &
                  + coeffs_s(4, 1)*u_s(i, 4, b) &
                  + coeffs_s(5, 1)*u(i, 1, b) &
                  + coeffs_s(6, 1)*u(i, 2, b) &
                  + coeffs_s(7, 1)*u(i, 3, b) &
                  + coeffs_s(8, 1)*u(i, 4, b) &
                  + coeffs_s(9, 1)*u(i, 5, b)
    du(i, 1, b) = du(i, 1, b)*faf(1)
    du(i, 2, b) = coeffs_s(1, 2)*u_s(i, 2, b) &
                  + coeffs_s(2, 2)*u_s(i, 3, b) &
                  + coeffs_s(3, 2)*u_s(i, 4, b) &
                  + coeffs_s(4, 2)*u(i, 1, b) &
                  + coeffs_s(5, 2)*u(i, 2, b) &
                  + coeffs_s(6, 2)*u(i, 3, b) &
                  + coeffs_s(7, 2)*u(i, 4, b) &
                  + coeffs_s(8, 2)*u(i, 5, b) &
                  + coeffs_s(9, 2)*u(i, 6, b)
    du(i, 2, b) = du(i, 2, b)*faf(2)
    du(i, 3, b) = coeffs_s(1, 3)*u_s(i, 3, b) &
                  + coeffs_s(2, 3)*u_s(i, 4, b) &
                  + coeffs_s(3, 3)*u(i, 1, b) &
                  + coeffs_s(4, 3)*u(i, 2, b) &
                  + coeffs_s(5, 3)*u(i, 3, b) &
                  + coeffs_s(6, 3)*u(i, 4, b) &
                  + coeffs_s(7, 3)*u(i, 5, b) &
                  + coeffs_s(8, 3)*u(i, 6, b) &
                  + coeffs_s(9, 3)*u(i, 7, b)
    du(i, 3, b) = ffr(3)*(du(i, 3, b) - faf(3)*du(i, 2, b))
    du(i, 4, b) = coeffs_s(1, 4)*u_s(i, 4, b) &
                  + coeffs_s(2, 4)*u(i, 1, b) &
                  + coeffs_s(3, 4)*u(i, 2, b) &
                  + coeffs_s(4, 4)*u(i, 3, b) &
                  + coeffs_s(5, 4)*u(i, 4, b) &
                  + coeffs_s(6, 4)*u(i, 5, b) &
                  + coeffs_s(7, 4)*u(i, 6, b) &
                  + coeffs_s(8, 4)*u(i, 7, b) &
                  + coeffs_s(9, 4)*u(i, 8, b)
    du(i, 4, b) = ffr(4)*(du(i, 4, b) - faf(3)*du(i, 3, b))

    alpha = faf(5)

    do j = 5, n_rhs - 4
      temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
      du(i, j, b) = ffr(j)*(temp_du - alpha*du(i, j - 1, b))
    end do

    j = n_rhs - 3
    du(i, j, b) = coeffs_e(1, 1)*u(i, j - 4, b) &
                  + coeffs_e(2, 1)*u(i, j - 3, b) &
                  + coeffs_e(3, 1)*u(i, j - 2, b) &
                  + coeffs_e(4, 1)*u(i, j - 1, b) &
                  + coeffs_e(5, 1)*u(i, j, b) &
                  + coeffs_e(6, 1)*u(i, j + 1, b) &
                  + coeffs_e(7, 1)*u(i, j + 2, b) &
                  + coeffs_e(8, 1)*u(i, j + 3, b) &
                  + coeffs_e(9, 1)*u_e(i, 1, b)
    du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
    j = n_rhs - 2
    du(i, j, b) = coeffs_e(1, 2)*u(i, j - 4, b) &
                  + coeffs_e(2, 2)*u(i, j - 3, b) &
                  + coeffs_e(3, 2)*u(i, j - 2, b) &
                  + coeffs_e(4, 2)*u(i, j - 1, b) &
                  + coeffs_e(5, 2)*u(i, j, b) &
                  + coeffs_e(6, 2)*u(i, j + 1, b) &
                  + coeffs_e(7, 2)*u(i, j + 2, b) &
                  + coeffs_e(8, 2)*u_e(i, 1, b) &
                  + coeffs_e(9, 2)*u_e(i, 2, b)
    du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
    j = n_rhs - 1
    du(i, j, b) = coeffs_e(1, 3)*u(i, j - 4, b) &
                  + coeffs_e(2, 3)*u(i, j - 3, b) &
                  + coeffs_e(3, 3)*u(i, j - 2, b) &
                  + coeffs_e(4, 3)*u(i, j - 1, b) &
                  + coeffs_e(5, 3)*u(i, j, b) &
                  + coeffs_e(6, 3)*u(i, j + 1, b) &
                  + coeffs_e(7, 3)*u_e(i, 1, b) &
                  + coeffs_e(8, 3)*u_e(i, 2, b) &
                  + coeffs_e(9, 3)*u_e(i, 3, b)
    du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
    j = n_rhs
    du(i, j, b) = coeffs_e(1, 4)*u(i, j - 4, b) &
                  + coeffs_e(2, 4)*u(i, j - 3, b) &
                  + coeffs_e(3, 4)*u(i, j - 2, b) &
                  + coeffs_e(4, 4)*u(i, j - 1, b) &
                  + coeffs_e(5, 4)*u(i, j, b) &
                  + coeffs_e(6, 4)*u_e(i, 1, b) &
                  + coeffs_e(7, 4)*u_e(i, 2, b) &
                  + coeffs_e(8, 4)*u_e(i, 3, b) &
                  + coeffs_e(9, 4)*u_e(i, 4, b)
    du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))

    send_u_e(i, 1, b) = du(i, n_tds, b)

    ! Backward pass of the hybrid algorithm
    do j = n_tds - 2, 2, -1
      du(i, j, b) = du(i, j, b) - fbc(j)*du(i, j + 1, b)
    end do
    du(i, 1, b) = last_r*(du(i, 1, b) - fbc(1)*du(i, 2, b))
    send_u_s(i, 1, b) = du(i, 1, b)

  end subroutine der_univ_dist

  attributes(global) subroutine der_univ_subs(du, recv_u_s, recv_u_e, &
                                              n, dist_sa, dist_sc, strch)
    implicit none

    ! Arguments
    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: recv_u_s, recv_u_e
    real(dp), device, intent(in), dimension(:) :: dist_sa, dist_sc, strch
    integer, value, intent(in) :: n

    ! Local variables
    integer :: i, j, b
    real(dp) :: ur, bl, recp, du_s, du_e

    i = threadIdx%x
    b = blockIdx%x

    ! A small trick we do here is valid for symmetric Toeplitz matrices.
    ! In our case our matrices satisfy this criteria in the (5:n-4) region
    ! and as long as a rank has around at least 20 entries the assumptions
    ! we make here are perfectly valid.

    ! bl is the bottom left entry in the 2x2 matrix
    ! ur is the upper right entry in the 2x2 matrix

    ! Start
    ! At the start we have the 'bl', and assume 'ur'
    bl = dist_sa(1)
    ur = dist_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    du_s = recp*(du(i, 1, b) - bl*recv_u_s(i, 1, b))

    ! End
    ! At the end we have the 'ur', and assume 'bl'
    bl = dist_sc(n)
    ur = dist_sc(n)
    recp = 1._dp/(1._dp - ur*bl)
    du_e = recp*(du(i, n, b) - ur*recv_u_e(i, 1, b))

    du(i, 1, b) = du_s*strch(1)
    do j = 2, n - 1
      du(i, j, b) = (du(i, j, b) - dist_sa(j)*du_s - dist_sc(j)*du_e)*strch(j)
    end do
    du(i, n, b) = du_e*strch(n)

  end subroutine der_univ_subs

  attributes(global) subroutine transeq_3fused_dist( &
    du, dud, d2u, &
    send_du_s, send_du_e, send_dud_s, send_dud_e, send_d2u_s, send_d2u_e, &
    u, u_s, u_e, v, v_s, v_e, n_tds, n_rhs, &
    du_coeffs_s, du_coeffs_e, du_coeffs, du_fw, du_bw, du_af, &
    dud_coeffs_s, dud_coeffs_e, dud_coeffs, dud_fw, dud_bw, dud_af, &
    d2u_coeffs_s, d2u_coeffs_e, d2u_coeffs, d2u_fw, d2u_bw, d2u_af &
    )
    implicit none

    ! Arguments
    real(dp), device, intent(out), dimension(:, :, :) :: du, dud, d2u
    real(dp), device, intent(out), dimension(:, :, :) :: &
      send_du_s, send_du_e, send_dud_s, send_dud_e, send_d2u_s, send_d2u_e
    real(dp), device, intent(in), dimension(:, :, :) :: u, u_s, u_e, &
                                                        v, v_s, v_e
    integer, value, intent(in) :: n_tds, n_rhs
    real(dp), device, intent(in) :: du_coeffs_s(:, :), du_coeffs_e(:, :), &
                                    du_coeffs(:)
    real(dp), device, intent(in) :: du_fw(:), du_bw(:), du_af(:)
    real(dp), device, intent(in) :: dud_coeffs_s(:, :), dud_coeffs_e(:, :), &
                                    dud_coeffs(:)
    real(dp), device, intent(in) :: dud_fw(:), dud_bw(:), dud_af(:)
    real(dp), device, intent(in) :: d2u_coeffs_s(:, :), d2u_coeffs_e(:, :), &
                                    d2u_coeffs(:)
    real(dp), device, intent(in) :: d2u_fw(:), d2u_bw(:), d2u_af(:)

    ! Local variables
    integer :: i, j, b

    real(dp) :: du_c_m4, du_c_m3, du_c_m2, du_c_m1, du_c_j, &
                du_c_p1, du_c_p2, du_c_p3, du_c_p4, &
                du_alpha, du_last_r
    real(dp) :: dud_c_m4, dud_c_m3, dud_c_m2, dud_c_m1, dud_c_j, &
                dud_c_p1, dud_c_p2, dud_c_p3, dud_c_p4, &
                dud_alpha, dud_last_r
    real(dp) :: d2u_c_m4, d2u_c_m3, d2u_c_m2, d2u_c_m1, d2u_c_j, &
                d2u_c_p1, d2u_c_p2, d2u_c_p3, d2u_c_p4, &
                d2u_alpha, d2u_last_r
    real(dp) :: temp_du, temp_dud, temp_d2u
    real(dp) :: u_m4, u_m3, u_m2, u_m1, u_j, u_p1, u_p2, u_p3, u_p4
    real(dp) :: v_m4, v_m3, v_m2, v_m1, v_j, v_p1, v_p2, v_p3, v_p4
    real(dp) :: old_du, old_dud, old_d2u

    i = threadIdx%x
    b = blockIdx%x

    du_last_r = du_fw(1)
    dud_last_r = dud_fw(1)
    d2u_last_r = d2u_fw(1)

    ! j = 1
    temp_du = du_coeffs_s(1, 1)*u_s(i, 1, b) &
              + du_coeffs_s(2, 1)*u_s(i, 2, b) &
              + du_coeffs_s(3, 1)*u_s(i, 3, b) &
              + du_coeffs_s(4, 1)*u_s(i, 4, b) &
              + du_coeffs_s(5, 1)*u(i, 1, b) &
              + du_coeffs_s(6, 1)*u(i, 2, b) &
              + du_coeffs_s(7, 1)*u(i, 3, b) &
              + du_coeffs_s(8, 1)*u(i, 4, b) &
              + du_coeffs_s(9, 1)*u(i, 5, b)
    du(i, 1, b) = temp_du*du_af(1)
    temp_dud = dud_coeffs_s(1, 1)*u_s(i, 1, b)*v_s(i, 1, b) &
               + dud_coeffs_s(2, 1)*u_s(i, 2, b)*v_s(i, 2, b) &
               + dud_coeffs_s(3, 1)*u_s(i, 3, b)*v_s(i, 3, b) &
               + dud_coeffs_s(4, 1)*u_s(i, 4, b)*v_s(i, 4, b) &
               + dud_coeffs_s(5, 1)*u(i, 1, b)*v(i, 1, b) &
               + dud_coeffs_s(6, 1)*u(i, 2, b)*v(i, 2, b) &
               + dud_coeffs_s(7, 1)*u(i, 3, b)*v(i, 3, b) &
               + dud_coeffs_s(8, 1)*u(i, 4, b)*v(i, 4, b) &
               + dud_coeffs_s(9, 1)*u(i, 5, b)*v(i, 5, b)
    dud(i, 1, b) = temp_dud*dud_af(1)
    temp_d2u = d2u_coeffs_s(1, 1)*u_s(i, 1, b) &
               + d2u_coeffs_s(2, 1)*u_s(i, 2, b) &
               + d2u_coeffs_s(3, 1)*u_s(i, 3, b) &
               + d2u_coeffs_s(4, 1)*u_s(i, 4, b) &
               + d2u_coeffs_s(5, 1)*u(i, 1, b) &
               + d2u_coeffs_s(6, 1)*u(i, 2, b) &
               + d2u_coeffs_s(7, 1)*u(i, 3, b) &
               + d2u_coeffs_s(8, 1)*u(i, 4, b) &
               + d2u_coeffs_s(9, 1)*u(i, 5, b)
    d2u(i, 1, b) = temp_d2u*d2u_af(1)
    ! j = 2
    temp_du = du_coeffs_s(1, 2)*u_s(i, 2, b) &
              + du_coeffs_s(2, 2)*u_s(i, 3, b) &
              + du_coeffs_s(3, 2)*u_s(i, 4, b) &
              + du_coeffs_s(4, 2)*u(i, 1, b) &
              + du_coeffs_s(5, 2)*u(i, 2, b) &
              + du_coeffs_s(6, 2)*u(i, 3, b) &
              + du_coeffs_s(7, 2)*u(i, 4, b) &
              + du_coeffs_s(8, 2)*u(i, 5, b) &
              + du_coeffs_s(9, 2)*u(i, 6, b)
    du(i, 2, b) = temp_du*du_af(2)
    temp_dud = dud_coeffs_s(1, 2)*u_s(i, 2, b)*v_s(i, 2, b) &
               + dud_coeffs_s(2, 2)*u_s(i, 3, b)*v_s(i, 3, b) &
               + dud_coeffs_s(3, 2)*u_s(i, 4, b)*v_s(i, 4, b) &
               + dud_coeffs_s(4, 2)*u(i, 1, b)*v(i, 1, b) &
               + dud_coeffs_s(5, 2)*u(i, 2, b)*v(i, 2, b) &
               + dud_coeffs_s(6, 2)*u(i, 3, b)*v(i, 3, b) &
               + dud_coeffs_s(7, 2)*u(i, 4, b)*v(i, 4, b) &
               + dud_coeffs_s(8, 2)*u(i, 5, b)*v(i, 5, b) &
               + dud_coeffs_s(9, 2)*u(i, 6, b)*v(i, 6, b)
    dud(i, 2, b) = temp_dud*dud_af(2)
    temp_d2u = d2u_coeffs_s(1, 2)*u_s(i, 2, b) &
               + d2u_coeffs_s(2, 2)*u_s(i, 3, b) &
               + d2u_coeffs_s(3, 2)*u_s(i, 4, b) &
               + d2u_coeffs_s(4, 2)*u(i, 1, b) &
               + d2u_coeffs_s(5, 2)*u(i, 2, b) &
               + d2u_coeffs_s(6, 2)*u(i, 3, b) &
               + d2u_coeffs_s(7, 2)*u(i, 4, b) &
               + d2u_coeffs_s(8, 2)*u(i, 5, b) &
               + d2u_coeffs_s(9, 2)*u(i, 6, b)
    d2u(i, 2, b) = temp_d2u*d2u_af(2)
    ! j = 3
    temp_du = du_coeffs_s(1, 3)*u_s(i, 3, b) &
              + du_coeffs_s(2, 3)*u_s(i, 4, b) &
              + du_coeffs_s(3, 3)*u(i, 1, b) &
              + du_coeffs_s(4, 3)*u(i, 2, b) &
              + du_coeffs_s(5, 3)*u(i, 3, b) &
              + du_coeffs_s(6, 3)*u(i, 4, b) &
              + du_coeffs_s(7, 3)*u(i, 5, b) &
              + du_coeffs_s(8, 3)*u(i, 6, b) &
              + du_coeffs_s(9, 3)*u(i, 7, b)
    du(i, 3, b) = du_fw(3)*(temp_du - du_af(3)*du(i, 2, b))
    temp_dud = dud_coeffs_s(1, 3)*u_s(i, 3, b)*v_s(i, 3, b) &
               + dud_coeffs_s(2, 3)*u_s(i, 4, b)*v_s(i, 4, b) &
               + dud_coeffs_s(3, 3)*u(i, 1, b)*v(i, 1, b) &
               + dud_coeffs_s(4, 3)*u(i, 2, b)*v(i, 2, b) &
               + dud_coeffs_s(5, 3)*u(i, 3, b)*v(i, 3, b) &
               + dud_coeffs_s(6, 3)*u(i, 4, b)*v(i, 4, b) &
               + dud_coeffs_s(7, 3)*u(i, 5, b)*v(i, 5, b) &
               + dud_coeffs_s(8, 3)*u(i, 6, b)*v(i, 6, b) &
               + dud_coeffs_s(9, 3)*u(i, 7, b)*v(i, 7, b)
    dud(i, 3, b) = dud_fw(3)*(temp_dud - dud_af(3)*dud(i, 2, b))
    temp_d2u = d2u_coeffs_s(1, 3)*u_s(i, 3, b) &
               + d2u_coeffs_s(2, 3)*u_s(i, 4, b) &
               + d2u_coeffs_s(3, 3)*u(i, 1, b) &
               + d2u_coeffs_s(4, 3)*u(i, 2, b) &
               + d2u_coeffs_s(5, 3)*u(i, 3, b) &
               + d2u_coeffs_s(6, 3)*u(i, 4, b) &
               + d2u_coeffs_s(7, 3)*u(i, 5, b) &
               + d2u_coeffs_s(8, 3)*u(i, 6, b) &
               + d2u_coeffs_s(9, 3)*u(i, 7, b)
    d2u(i, 3, b) = d2u_fw(3)*(temp_d2u - d2u_af(3)*d2u(i, 2, b))
    ! j = 4
    temp_du = du_coeffs_s(1, 4)*u_s(i, 4, b) &
              + du_coeffs_s(2, 4)*u(i, 1, b) &
              + du_coeffs_s(3, 4)*u(i, 2, b) &
              + du_coeffs_s(4, 4)*u(i, 3, b) &
              + du_coeffs_s(5, 4)*u(i, 4, b) &
              + du_coeffs_s(6, 4)*u(i, 5, b) &
              + du_coeffs_s(7, 4)*u(i, 6, b) &
              + du_coeffs_s(8, 4)*u(i, 7, b) &
              + du_coeffs_s(9, 4)*u(i, 8, b)
    du(i, 4, b) = du_fw(4)*(temp_du - du_af(3)*du(i, 3, b))
    temp_dud = dud_coeffs_s(1, 4)*u_s(i, 4, b)*v_s(i, 4, b) &
               + dud_coeffs_s(2, 4)*u(i, 1, b)*v(i, 1, b) &
               + dud_coeffs_s(3, 4)*u(i, 2, b)*v(i, 2, b) &
               + dud_coeffs_s(4, 4)*u(i, 3, b)*v(i, 3, b) &
               + dud_coeffs_s(5, 4)*u(i, 4, b)*v(i, 4, b) &
               + dud_coeffs_s(6, 4)*u(i, 5, b)*v(i, 5, b) &
               + dud_coeffs_s(7, 4)*u(i, 6, b)*v(i, 6, b) &
               + dud_coeffs_s(8, 4)*u(i, 7, b)*v(i, 7, b) &
               + dud_coeffs_s(9, 4)*u(i, 8, b)*v(i, 8, b)
    dud(i, 4, b) = dud_fw(4)*(temp_dud - dud_af(3)*dud(i, 3, b))
    temp_d2u = d2u_coeffs_s(1, 4)*u_s(i, 4, b) &
               + d2u_coeffs_s(2, 4)*u(i, 1, b) &
               + d2u_coeffs_s(3, 4)*u(i, 2, b) &
               + d2u_coeffs_s(4, 4)*u(i, 3, b) &
               + d2u_coeffs_s(5, 4)*u(i, 4, b) &
               + d2u_coeffs_s(6, 4)*u(i, 5, b) &
               + d2u_coeffs_s(7, 4)*u(i, 6, b) &
               + d2u_coeffs_s(8, 4)*u(i, 7, b) &
               + d2u_coeffs_s(9, 4)*u(i, 8, b)
    d2u(i, 4, b) = d2u_fw(4)*(temp_d2u - d2u_af(3)*d2u(i, 3, b))

    du_alpha = du_af(5)
    dud_alpha = dud_af(5)
    d2u_alpha = d2u_af(5)

    ! store bulk coeffs in the registers
    du_c_m4 = du_coeffs(1); du_c_m3 = du_coeffs(2)
    du_c_m2 = du_coeffs(3); du_c_m1 = du_coeffs(4)
    du_c_j = du_coeffs(5)
    du_c_p1 = du_coeffs(6); du_c_p2 = du_coeffs(7)
    du_c_p3 = du_coeffs(8); du_c_p4 = du_coeffs(9)

    dud_c_m4 = dud_coeffs(1); dud_c_m3 = dud_coeffs(2)
    dud_c_m2 = dud_coeffs(3); dud_c_m1 = dud_coeffs(4)
    dud_c_j = dud_coeffs(5)
    dud_c_p1 = dud_coeffs(6); dud_c_p2 = dud_coeffs(7)
    dud_c_p3 = dud_coeffs(8); dud_c_p4 = dud_coeffs(9)

    d2u_c_m4 = d2u_coeffs(1); d2u_c_m3 = d2u_coeffs(2)
    d2u_c_m2 = d2u_coeffs(3); d2u_c_m1 = d2u_coeffs(4)
    d2u_c_j = d2u_coeffs(5)
    d2u_c_p1 = d2u_coeffs(6); d2u_c_p2 = d2u_coeffs(7)
    d2u_c_p3 = d2u_coeffs(8); d2u_c_p4 = d2u_coeffs(9)

    ! It is better to access d?(i, j - 1, b) via old_d?
    old_du = du(i, 4, b)
    old_dud = dud(i, 4, b)
    old_d2u = d2u(i, 4, b)

    ! Populate registers with the u and v stencils
    u_m4 = u(i, 1, b); u_m3 = u(i, 2, b)
    u_m2 = u(i, 3, b); u_m1 = u(i, 4, b)
    u_j = u(i, 5, b); u_p1 = u(i, 6, b)
    u_p2 = u(i, 7, b); u_p3 = u(i, 8, b)
    v_m4 = v(i, 1, b); v_m3 = v(i, 2, b)
    v_m2 = v(i, 3, b); v_m1 = v(i, 4, b)
    v_j = v(i, 5, b); v_p1 = v(i, 6, b)
    v_p2 = v(i, 7, b); v_p3 = v(i, 8, b)

    do j = 5, n_rhs - 4
      u_p4 = u(i, j + 4, b); v_p4 = v(i, j + 4, b)

      ! du
      temp_du = du_c_m4*u_m4 + du_c_m3*u_m3 + du_c_m2*u_m2 + du_c_m1*u_m1 &
                + du_c_j*u_j &
                + du_c_p1*u_p1 + du_c_p2*u_p2 + du_c_p3*u_p3 + du_c_p4*u_p4
      du(i, j, b) = du_fw(j)*(temp_du - du_alpha*old_du)
      old_du = du(i, j, b)

      ! dud
      temp_dud = dud_c_m4*u_m4*v_m4 + dud_c_m3*u_m3*v_m3 &
                 + dud_c_m2*u_m2*v_m2 + dud_c_m1*u_m1*v_m1 &
                 + dud_c_j*u_j*v_j &
                 + dud_c_p1*u_p1*v_p1 + dud_c_p2*u_p2*v_p2 &
                 + dud_c_p3*u_p3*v_p3 + dud_c_p4*u_p4*v_p4
      dud(i, j, b) = dud_fw(j)*(temp_dud - dud_alpha*old_dud)
      old_dud = dud(i, j, b)

      ! d2u
      temp_d2u = d2u_c_m4*u_m4 + d2u_c_m3*u_m3 &
                 + d2u_c_m2*u_m2 + d2u_c_m1*u_m1 &
                 + d2u_c_j*u_j &
                 + d2u_c_p1*u_p1 + d2u_c_p2*u_p2 &
                 + d2u_c_p3*u_p3 + d2u_c_p4*u_p4
      d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_alpha*old_d2u)
      old_d2u = d2u(i, j, b)

      ! Prepare registers for the next step
      u_m4 = u_m3; u_m3 = u_m2; u_m2 = u_m1; u_m1 = u_j
      u_j = u_p1; u_p1 = u_p2; u_p2 = u_p3; u_p3 = u_p4
      v_m4 = v_m3; v_m3 = v_m2; v_m2 = v_m1; v_m1 = v_j
      v_j = v_p1; v_p1 = v_p2; v_p2 = v_p3; v_p3 = v_p4
    end do

    j = n_rhs - 3
    temp_du = du_coeffs_e(1, 1)*u(i, j - 4, b) &
              + du_coeffs_e(2, 1)*u(i, j - 3, b) &
              + du_coeffs_e(3, 1)*u(i, j - 2, b) &
              + du_coeffs_e(4, 1)*u(i, j - 1, b) &
              + du_coeffs_e(5, 1)*u(i, j, b) &
              + du_coeffs_e(6, 1)*u(i, j + 1, b) &
              + du_coeffs_e(7, 1)*u(i, j + 2, b) &
              + du_coeffs_e(8, 1)*u(i, j + 3, b) &
              + du_coeffs_e(9, 1)*u_e(i, 1, b)
    du(i, j, b) = du_fw(j)*(temp_du - du_af(j)*du(i, j - 1, b))
    temp_dud = dud_coeffs_e(1, 1)*u(i, j - 4, b)*v(i, j - 4, b) &
               + dud_coeffs_e(2, 1)*u(i, j - 3, b)*v(i, j - 3, b) &
               + dud_coeffs_e(3, 1)*u(i, j - 2, b)*v(i, j - 2, b) &
               + dud_coeffs_e(4, 1)*u(i, j - 1, b)*v(i, j - 1, b) &
               + dud_coeffs_e(5, 1)*u(i, j, b)*v(i, j, b) &
               + dud_coeffs_e(6, 1)*u(i, j + 1, b)*v(i, j + 1, b) &
               + dud_coeffs_e(7, 1)*u(i, j + 2, b)*v(i, j + 2, b) &
               + dud_coeffs_e(8, 1)*u(i, j + 3, b)*v(i, j + 3, b) &
               + dud_coeffs_e(9, 1)*u_e(i, 1, b)*v_e(i, 1, b)
    dud(i, j, b) = dud_fw(j)*(temp_dud - dud_af(j)*dud(i, j - 1, b))
    temp_d2u = d2u_coeffs_e(1, 1)*u(i, j - 4, b) &
               + d2u_coeffs_e(2, 1)*u(i, j - 3, b) &
               + d2u_coeffs_e(3, 1)*u(i, j - 2, b) &
               + d2u_coeffs_e(4, 1)*u(i, j - 1, b) &
               + d2u_coeffs_e(5, 1)*u(i, j, b) &
               + d2u_coeffs_e(6, 1)*u(i, j + 1, b) &
               + d2u_coeffs_e(7, 1)*u(i, j + 2, b) &
               + d2u_coeffs_e(8, 1)*u(i, j + 3, b) &
               + d2u_coeffs_e(9, 1)*u_e(i, 1, b)
    d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_af(j)*d2u(i, j - 1, b))
    j = n_rhs - 2
    temp_du = du_coeffs_e(1, 2)*u(i, j - 4, b) &
              + du_coeffs_e(2, 2)*u(i, j - 3, b) &
              + du_coeffs_e(3, 2)*u(i, j - 2, b) &
              + du_coeffs_e(4, 2)*u(i, j - 1, b) &
              + du_coeffs_e(5, 2)*u(i, j, b) &
              + du_coeffs_e(6, 2)*u(i, j + 1, b) &
              + du_coeffs_e(7, 2)*u(i, j + 2, b) &
              + du_coeffs_e(8, 2)*u_e(i, 1, b) &
              + du_coeffs_e(9, 2)*u_e(i, 2, b)
    du(i, j, b) = du_fw(j)*(temp_du - du_af(j)*du(i, j - 1, b))
    temp_dud = dud_coeffs_e(1, 2)*u(i, j - 4, b)*v(i, j - 4, b) &
               + dud_coeffs_e(2, 2)*u(i, j - 3, b)*v(i, j - 3, b) &
               + dud_coeffs_e(3, 2)*u(i, j - 2, b)*v(i, j - 2, b) &
               + dud_coeffs_e(4, 2)*u(i, j - 1, b)*v(i, j - 1, b) &
               + dud_coeffs_e(5, 2)*u(i, j, b)*v(i, j, b) &
               + dud_coeffs_e(6, 2)*u(i, j + 1, b)*v(i, j + 1, b) &
               + dud_coeffs_e(7, 2)*u(i, j + 2, b)*v(i, j + 2, b) &
               + dud_coeffs_e(8, 2)*u_e(i, 1, b)*v_e(i, 1, b) &
               + dud_coeffs_e(9, 2)*u_e(i, 2, b)*v_e(i, 2, b)
    dud(i, j, b) = dud_fw(j)*(temp_dud - dud_af(j)*dud(i, j - 1, b))
    temp_d2u = d2u_coeffs_e(1, 2)*u(i, j - 4, b) &
               + d2u_coeffs_e(2, 2)*u(i, j - 3, b) &
               + d2u_coeffs_e(3, 2)*u(i, j - 2, b) &
               + d2u_coeffs_e(4, 2)*u(i, j - 1, b) &
               + d2u_coeffs_e(5, 2)*u(i, j, b) &
               + d2u_coeffs_e(6, 2)*u(i, j + 1, b) &
               + d2u_coeffs_e(7, 2)*u(i, j + 2, b) &
               + d2u_coeffs_e(8, 2)*u_e(i, 1, b) &
               + d2u_coeffs_e(9, 2)*u_e(i, 2, b)
    d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_af(j)*d2u(i, j - 1, b))
    j = n_rhs - 1
    temp_du = du_coeffs_e(1, 3)*u(i, j - 4, b) &
              + du_coeffs_e(2, 3)*u(i, j - 3, b) &
              + du_coeffs_e(3, 3)*u(i, j - 2, b) &
              + du_coeffs_e(4, 3)*u(i, j - 1, b) &
              + du_coeffs_e(5, 3)*u(i, j, b) &
              + du_coeffs_e(6, 3)*u(i, j + 1, b) &
              + du_coeffs_e(7, 3)*u_e(i, 1, b) &
              + du_coeffs_e(8, 3)*u_e(i, 2, b) &
              + du_coeffs_e(9, 3)*u_e(i, 3, b)
    du(i, j, b) = du_fw(j)*(temp_du - du_af(j)*du(i, j - 1, b))
    temp_dud = dud_coeffs_e(1, 3)*u(i, j - 4, b)*v(i, j - 4, b) &
               + dud_coeffs_e(2, 3)*u(i, j - 3, b)*v(i, j - 3, b) &
               + dud_coeffs_e(3, 3)*u(i, j - 2, b)*v(i, j - 2, b) &
               + dud_coeffs_e(4, 3)*u(i, j - 1, b)*v(i, j - 1, b) &
               + dud_coeffs_e(5, 3)*u(i, j, b)*v(i, j, b) &
               + dud_coeffs_e(6, 3)*u(i, j + 1, b)*v(i, j + 1, b) &
               + dud_coeffs_e(7, 3)*u_e(i, 1, b)*v_e(i, 1, b) &
               + dud_coeffs_e(8, 3)*u_e(i, 2, b)*v_e(i, 2, b) &
               + dud_coeffs_e(9, 3)*u_e(i, 3, b)*v_e(i, 3, b)
    dud(i, j, b) = dud_fw(j)*(temp_dud - dud_af(j)*dud(i, j - 1, b))
    temp_d2u = d2u_coeffs_e(1, 3)*u(i, j - 4, b) &
               + d2u_coeffs_e(2, 3)*u(i, j - 3, b) &
               + d2u_coeffs_e(3, 3)*u(i, j - 2, b) &
               + d2u_coeffs_e(4, 3)*u(i, j - 1, b) &
               + d2u_coeffs_e(5, 3)*u(i, j, b) &
               + d2u_coeffs_e(6, 3)*u(i, j + 1, b) &
               + d2u_coeffs_e(7, 3)*u_e(i, 1, b) &
               + d2u_coeffs_e(8, 3)*u_e(i, 2, b) &
               + d2u_coeffs_e(9, 3)*u_e(i, 3, b)
    d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_af(j)*d2u(i, j - 1, b))
    j = n_rhs
    temp_du = du_coeffs_e(1, 4)*u(i, j - 4, b) &
              + du_coeffs_e(2, 4)*u(i, j - 3, b) &
              + du_coeffs_e(3, 4)*u(i, j - 2, b) &
              + du_coeffs_e(4, 4)*u(i, j - 1, b) &
              + du_coeffs_e(5, 4)*u(i, j, b) &
              + du_coeffs_e(6, 4)*u_e(i, 1, b) &
              + du_coeffs_e(7, 4)*u_e(i, 2, b) &
              + du_coeffs_e(8, 4)*u_e(i, 3, b) &
              + du_coeffs_e(9, 4)*u_e(i, 4, b)
    du(i, j, b) = du_fw(j)*(temp_du - du_af(j)*du(i, j - 1, b))
    temp_dud = dud_coeffs_e(1, 4)*u(i, j - 4, b)*v(i, j - 4, b) &
               + dud_coeffs_e(2, 4)*u(i, j - 3, b)*v(i, j - 3, b) &
               + dud_coeffs_e(3, 4)*u(i, j - 2, b)*v(i, j - 2, b) &
               + dud_coeffs_e(4, 4)*u(i, j - 1, b)*v(i, j - 1, b) &
               + dud_coeffs_e(5, 4)*u(i, j, b)*v(i, j, b) &
               + dud_coeffs_e(6, 4)*u_e(i, 1, b)*v_e(i, 1, b) &
               + dud_coeffs_e(7, 4)*u_e(i, 2, b)*v_e(i, 2, b) &
               + dud_coeffs_e(8, 4)*u_e(i, 3, b)*v_e(i, 3, b) &
               + dud_coeffs_e(9, 4)*u_e(i, 4, b)*v_e(i, 4, b)
    dud(i, j, b) = dud_fw(j)*(temp_dud - dud_af(j)*dud(i, j - 1, b))
    temp_d2u = d2u_coeffs_e(1, 4)*u(i, j - 4, b) &
               + d2u_coeffs_e(2, 4)*u(i, j - 3, b) &
               + d2u_coeffs_e(3, 4)*u(i, j - 2, b) &
               + d2u_coeffs_e(4, 4)*u(i, j - 1, b) &
               + d2u_coeffs_e(5, 4)*u(i, j, b) &
               + d2u_coeffs_e(6, 4)*u_e(i, 1, b) &
               + d2u_coeffs_e(7, 4)*u_e(i, 2, b) &
               + d2u_coeffs_e(8, 4)*u_e(i, 3, b) &
               + d2u_coeffs_e(9, 4)*u_e(i, 4, b)
    d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_af(j)*d2u(i, j - 1, b))

    send_du_e(i, 1, b) = du(i, n_tds, b)
    send_dud_e(i, 1, b) = dud(i, n_tds, b)
    send_d2u_e(i, 1, b) = d2u(i, n_tds, b)

    ! Backward pass of the hybrid algorithm
    do j = n_tds - 2, 2, -1
      du(i, j, b) = du(i, j, b) - du_bw(j)*du(i, j + 1, b)
      dud(i, j, b) = dud(i, j, b) - dud_bw(j)*dud(i, j + 1, b)
      d2u(i, j, b) = d2u(i, j, b) - d2u_bw(j)*d2u(i, j + 1, b)
    end do
    du(i, 1, b) = du_last_r*(du(i, 1, b) - du_bw(1)*du(i, 2, b))
    dud(i, 1, b) = dud_last_r*(dud(i, 1, b) - dud_bw(1)*dud(i, 2, b))
    d2u(i, 1, b) = d2u_last_r*(d2u(i, 1, b) - d2u_bw(1)*d2u(i, 2, b))

    send_du_s(i, 1, b) = du(i, 1, b)
    send_dud_s(i, 1, b) = dud(i, 1, b)
    send_d2u_s(i, 1, b) = d2u(i, 1, b)

  end subroutine transeq_3fused_dist

  attributes(global) subroutine transeq_3fused_subs( &
    r_du, conv, dud, d2u, &
    recv_du_s, recv_du_e, recv_dud_s, recv_dud_e, recv_d2u_s, recv_d2u_e, &
    n, nu, du_sa, du_sc, du_strch, dud_sa, dud_sc, dud_strch, &
    d2u_sa, d2u_sc, d2u_strch, d2u_strch_cor &
    )
    implicit none

    ! Arguments
    !> The result array, it stores 'du' first then its overwritten
    real(dp), device, intent(inout), dimension(:, :, :) :: r_du
    real(dp), device, intent(in), dimension(:, :, :) :: conv, dud, d2u
    real(dp), device, intent(in), dimension(:, :, :) :: &
      recv_du_s, recv_du_e, recv_dud_s, recv_dud_e, recv_d2u_s, recv_d2u_e
    integer, value, intent(in) :: n
    real(dp), value, intent(in) :: nu
    real(dp), device, intent(in), dimension(:) :: du_sa, du_sc, du_strch, &
                                                  dud_sa, dud_sc, dud_strch, &
                                                  d2u_sa, d2u_sc, d2u_strch, &
                                                  d2u_strch_cor

    ! Local variables
    integer :: i, j, b
    real(dp) :: ur, bl, recp
    real(dp) :: du_temp, dud_temp, d2u_temp
    real(dp) :: du_s, du_e, dud_s, dud_e, d2u_s, d2u_e

    i = threadIdx%x
    b = blockIdx%x

    ! A small trick we do here is valid for symmetric Toeplitz matrices.
    ! In our case our matrices satisfy this criteria in the (5:n-4) region
    ! and as long as a rank has around at least 20 entries the assumptions
    ! we make here are perfectly valid.

    ! bl is the bottom left entry in the 2x2 matrix
    ! ur is the upper right entry in the 2x2 matrix

    ! Start
    ! At the start we have the 'bl', and assume 'ur'
    ! first derivative
    bl = du_sa(1)
    ur = du_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    du_s = recp*(r_du(i, 1, b) - bl*recv_du_s(i, 1, b))

    ! first derivative (u*v)
    bl = dud_sa(1)
    ur = dud_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    dud_s = recp*(dud(i, 1, b) - bl*recv_dud_s(i, 1, b))

    ! second derivative
    bl = d2u_sa(1)
    ur = d2u_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    d2u_s = recp*(d2u(i, 1, b) - bl*recv_d2u_s(i, 1, b))

    ! End
    ! At the end we have the 'ur', and assume 'bl'
    ! first derivative
    bl = du_sc(n)
    ur = du_sc(n)
    recp = 1._dp/(1._dp - ur*bl)
    du_e = recp*(r_du(i, n, b) - ur*recv_du_e(i, 1, b))

    bl = dud_sc(n)
    ur = dud_sc(n)
    recp = 1._dp/(1._dp - ur*bl)
    dud_e = recp*(dud(i, n, b) - ur*recv_dud_e(i, 1, b))

    ! second derivative
    bl = d2u_sc(n)
    ur = d2u_sc(n)
    recp = 1._dp/(1._dp - ur*bl)
    d2u_e = recp*(d2u(i, n, b) - ur*recv_d2u_e(i, 1, b))

    ! final substitution
    r_du(i, 1, b) = -0.5_dp*(conv(i, 1, b)*du_s*du_strch(1) &
                             + dud_s*dud_strch(1)) &
                    + nu*(d2u_s*d2u_strch(1) &
                          + du_s*du_strch(1)*d2u_strch_cor(1))
    do j = 2, n - 1
      du_temp = (r_du(i, j, b) - du_sa(j)*du_s - du_sc(j)*du_e)*du_strch(j)
      dud_temp = (dud(i, j, b) - dud_sa(j)*dud_s - dud_sc(j)*dud_e) &
                 *dud_strch(j)
      d2u_temp = (d2u(i, j, b) - d2u_sa(j)*d2u_s - d2u_sc(j)*d2u_e) &
                 *d2u_strch(j) + du_temp*d2u_strch_cor(j)
      r_du(i, j, b) = -0.5_dp*(conv(i, j, b)*du_temp + dud_temp) + nu*d2u_temp
    end do
    r_du(i, n, b) = -0.5_dp*(conv(i, n, b)*du_e*du_strch(n) &
                             + dud_e*dud_strch(n)) &
                    + nu*(d2u_e*d2u_strch(n) &
                          + du_s*du_strch(n)*d2u_strch_cor(n))

  end subroutine transeq_3fused_subs

end module m_cuda_kernels_dist
