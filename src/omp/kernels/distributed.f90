module m_omp_kernels_dist
  use omp_lib

  use m_common, only: dp
  use m_omp_common, only: SZ

  implicit none

contains

  subroutine der_univ_dist( &
    du, send_u_s, send_u_e, u, u_s, u_e, coeffs_s, coeffs_e, coeffs, n, &
    ffr, fbc, faf &
    )
    implicit none

    ! Arguments
    real(dp), intent(out), dimension(:, :) :: du, send_u_s, send_u_e
    real(dp), intent(in), dimension(:, :) :: u, u_s, u_e
    real(dp), intent(in), dimension(:, :) :: coeffs_s, coeffs_e ! start/end
    real(dp), intent(in), dimension(:) :: coeffs
    integer, intent(in) :: n
    real(dp), intent(in), dimension(:) :: ffr, fbc, faf

    ! Local variables
    integer :: i, j!, b

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, &
                alpha, last_r

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)
    last_r = ffr(1)

    !$omp simd
    do i = 1, SZ
      du(i, 1) = coeffs_s(1, 1)*u_s(i, 1) &
                 + coeffs_s(2, 1)*u_s(i, 2) &
                 + coeffs_s(3, 1)*u_s(i, 3) &
                 + coeffs_s(4, 1)*u_s(i, 4) &
                 + coeffs_s(5, 1)*u(i, 1) &
                 + coeffs_s(6, 1)*u(i, 2) &
                 + coeffs_s(7, 1)*u(i, 3) &
                 + coeffs_s(8, 1)*u(i, 4) &
                 + coeffs_s(9, 1)*u(i, 5)
      du(i, 1) = du(i, 1)*faf(1)
      du(i, 2) = coeffs_s(1, 2)*u_s(i, 2) &
                 + coeffs_s(2, 2)*u_s(i, 3) &
                 + coeffs_s(3, 2)*u_s(i, 4) &
                 + coeffs_s(4, 2)*u(i, 1) &
                 + coeffs_s(5, 2)*u(i, 2) &
                 + coeffs_s(6, 2)*u(i, 3) &
                 + coeffs_s(7, 2)*u(i, 4) &
                 + coeffs_s(8, 2)*u(i, 5) &
                 + coeffs_s(9, 2)*u(i, 6)
      du(i, 2) = du(i, 2)*faf(2)
      du(i, 3) = coeffs_s(1, 3)*u_s(i, 3) &
                 + coeffs_s(2, 3)*u_s(i, 4) &
                 + coeffs_s(3, 3)*u(i, 1) &
                 + coeffs_s(4, 3)*u(i, 2) &
                 + coeffs_s(5, 3)*u(i, 3) &
                 + coeffs_s(6, 3)*u(i, 4) &
                 + coeffs_s(7, 3)*u(i, 5) &
                 + coeffs_s(8, 3)*u(i, 6) &
                 + coeffs_s(9, 3)*u(i, 7)
      du(i, 3) = ffr(3)*(du(i, 3) - faf(3)*du(i, 2))
      du(i, 4) = coeffs_s(1, 4)*u_s(i, 4) &
                 + coeffs_s(2, 4)*u(i, 1) &
                 + coeffs_s(3, 4)*u(i, 2) &
                 + coeffs_s(4, 4)*u(i, 3) &
                 + coeffs_s(5, 4)*u(i, 4) &
                 + coeffs_s(6, 4)*u(i, 5) &
                 + coeffs_s(7, 4)*u(i, 6) &
                 + coeffs_s(8, 4)*u(i, 7) &
                 + coeffs_s(9, 4)*u(i, 8)
      du(i, 4) = ffr(4)*(du(i, 4) - faf(4)*du(i, 3))
    end do
    !$omp end simd

    ! alpha is always the same in the bulk region for us
    alpha = faf(5)

    do j = 5, n - 4
      !$omp simd
      do i = 1, SZ
        du(i, j) = c_m4*u(i, j - 4) + c_m3*u(i, j - 3) &
                   + c_m2*u(i, j - 2) + c_m1*u(i, j - 1) &
                   + c_j*u(i, j) &
                   + c_p1*u(i, j + 1) + c_p2*u(i, j + 2) &
                   + c_p3*u(i, j + 3) + c_p4*u(i, j + 4)
        du(i, j) = ffr(j)*(du(i, j) - alpha*du(i, j - 1))
      end do
      !$omp end simd
    end do

    !$omp simd
    do i = 1, SZ
      j = n - 3
      du(i, j) = coeffs_e(1, 1)*u(i, j - 4) &
                 + coeffs_e(2, 1)*u(i, j - 3) &
                 + coeffs_e(3, 1)*u(i, j - 2) &
                 + coeffs_e(4, 1)*u(i, j - 1) &
                 + coeffs_e(5, 1)*u(i, j) &
                 + coeffs_e(6, 1)*u(i, j + 1) &
                 + coeffs_e(7, 1)*u(i, j + 2) &
                 + coeffs_e(8, 1)*u(i, j + 3) &
                 + coeffs_e(9, 1)*u_e(i, 1)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n - 2
      du(i, j) = coeffs_e(1, 2)*u(i, j - 4) &
                 + coeffs_e(2, 2)*u(i, j - 3) &
                 + coeffs_e(3, 2)*u(i, j - 2) &
                 + coeffs_e(4, 2)*u(i, j - 1) &
                 + coeffs_e(5, 2)*u(i, j) &
                 + coeffs_e(6, 2)*u(i, j + 1) &
                 + coeffs_e(7, 2)*u(i, j + 2) &
                 + coeffs_e(8, 2)*u_e(i, 1) &
                 + coeffs_e(9, 2)*u_e(i, 2)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n - 1
      du(i, j) = coeffs_e(1, 3)*u(i, j - 4) &
                 + coeffs_e(2, 3)*u(i, j - 3) &
                 + coeffs_e(3, 3)*u(i, j - 2) &
                 + coeffs_e(4, 3)*u(i, j - 1) &
                 + coeffs_e(5, 3)*u(i, j) &
                 + coeffs_e(6, 3)*u(i, j + 1) &
                 + coeffs_e(7, 3)*u_e(i, 1) &
                 + coeffs_e(8, 3)*u_e(i, 2) &
                 + coeffs_e(9, 3)*u_e(i, 3)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n
      du(i, j) = coeffs_e(1, 4)*u(i, j - 4) &
                 + coeffs_e(2, 4)*u(i, j - 3) &
                 + coeffs_e(3, 4)*u(i, j - 2) &
                 + coeffs_e(4, 4)*u(i, j - 1) &
                 + coeffs_e(5, 4)*u(i, j) &
                 + coeffs_e(6, 4)*u_e(i, 1) &
                 + coeffs_e(7, 4)*u_e(i, 2) &
                 + coeffs_e(8, 4)*u_e(i, 3) &
                 + coeffs_e(9, 4)*u_e(i, 4)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
    end do
    !$omp end simd

    !$omp simd
    do i = 1, SZ
      send_u_e(i, 1) = du(i, n)
    end do
    !$omp end simd

    ! Backward pass of the hybrid algorithm
    do j = n - 2, 2, -1
      !$omp simd
      do i = 1, SZ
        du(i, j) = du(i, j) - fbc(j)*du(i, j + 1)
      end do
      !$omp end simd
    end do
    !$omp simd
    do i = 1, SZ
      du(i, 1) = last_r*(du(i, 1) - fbc(1)*du(i, 2))
      send_u_s(i, 1) = du(i, 1)
    end do
    !$omp end simd

  end subroutine der_univ_dist

  subroutine der_univ_subs(du, recv_u_s, recv_u_e, n, dist_sa, dist_sc)
    implicit none

    ! Arguments
    real(dp), intent(out), dimension(:, :) :: du
    real(dp), intent(in), dimension(:, :) :: recv_u_s, recv_u_e
    real(dp), intent(in), dimension(:) :: dist_sa, dist_sc
    integer, intent(in) :: n

    ! Local variables
    integer :: i, j!, b
    real(dp) :: ur, bl, recp
    real(dp), dimension(SZ) :: du_s, du_e

    !$omp simd
    do i = 1, SZ
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
      du_s(i) = recp*(du(i, 1) - bl*recv_u_s(i, 1))

      ! End
      ! At the end we have the 'ur', and assume 'bl'
      bl = dist_sc(n)
      ur = dist_sc(n)
      recp = 1._dp/(1._dp - ur*bl)
      du_e(i) = recp*(du(i, n) - ur*recv_u_e(i, 1))
    end do
    !$omp end simd

    !$omp simd
    do i = 1, SZ
      du(i, 1) = du_s(i)
    end do
    !$omp end simd
    do j = 2, n - 1
      !$omp simd
      do i = 1, SZ
        du(i, j) = (du(i, j) - dist_sa(j)*du_s(i) - dist_sc(j)*du_e(i))
      end do
      !$omp end simd
    end do
    !$omp simd
    do i = 1, SZ
      du(i, n) = du_e(i)
    end do
    !$omp end simd

  end subroutine der_univ_subs

end module m_omp_kernels_dist
