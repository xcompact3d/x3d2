module m_omp_kernels_thom
  use m_common, only: dp
  use m_omp_common, only: SZ

  implicit none

contains

  subroutine der_univ_thom(du, u, n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, &
                           thom_f, thom_s, thom_w, strch)
    implicit none

    real(dp), dimension(:, :), intent(out) :: du
    real(dp), dimension(:, :), intent(in) :: u
    integer, intent(in) :: n_tds, n_rhs
    real(dp), intent(in), dimension(:, :) :: coeffs_s, coeffs_e ! start/end
    real(dp), intent(in), dimension(:) :: coeffs
    real(dp), intent(in), dimension(:) :: thom_f, thom_s, thom_w, strch

    integer :: i, j
    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    ! Forward pass
    !$omp simd
    do i = 1, SZ
      du(i, 1) = coeffs_s(5, 1)*u(i, 1) &
                 + coeffs_s(6, 1)*u(i, 2) &
                 + coeffs_s(7, 1)*u(i, 3) &
                 + coeffs_s(8, 1)*u(i, 4) &
                 + coeffs_s(9, 1)*u(i, 5)
      du(i, 2) = coeffs_s(4, 2)*u(i, 1) &
                 + coeffs_s(5, 2)*u(i, 2) &
                 + coeffs_s(6, 2)*u(i, 3) &
                 + coeffs_s(7, 2)*u(i, 4) &
                 + coeffs_s(8, 2)*u(i, 5) &
                 + coeffs_s(9, 2)*u(i, 6) &
                 - du(i, 1)*thom_s(2)
      du(i, 3) = coeffs_s(3, 3)*u(i, 1) &
                 + coeffs_s(4, 3)*u(i, 2) &
                 + coeffs_s(5, 3)*u(i, 3) &
                 + coeffs_s(6, 3)*u(i, 4) &
                 + coeffs_s(7, 3)*u(i, 5) &
                 + coeffs_s(8, 3)*u(i, 6) &
                 + coeffs_s(9, 3)*u(i, 7) &
                 - du(i, 2)*thom_s(3)
      du(i, 4) = coeffs_s(2, 4)*u(i, 1) &
                 + coeffs_s(3, 4)*u(i, 2) &
                 + coeffs_s(4, 4)*u(i, 3) &
                 + coeffs_s(5, 4)*u(i, 4) &
                 + coeffs_s(6, 4)*u(i, 5) &
                 + coeffs_s(7, 4)*u(i, 6) &
                 + coeffs_s(8, 4)*u(i, 7) &
                 + coeffs_s(9, 4)*u(i, 8) &
                 - du(i, 3)*thom_s(4)
    end do
    !$omp end simd

    do j = 5, n_rhs - 4
      !$omp simd
      do i = 1, SZ
        du(i, j) = c_m4*u(i, j - 4) + c_m3*u(i, j - 3) &
                   + c_m2*u(i, j - 2) + c_m1*u(i, j - 1) &
                   + c_j*u(i, j) &
                   + c_p1*u(i, j + 1) + c_p2*u(i, j + 2) &
                   + c_p3*u(i, j + 3) + c_p4*u(i, j + 4) &
                   - du(i, j - 1)*thom_s(j)
      end do
      !$omp end simd
    end do

    !$omp simd
    do i = 1, SZ
      j = n_rhs - 3
      du(i, j) = coeffs_e(1, 1)*u(i, j - 4) &
                 + coeffs_e(2, 1)*u(i, j - 3) &
                 + coeffs_e(3, 1)*u(i, j - 2) &
                 + coeffs_e(4, 1)*u(i, j - 1) &
                 + coeffs_e(5, 1)*u(i, j) &
                 + coeffs_e(6, 1)*u(i, j + 1) &
                 + coeffs_e(7, 1)*u(i, j + 2) &
                 + coeffs_e(8, 1)*u(i, j + 3) &
                 - du(i, j - 1)*thom_s(j)
      j = n_rhs - 2
      du(i, j) = coeffs_e(1, 2)*u(i, j - 4) &
                 + coeffs_e(2, 2)*u(i, j - 3) &
                 + coeffs_e(3, 2)*u(i, j - 2) &
                 + coeffs_e(4, 2)*u(i, j - 1) &
                 + coeffs_e(5, 2)*u(i, j) &
                 + coeffs_e(6, 2)*u(i, j + 1) &
                 + coeffs_e(7, 2)*u(i, j + 2) &
                 - du(i, j - 1)*thom_s(j)
      j = n_rhs - 1
      du(i, j) = coeffs_e(1, 3)*u(i, j - 4) &
                 + coeffs_e(2, 3)*u(i, j - 3) &
                 + coeffs_e(3, 3)*u(i, j - 2) &
                 + coeffs_e(4, 3)*u(i, j - 1) &
                 + coeffs_e(5, 3)*u(i, j) &
                 + coeffs_e(6, 3)*u(i, j + 1) &
                 - du(i, j - 1)*thom_s(j)
      j = n_rhs
      du(i, j) = coeffs_e(1, 4)*u(i, j - 4) &
                 + coeffs_e(2, 4)*u(i, j - 3) &
                 + coeffs_e(3, 4)*u(i, j - 2) &
                 + coeffs_e(4, 4)*u(i, j - 1) &
                 + coeffs_e(5, 4)*u(i, j) &
                 - du(i, j - 1)*thom_s(j)
    end do
    !$omp end simd

    ! Backward pass
    !$omp simd
    do i = 1, SZ
      du(i, n_tds) = du(i, n_tds)*thom_w(n_tds)*strch(n_tds)
    end do
    !$omp end simd
    do j = n_tds - 1, 1, -1
      !$omp simd
      do i = 1, SZ
        ! du(j) = (du(j) - f*du(j+1)/strch(j))*w*strch(j)
        du(i, j) = (du(i, j)*strch(j) - thom_f(j)*du(i, j + 1))*thom_w(j)
      end do
      !$omp end simd
    end do

  end subroutine der_univ_thom

  subroutine der_univ_thom_per( &
    du, u, n, coeffs, alpha, thom_f, thom_s, thom_w, thom_p, strch &
    )
    implicit none

    real(dp), dimension(:, :), intent(out) :: du
    real(dp), dimension(:, :), intent(in) :: u
    integer, intent(in) :: n
    real(dp), intent(in), dimension(:) :: coeffs
    real(dp), intent(in) :: alpha
    real(dp), intent(in), dimension(:) :: thom_f, thom_s, thom_w, thom_p, strch

    integer :: i, j
    integer :: jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4
    real(dp), dimension(SZ) :: ss

    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    ! Forward pass
    do j = 1, n
      jm4 = modulo(j - 5, n) + 1
      jm3 = modulo(j - 4, n) + 1
      jm2 = modulo(j - 3, n) + 1
      jm1 = modulo(j - 2, n) + 1
      jp1 = modulo(j - n, n) + 1
      jp2 = modulo(j - n + 1, n) + 1
      jp3 = modulo(j - n + 2, n) + 1
      jp4 = modulo(j - n + 3, n) + 1

      !$omp simd
      do i = 1, SZ
        du(i, j) = c_m4*u(i, jm4) + c_m3*u(i, jm3) &
                   + c_m2*u(i, jm2) + c_m1*u(i, jm1) &
                   + c_j*u(i, j) &
                   + c_p1*u(i, jp1) + c_p2*u(i, jp2) &
                   + c_p3*u(i, jp3) + c_p4*u(i, jp4) &
                   - du(i, jm1)*thom_s(j)
      end do
      !$omp end simd
    end do

    ! Backward pass
    !$omp simd
    do i = 1, SZ
      du(i, n) = du(i, n)*thom_w(n)
    end do
    !$omp end simd
    do j = n - 1, 1, -1
      !$omp simd
      do i = 1, SZ
        du(i, j) = (du(i, j) - thom_f(j)*du(i, j + 1))*thom_w(j)
      end do
      !$omp end simd
    end do

    ! Periodic final pass
    !$omp simd
    do i = 1, SZ
      ss(i) = (du(i, 1) - alpha*du(i, n)) &
              /(1.0_dp + thom_p(1) - alpha*thom_p(n))
    end do
    !$omp end simd
    do j = 1, n
      !$omp simd
      do i = 1, SZ
        du(i, j) = (du(i, j) - ss(i)*thom_p(j))*strch(j)
      end do
      !$omp end simd
    end do

  end subroutine der_univ_thom_per

end module m_omp_kernels_thom
