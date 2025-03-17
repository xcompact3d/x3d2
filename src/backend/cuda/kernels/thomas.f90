module m_cuda_kernels_thom
  use cudafor

  use m_common, only: dp

  implicit none

contains

  attributes(global) subroutine der_univ_thom( &
    du, u, n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, &
    thom_f, thom_s, thom_w, strch &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    integer, value, intent(in) :: n_tds, n_rhs
    real(dp), device, intent(in), dimension(:, :) :: coeffs_s, coeffs_e
    real(dp), device, intent(in), dimension(:) :: coeffs
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w, strch

    integer :: i, j, b

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, temp_du

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    du(i, 1, b) = coeffs_s(5, 1)*u(i, 1, b) &
                  + coeffs_s(6, 1)*u(i, 2, b) &
                  + coeffs_s(7, 1)*u(i, 3, b) &
                  + coeffs_s(8, 1)*u(i, 4, b) &
                  + coeffs_s(9, 1)*u(i, 5, b)
    du(i, 1, b) = du(i, 1, b)
    du(i, 2, b) = coeffs_s(4, 2)*u(i, 1, b) &
                  + coeffs_s(5, 2)*u(i, 2, b) &
                  + coeffs_s(6, 2)*u(i, 3, b) &
                  + coeffs_s(7, 2)*u(i, 4, b) &
                  + coeffs_s(8, 2)*u(i, 5, b) &
                  + coeffs_s(9, 2)*u(i, 6, b)
    du(i, 2, b) = du(i, 2, b) - du(i, 1, b)*thom_s(2)
    du(i, 3, b) = coeffs_s(3, 3)*u(i, 1, b) &
                  + coeffs_s(4, 3)*u(i, 2, b) &
                  + coeffs_s(5, 3)*u(i, 3, b) &
                  + coeffs_s(6, 3)*u(i, 4, b) &
                  + coeffs_s(7, 3)*u(i, 5, b) &
                  + coeffs_s(8, 3)*u(i, 6, b) &
                  + coeffs_s(9, 3)*u(i, 7, b)
    du(i, 3, b) = du(i, 3, b) - du(i, 2, b)*thom_s(3)
    du(i, 4, b) = coeffs_s(2, 4)*u(i, 1, b) &
                  + coeffs_s(3, 4)*u(i, 2, b) &
                  + coeffs_s(4, 4)*u(i, 3, b) &
                  + coeffs_s(5, 4)*u(i, 4, b) &
                  + coeffs_s(6, 4)*u(i, 5, b) &
                  + coeffs_s(7, 4)*u(i, 6, b) &
                  + coeffs_s(8, 4)*u(i, 7, b) &
                  + coeffs_s(9, 4)*u(i, 8, b)
    du(i, 4, b) = du(i, 4, b) - du(i, 3, b)*thom_s(4)

    do j = 5, n_rhs - 4
      temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
      du(i, j, b) = temp_du - du(i, j - 1, b)*thom_s(j)
    end do

    j = n_rhs - 3
    du(i, j, b) = coeffs_e(1, 1)*u(i, j - 4, b) &
                  + coeffs_e(2, 1)*u(i, j - 3, b) &
                  + coeffs_e(3, 1)*u(i, j - 2, b) &
                  + coeffs_e(4, 1)*u(i, j - 1, b) &
                  + coeffs_e(5, 1)*u(i, j, b) &
                  + coeffs_e(6, 1)*u(i, j + 1, b) &
                  + coeffs_e(7, 1)*u(i, j + 2, b) &
                  + coeffs_e(8, 1)*u(i, j + 3, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n_rhs - 2
    du(i, j, b) = coeffs_e(1, 2)*u(i, j - 4, b) &
                  + coeffs_e(2, 2)*u(i, j - 3, b) &
                  + coeffs_e(3, 2)*u(i, j - 2, b) &
                  + coeffs_e(4, 2)*u(i, j - 1, b) &
                  + coeffs_e(5, 2)*u(i, j, b) &
                  + coeffs_e(6, 2)*u(i, j + 1, b) &
                  + coeffs_e(7, 2)*u(i, j + 2, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n_rhs - 1
    du(i, j, b) = coeffs_e(1, 3)*u(i, j - 4, b) &
                  + coeffs_e(2, 3)*u(i, j - 3, b) &
                  + coeffs_e(3, 3)*u(i, j - 2, b) &
                  + coeffs_e(4, 3)*u(i, j - 1, b) &
                  + coeffs_e(5, 3)*u(i, j, b) &
                  + coeffs_e(6, 3)*u(i, j + 1, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n_rhs
    du(i, j, b) = coeffs_e(1, 4)*u(i, j - 4, b) &
                  + coeffs_e(2, 4)*u(i, j - 3, b) &
                  + coeffs_e(3, 4)*u(i, j - 2, b) &
                  + coeffs_e(4, 4)*u(i, j - 1, b) &
                  + coeffs_e(5, 4)*u(i, j, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)

    ! Backward pass of the Thomas algorithm
    du(i, n_tds, b) = du(i, n_tds, b)*thom_w(n_tds)*strch(n_tds)
    do j = n_tds - 1, 1, -1
      ! du(j) = (du(j) - f*du(j+1)/strch(j))*w*strch(j)
      du(i, j, b) = (du(i, j, b)*strch(j) - thom_f(j)*du(i, j + 1, b)) &
                    *thom_w(j)
    end do

  end subroutine der_univ_thom

  attributes(global) subroutine der_univ_thom_per( &
    du, u, n, coeffs, alpha, thom_f, thom_s, thom_w, thom_p, strch &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    integer, value, intent(in) :: n
    real(dp), device, intent(in), dimension(:) :: coeffs
    real(dp), value, intent(in) :: alpha
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w, &
                                                  thom_p, strch

    integer :: i, j, b
    integer :: jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4
    real(dp) :: temp_du, ss

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    do j = 1, n
      jm4 = modulo(j - 5, n) + 1
      jm3 = modulo(j - 4, n) + 1
      jm2 = modulo(j - 3, n) + 1
      jm1 = modulo(j - 2, n) + 1
      jp1 = modulo(j - n, n) + 1
      jp2 = modulo(j - n + 1, n) + 1
      jp3 = modulo(j - n + 2, n) + 1
      jp4 = modulo(j - n + 3, n) + 1

      temp_du = c_m4*u(i, jm4, b) + c_m3*u(i, jm3, b) &
                + c_m2*u(i, jm2, b) + c_m1*u(i, jm1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, jp1, b) + c_p2*u(i, jp2, b) &
                + c_p3*u(i, jp3, b) + c_p4*u(i, jp4, b)
      du(i, j, b) = temp_du - du(i, jm1, b)*thom_s(j)
    end do

    ! Backward pass of the Thomas algorithm
    du(i, n, b) = du(i, n, b)*thom_w(n)
    do j = n - 1, 1, -1
      du(i, j, b) = (du(i, j, b) - thom_f(j)*du(i, j + 1, b))*thom_w(j)
    end do

    ! Periodic final pass
    ss = (du(i, 1, b) - alpha*du(i, n, b)) &
         /(1._dp + thom_p(1) - alpha*thom_p(n))
    do j = 1, n
      du(i, j, b) = (du(i, j, b) - ss*thom_p(j))*strch(j)
    end do

  end subroutine der_univ_thom_per

end module m_cuda_kernels_thom
