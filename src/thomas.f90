module m_omptgt_kernels_thom
  !! Thomas algorithm kernels for the OpenMP target backend, following
  !! cuda/kernels/thomas.f90. See m_omptgt_kernels_dist for the offload
  !! pattern: one `declare target` line solve per (i, b).
  use m_common, only: dp
  use m_omptgt_common, only: SZ

  implicit none

  private
  public :: der_univ_thom, der_univ_thom_per

contains

  subroutine der_univ_thom(du, u, n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, &
                           thom_f, thom_s, thom_w, strch, n, n_groups)
    !! Non-periodic Thomas solve along each line of a single subdomain.
    integer, intent(in) :: n_tds, n_rhs, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(in) :: u(SZ, n, n_groups)
    real(dp), intent(in) :: coeffs_s(9, 4), coeffs_e(9, 4), coeffs(9)
    real(dp), intent(in) :: thom_f(n_rhs), thom_s(n_rhs), thom_w(n_rhs), &
                            strch(n_tds)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, u, coeffs_s, coeffs_e, coeffs, thom_f, thom_s, &
    !$omp                 thom_w, strch)
    do b = 1, n_groups
      do i = 1, SZ
        call der_univ_thom_line(i, b, n, du, u, n_tds, n_rhs, coeffs_s, &
                                coeffs_e, coeffs, thom_f, thom_s, thom_w, &
                                strch)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_univ_thom

  subroutine der_univ_thom_line(i, b, n, du, u, n_tds, n_rhs, coeffs_s, &
                                coeffs_e, coeffs, thom_f, thom_s, thom_w, &
                                strch)
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds, n_rhs
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(in) :: u(SZ, n, *)
    real(dp), intent(in) :: coeffs_s(9, *), coeffs_e(9, *), coeffs(*)
    real(dp), intent(in) :: thom_f(*), thom_s(*), thom_w(*), strch(*)

    integer :: j, k, row, m
    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, temp_du

    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    ! Start rows: there is no halo, so row j only takes the stencil entries
    ! that land inside the line, k = 6 - j..9.
    do j = 1, 4
      temp_du = 0._dp
      do k = 6 - j, 9
        temp_du = temp_du + coeffs_s(k, j)*u(i, j - 5 + k, b)
      end do
      if (j == 1) then
        du(i, j, b) = temp_du
      else
        du(i, j, b) = temp_du - du(i, j - 1, b)*thom_s(j)
      end if
    end do

    do j = 5, n_rhs - 4
      temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
      du(i, j, b) = temp_du - du(i, j - 1, b)*thom_s(j)
    end do

    ! End rows take the entries up to n_rhs, k = 1..9 - row.
    do row = 1, 4
      j = n_rhs - 4 + row
      temp_du = 0._dp
      do k = 1, 9 - row
        m = j - 5 + k
        temp_du = temp_du + coeffs_e(k, row)*u(i, m, b)
      end do
      du(i, j, b) = temp_du - du(i, j - 1, b)*thom_s(j)
    end do

    ! Backward pass of the Thomas algorithm
    du(i, n_tds, b) = du(i, n_tds, b)*thom_w(n_tds)*strch(n_tds)
    do j = n_tds - 1, 1, -1
      du(i, j, b) = (du(i, j, b)*strch(j) - thom_f(j)*du(i, j + 1, b)) &
                    *thom_w(j)
    end do

  end subroutine der_univ_thom_line

  subroutine der_univ_thom_per(du, u, n_tds, coeffs, alpha, thom_f, thom_s, &
                               thom_w, thom_p, strch, n, n_groups)
    !! Periodic Thomas solve along each line of a single subdomain.
    integer, intent(in) :: n_tds, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(in) :: u(SZ, n, n_groups)
    real(dp), intent(in) :: coeffs(9)
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: thom_f(n_tds), thom_s(n_tds), thom_w(n_tds), &
                            thom_p(n_tds), strch(n_tds)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, u, coeffs, thom_f, thom_s, thom_w, thom_p, strch)
    do b = 1, n_groups
      do i = 1, SZ
        call der_univ_thom_per_line(i, b, n, du, u, n_tds, coeffs, alpha, &
                                    thom_f, thom_s, thom_w, thom_p, strch)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_univ_thom_per

  subroutine der_univ_thom_per_line(i, b, n, du, u, n_tds, coeffs, alpha, &
                                    thom_f, thom_s, thom_w, thom_p, strch)
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(in) :: u(SZ, n, *)
    real(dp), intent(in) :: coeffs(*)
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: thom_f(*), thom_s(*), thom_w(*), thom_p(*), &
                            strch(*)

    integer :: j, jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4
    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4
    real(dp) :: temp_du, ss

    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    ! j = 1 is peeled so its j - 1 does not wrap onto du(i, n_tds, b) before
    ! that has been written (see cuda/kernels/thomas.f90).
    temp_du = c_m4*u(i, n_tds - 3, b) + c_m3*u(i, n_tds - 2, b) &
              + c_m2*u(i, n_tds - 1, b) + c_m1*u(i, n_tds, b) &
              + c_j*u(i, 1, b) &
              + c_p1*u(i, 2, b) + c_p2*u(i, 3, b) &
              + c_p3*u(i, 4, b) + c_p4*u(i, 5, b)
    du(i, 1, b) = temp_du

    do j = 2, n_tds
      jm4 = modulo(j - 5, n_tds) + 1
      jm3 = modulo(j - 4, n_tds) + 1
      jm2 = modulo(j - 3, n_tds) + 1
      jm1 = j - 1
      jp1 = modulo(j - n_tds, n_tds) + 1
      jp2 = modulo(j - n_tds + 1, n_tds) + 1
      jp3 = modulo(j - n_tds + 2, n_tds) + 1
      jp4 = modulo(j - n_tds + 3, n_tds) + 1

      temp_du = c_m4*u(i, jm4, b) + c_m3*u(i, jm3, b) &
                + c_m2*u(i, jm2, b) + c_m1*u(i, jm1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, jp1, b) + c_p2*u(i, jp2, b) &
                + c_p3*u(i, jp3, b) + c_p4*u(i, jp4, b)
      du(i, j, b) = temp_du - du(i, jm1, b)*thom_s(j)
    end do

    ! Backward pass of the Thomas algorithm
    du(i, n_tds, b) = du(i, n_tds, b)*thom_w(n_tds)
    do j = n_tds - 1, 1, -1
      du(i, j, b) = (du(i, j, b) - thom_f(j)*du(i, j + 1, b))*thom_w(j)
    end do

    ! Periodic final pass
    ss = (du(i, 1, b) - alpha*du(i, n_tds, b)) &
         /(1._dp + thom_p(1) - alpha*thom_p(n_tds))
    do j = 1, n_tds
      du(i, j, b) = (du(i, j, b) - ss*thom_p(j))*strch(j)
    end do

  end subroutine der_univ_thom_per_line

end module m_omptgt_kernels_thom
