module m_exec_thom

  use m_common, only: dp
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t

  implicit none

  private
  public :: exec_thom_tds_compact

contains

  subroutine exec_thom_tds_compact(du, u, tdsops)

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    type(tdsops_t), intent(in) :: tdsops

    if (tdsops%periodic) then
      call der_univ_thom_per(du, u, tdsops)
    else
      call der_univ_thom(du, u, tdsops)
    end if
    
  end subroutine exec_thom_tds_compact

  subroutine der_univ_thom(du, u, tdsops)

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    type(tdsops_t), intent(in) :: tdsops

    integer :: i, j, k
    integer :: nk
    integer :: ctr

    nk = size(u, dim=3)
    associate(coeffs_s => tdsops%coeffs_s, coeffs_e => tdsops%coeffs_e, coeffs => tdsops%coeffs, &
              n => tdsops%tds_n, &
              thom_f => tdsops%thom_f, thom_s => tdsops%thom_s, thom_w => tdsops%thom_w)
      !$omp parallel do private(i, j, ctr)
      do k = 1, nk
        !! Forward pass
        ! XXX: Extract as PURE helper?

        !$omp simd
        do i = 1, SZ
          du(i, 1, k) = sum(coeffs_s(5:9, 1) * u(i, 1:5, k))
        end do
        !$omp end simd
        do j = 2, 4
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = sum(coeffs_s(4-(j-2):9, j) * u(i, 1:4+j, k)) &
              - du(i, j - 1, k) * thom_s(j)
          end do
          !$omp end simd
        end do
       
        do j = 5, n - 4
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = sum(coeffs(1:9) * u(i, j-4:j+4, k)) &
              - du(i, j - 1, k) * thom_s(j)
          end do
          !$omp end simd
        end do
       
        ctr = 0
        do j = n - 3, n
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = sum(coeffs_e(1:8-ctr, ctr + 1) * u(i, j-4:j+(3-ctr), k)) &
              - du(i, j - 1, k) * thom_s(j)
          end do
          !$omp end simd
          ctr = ctr + 1
        end do

        !! Backward pass
        ! XXX: Extract as PURE helper?

        !$omp simd
        do i = 1, SZ
          du(i, n, k) = du(i, n, k) * thom_w(n)
        end do
        !$omp end simd
        do j = n - 1, 1, -1
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = (du(i, j, k) - thom_f(j) * du(i, j + 1, k)) * thom_w(j)
          end do
          !$omp end simd
        end do
      end do
      !$omp end parallel do
    end associate
    
  end subroutine der_univ_thom

  subroutine der_univ_thom_per(du, u, tdsops)

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    type(tdsops_t), intent(in) :: tdsops

    integer :: nk
    integer :: i, j, k

    real(dp), dimension(SZ) :: ss
    integer, dimension(4) :: jm, jp

    nk = size(u, dim=3)
    associate(coeffs_s => tdsops%coeffs_s, coeffs_e => tdsops%coeffs_e, coeffs => tdsops%coeffs, &
              n => tdsops%tds_n, &
              alpha => tdsops%alpha, &
              thom_f => tdsops%thom_f, thom_s => tdsops%thom_s, thom_w => tdsops%thom_w, &
              thom_p => tdsops%thom_p)
      !$omp parallel do private(i, j, jm, jp, ss)
      do k = 1, nk
        !! Forward pass
        do j = 1, n
          do i = 1, 4
            jm(i) = modulo(j - (5 - (i - 1)), n) + 1
            jp(i) = modulo(j - (n - (i - 1)), n) + 1
          end do

          !$omp simd
          do i = 1, SZ
            du(i, j, k) = sum(coeffs(1:4) * u(i, jm(1:4), k))
            du(i, j, k) = du(i, j, k) + coeffs(5) * u(i, j, k)
            du(i, j, k) = du(i, j, k) + sum(coeffs(6:9) * u(i, jp(1:4), k))
            du(i, j, k) = du(i, j, k) - du(i, jm(4), k) * thom_s(j)
          end do
          !$omp end simd
        end do

        !! Backward pass

        !$omp simd
        do i = 1, SZ
          du(i, n, k) = du(i, n, k) * thom_w(n)
        end do
        !$omp end simd
        do j = n - 1, 1, -1
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = (du(i, j, k) - thom_f(j) * du(i, j + 1, k)) * thom_w(j)
          end do
          !$omp end simd
        end do

        !! Periodic final pass

        !$omp simd
        do i = 1, SZ
          ss(i) = (du(i, 1, k) - alpha * du(i, n, k)) &
            / (1.0_dp + thom_p(1) - alpha * thom_p(n))
        end do
        !$omp end simd
        do j = 1, n
          !$omp simd
          do i = 1, SZ
            du(i, j, k) = du(i, j, k) - ss(i) * thom_p(j)
          end do
          !$omp end simd
        end do
      end do
      !$omp end parallel do
    end associate
    
  end subroutine der_univ_thom_per
  
end module m_exec_thom
