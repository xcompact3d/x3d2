module m_exec_thom

  use m_common, only: dp
  use m_tdsops, only: tdsops_t

  use m_omp_kernels_thom, only: der_univ_thom, der_univ_thom_per

  implicit none

  private
  public :: exec_thom_tds_compact

contains

  subroutine exec_thom_tds_compact(du, u, tdsops, n_groups)

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    type(tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n_groups

    integer :: k

    if (tdsops%periodic) then
      !$omp parallel do
      do k = 1, n_groups
        call der_univ_thom_per( &
          du(:, :, k), u(:, :, k), tdsops%n_tds, tdsops%coeffs, tdsops%alpha, &
          tdsops%thom_f, tdsops%thom_s, tdsops%thom_w, &
          tdsops%thom_p, tdsops%stretch &
          )
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do k = 1, n_groups
        call der_univ_thom( &
          du(:, :, k), u(:, :, k), tdsops%n_tds, tdsops%n_rhs, &
          tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, &
          tdsops%thom_f, tdsops%thom_s, tdsops%thom_w, &
          tdsops%stretch &
          )
      end do
      !$omp end parallel do
    end if

  end subroutine exec_thom_tds_compact

end module m_exec_thom
