module m_omptgt_exec_thom
  use m_common, only: dp
  use m_omptgt_kernels_thom, only: der_univ_thom, der_univ_thom_per
  use m_omptgt_tdsops, only: omptgt_tdsops_t

  implicit none

contains

  subroutine exec_thom_tds_compact(du, u, tdsops, n_groups)
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: du, u
    type(omptgt_tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n_groups

    if (tdsops%periodic) then
      call der_univ_thom_per( &
        du, u, tdsops%n_tds, tdsops%coeffs_dev, tdsops%alpha, &
        tdsops%thom_f_dev, tdsops%thom_s_dev, tdsops%thom_w_dev, &
        tdsops%thom_p_dev, tdsops%stretch_dev, size(u, 2), n_groups &
        )
    else
      call der_univ_thom( &
        du, u, tdsops%n_tds, tdsops%n_rhs, &
        tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
        tdsops%thom_f_dev, tdsops%thom_s_dev, tdsops%thom_w_dev, &
        tdsops%stretch_dev, size(u, 2), n_groups &
        )
    end if

  end subroutine exec_thom_tds_compact

end module m_omptgt_exec_thom
