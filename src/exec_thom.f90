module m_cuda_exec_thom
  use cudafor

  use m_common, only: dp
  use m_cuda_kernels_thom, only: der_univ_thom, der_univ_thom_per
  use m_cuda_tdsops, only: cuda_tdsops_t

  implicit none

contains

  subroutine exec_thom_tds_compact(du, u, tdsops, blocks, threads)
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: du
    real(dp), device, dimension(:, :, :), intent(in) :: u
    type(cuda_tdsops_t), intent(in) :: tdsops
    type(dim3), intent(in) :: blocks, threads

    if (tdsops%periodic) then
      call der_univ_thom_per<<<blocks, threads>>>( & !&
        du, u, tdsops%n_tds, tdsops%coeffs_dev, tdsops%alpha, &
        tdsops%thom_f_dev, tdsops%thom_s_dev, tdsops%thom_w_dev, &
        tdsops%thom_p_dev, tdsops%stretch_dev &
        )
    else
      call der_univ_thom<<<blocks, threads>>>( & !&
        du, u, &
        tdsops%n_tds, tdsops%n_rhs, &
        tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
        tdsops%thom_f_dev, tdsops%thom_s_dev, tdsops%thom_w_dev, &
        tdsops%stretch_dev &
        )
    end if

  end subroutine exec_thom_tds_compact

end module m_cuda_exec_thom
