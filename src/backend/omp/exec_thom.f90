module m_exec_thom
  !! Local Thomas algorithm execution for OMP backend.
  !!
  !! Provides parallel execution of compact finite difference schemes using
  !! standard Thomas algorithm (tridiagonal solver). Used when domain is not
  !! decomposed in the derivative direction (all data local to process).
  !!
  !! **Two variants:**
  !!
  !! - **Non-periodic:** Standard Thomas with arbitrary boundary conditions
  !! - **Periodic:** Modified Thomas for cyclic tridiagonal systems
  !!
  !! **Parallelisation:** OpenMP over pencil groups (no MPI needed)
  !!
  !! **Contrast with distributed:** exec_dist handles multi-process case,
  !! this module handles single-process-per-direction case.

  use m_common, only: dp
  use m_tdsops, only: tdsops_t

  use m_omp_kernels_thom, only: der_univ_thom, der_univ_thom_per

  implicit none

  private
  public :: exec_thom_tds_compact

contains

  subroutine exec_thom_tds_compact(du, u, tdsops, n_groups)
    !! Execute local Thomas algorithm for compact scheme.
    !!
    !! Applies compact finite difference operator using tridiagonal solver.
    !! Chooses periodic or non-periodic variant based on operator configuration.
    !! All computation local to process (no MPI communication).
    !!
    !! **Algorithm selection:**
    !! - `periodic=.true.`: Sherman-Morrison formula for cyclic system
    !! - `periodic=.false.`: Standard forward/backward Thomas algorithm
    !!
    !! **Parallelisation:** OpenMP parallel loop over pencil groups

    real(dp), dimension(:, :, :), intent(out) :: du  !! Derivative output
    real(dp), dimension(:, :, :), intent(in) :: u    !! Input field
    type(tdsops_t), intent(in) :: tdsops             !! Compact scheme operator
    integer, intent(in) :: n_groups                  !! Number of pencil groups

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
