module m_cuda_tdsops
  use iso_fortran_env, only: stderr => error_unit

  use m_common, only: dp
  use m_tdsops, only: tdsops_t, tdsops_init

  implicit none

  type, extends(tdsops_t) :: cuda_tdsops_t
    !! CUDA extension of the Tridiagonal Solver Operators class.
    !!
    !! Regular tdsops_t class is initiated and the coefficient arrays are
    !! copied into device arrays so that cuda kernels can use them.
    real(dp), device, allocatable :: dist_fw_dev(:), dist_bw_dev(:), &
                                     dist_sa_dev(:), dist_sc_dev(:), &
                                     dist_af_dev(:)
    real(dp), device, allocatable :: thom_f_dev(:), thom_s_dev(:), &
                                     thom_w_dev(:), thom_p_dev(:)
    real(dp), device, allocatable :: stretch_dev(:), stretch_correct_dev(:)
    real(dp), device, allocatable :: coeffs_dev(:), &
                                     coeffs_s_dev(:, :), coeffs_e_dev(:, :)
  contains
  end type cuda_tdsops_t

  interface cuda_tdsops_t
    module procedure cuda_tdsops_init
  end interface cuda_tdsops_t

contains

  function cuda_tdsops_init( &
    n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    ) result(tdsops)
    !! Constructor function for the cuda_tdsops_t class.
    !! See tdsops_t for details.
    implicit none

    type(cuda_tdsops_t) :: tdsops !! return value of the function

    integer, intent(in) :: n_tds
    real(dp), intent(in) :: delta
    character(*), intent(in) :: operation, scheme
    integer, intent(in) :: bc_start, bc_end
    real(dp), optional, intent(in) :: stretch(:), stretch_correct(:)
    integer, optional, intent(in) :: n_halo
    character(*), optional, intent(in) :: from_to
    logical, optional, intent(in) :: sym
    real(dp), optional, intent(in) :: c_nu, nu0_nu

    integer :: n, n_stencil

    tdsops%tdsops_t = tdsops_init(n_tds, delta, operation, scheme, bc_start, &
                                  bc_end, stretch, stretch_correct, n_halo, &
                                  from_to, sym, c_nu, nu0_nu)

    n = tdsops%n_rhs
    allocate (tdsops%dist_fw_dev(n), tdsops%dist_bw_dev(n))
    allocate (tdsops%dist_sa_dev(n), tdsops%dist_sc_dev(n))
    allocate (tdsops%dist_af_dev(n))
    allocate (tdsops%thom_f_dev(n), tdsops%thom_s_dev(n))
    allocate (tdsops%thom_w_dev(n), tdsops%thom_p_dev(n))

    allocate (tdsops%stretch_dev(tdsops%n_tds))
    allocate (tdsops%stretch_correct_dev(tdsops%n_tds))

    n_stencil = 2*tdsops%n_halo + 1
    allocate (tdsops%coeffs_dev(n_stencil))
    allocate (tdsops%coeffs_s_dev(n_stencil, tdsops%n_halo))
    allocate (tdsops%coeffs_e_dev(n_stencil, tdsops%n_halo))

    tdsops%dist_fw_dev(:) = tdsops%dist_fw(:)
    tdsops%dist_bw_dev(:) = tdsops%dist_bw(:)
    tdsops%dist_sa_dev(:) = tdsops%dist_sa(:)
    tdsops%dist_sc_dev(:) = tdsops%dist_sc(:)
    tdsops%dist_af_dev(:) = tdsops%dist_af(:)

    tdsops%thom_f_dev(:) = tdsops%thom_f(:)
    tdsops%thom_s_dev(:) = tdsops%thom_s(:)
    tdsops%thom_w_dev(:) = tdsops%thom_w(:)
    tdsops%thom_p_dev(:) = tdsops%thom_p(:)

    tdsops%stretch_dev(:) = tdsops%stretch(:)
    tdsops%stretch_correct_dev(:) = tdsops%stretch_correct(:)

    tdsops%coeffs_dev(:) = tdsops%coeffs(:)
    tdsops%coeffs_s_dev(:, :) = tdsops%coeffs_s(:, :)
    tdsops%coeffs_e_dev(:, :) = tdsops%coeffs_e(:, :)

  end function cuda_tdsops_init

end module m_cuda_tdsops

