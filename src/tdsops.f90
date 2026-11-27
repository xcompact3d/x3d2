module m_omptgt_tdsops
  !! OpenMP target extension of the tridiagonal solver operators.
  !!
  !! The host `tdsops_t` is built as usual and its coefficient arrays are
  !! copied once into device arrays, as `cuda_tdsops_t` does, so the kernels
  !! read them in place instead of mapping them on every launch.
  !!
  !! The device arrays are not freed: there is no finaliser, because the
  !! intrinsic assignment of a structure constructor would finalise the
  !! temporary and release storage the copy still points at. Operators are
  !! built once at start-up and live for the whole run, so construct them in
  !! place with `omptgt_tdsops_init` rather than by assignment.
  use m_common, only: dp
  use m_tdsops, only: tdsops_t, tdsops_init

  use m_omptgt_common, only: dev_alloc, dev_clone, dev_upload

  implicit none

  type, extends(tdsops_t) :: omptgt_tdsops_t
    real(dp), pointer, contiguous :: dist_fw_dev(:) => null(), &
                                     dist_bw_dev(:) => null(), &
                                     dist_sa_dev(:) => null(), &
                                     dist_sc_dev(:) => null(), &
                                     dist_af_dev(:) => null()
    real(dp), pointer, contiguous :: thom_f_dev(:) => null(), &
                                     thom_s_dev(:) => null(), &
                                     thom_w_dev(:) => null(), &
                                     thom_p_dev(:) => null()
    real(dp), pointer, contiguous :: stretch_dev(:) => null(), &
                                     stretch_correct_dev(:) => null()
    real(dp), pointer, contiguous :: coeffs_dev(:) => null(), &
                                     coeffs_s_dev(:, :) => null(), &
                                     coeffs_e_dev(:, :) => null()
    !> Columns of A_np^{-1}*[e_1|e_2|e_{n-1}|e_n] for the cyclic
    !! pentadiagonal solve. They depend on the operator alone, so they are
    !! built once here where the CUDA kernel rebuilds them in every block.
    real(dp), pointer, contiguous :: penta_z_dev(:, :) => null()
  end type omptgt_tdsops_t

contains

  subroutine omptgt_tdsops_init( &
    tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu, &
    filter_alpha &
    )
    !! Initialises `tdsops` in place. See tdsops_t for the arguments.
    implicit none

    type(omptgt_tdsops_t), intent(inout) :: tdsops
    integer, intent(in) :: n_tds
    real(dp), intent(in) :: delta
    character(*), intent(in) :: operation, scheme
    integer, intent(in) :: bc_start, bc_end
    real(dp), optional, intent(in) :: stretch(:), stretch_correct(:)
    integer, optional, intent(in) :: n_halo
    character(*), optional, intent(in) :: from_to
    logical, optional, intent(in) :: sym
    real(dp), optional, intent(in) :: c_nu, nu0_nu
    real(dp), optional, intent(in) :: filter_alpha

    tdsops%tdsops_t = tdsops_init(n_tds, delta, operation, scheme, bc_start, &
                                  bc_end, stretch, stretch_correct, n_halo, &
                                  from_to, sym, c_nu, nu0_nu, filter_alpha)

    ! Each copy takes its host array's own size: preprocessing may
    ! reallocate an array away from the n_rhs it was first given.
    call dev_clone(tdsops%dist_fw_dev, tdsops%dist_fw)
    call dev_clone(tdsops%dist_bw_dev, tdsops%dist_bw)
    call dev_clone(tdsops%dist_sa_dev, tdsops%dist_sa)
    call dev_clone(tdsops%dist_sc_dev, tdsops%dist_sc)
    call dev_clone(tdsops%dist_af_dev, tdsops%dist_af)
    call dev_clone(tdsops%thom_f_dev, tdsops%thom_f)
    call dev_clone(tdsops%thom_s_dev, tdsops%thom_s)
    call dev_clone(tdsops%thom_w_dev, tdsops%thom_w)
    call dev_clone(tdsops%thom_p_dev, tdsops%thom_p)
    call dev_clone(tdsops%stretch_dev, tdsops%stretch)
    call dev_clone(tdsops%stretch_correct_dev, tdsops%stretch_correct)
    call dev_clone(tdsops%coeffs_dev, tdsops%coeffs)
    call dev_clone(tdsops%coeffs_s_dev, tdsops%coeffs_s)
    call dev_clone(tdsops%coeffs_e_dev, tdsops%coeffs_e)

    if (tdsops%pentadiag .and. tdsops%periodic) then
      call dev_alloc(tdsops%penta_z_dev, tdsops%n_tds, 4)
      call dev_upload(tdsops%penta_z_dev, penta_periodic_columns(tdsops))
    end if

  end subroutine omptgt_tdsops_init

  function penta_periodic_columns(tdsops) result(z)
    !! Solves A_np*z_k = e_{p_k} for p = 1, 2, n-1, n with the non-periodic
    !! pentadiagonal LU factors, as thread 1 of der_penta_periodic does in
    !! the CUDA backend.
    implicit none

    class(tdsops_t), intent(in) :: tdsops
    real(dp) :: z(tdsops%n_tds, 4)

    integer :: pos(4), k, j, n

    n = tdsops%n_tds
    pos = [1, 2, n - 1, n]
    associate (ffr => tdsops%dist_fw, faf => tdsops%dist_af, &
               fsa => tdsops%dist_sa, fbw => tdsops%dist_bw, &
               bet => tdsops%beta, bet_s => tdsops%beta_lhs_s)
      do k = 1, 4
        z(:, k) = 0._dp
        z(pos(k), k) = 1._dp
        z(2, k) = z(2, k) - faf(2)*z(1, k)
        do j = 3, n
          z(j, k) = z(j, k) - faf(j)*z(j - 1, k) - fsa(j)*z(j - 2, k)
        end do
        z(n, k) = z(n, k)*ffr(n)
        z(n - 1, k) = (z(n - 1, k) - fbw(n - 1)*z(n, k))*ffr(n - 1)
        do j = n - 2, 2, -1
          z(j, k) = (z(j, k) - fbw(j)*z(j + 1, k) - bet*z(j + 2, k))*ffr(j)
        end do
        z(1, k) = (z(1, k) - fbw(1)*z(2, k) - bet_s*z(3, k))*ffr(1)
      end do
    end associate

  end function penta_periodic_columns

end module m_omptgt_tdsops
