module m_omptgt_kernels_dist
  !! Distributed compact-scheme kernels for the OpenMP target backend,
  !! following cuda/kernels/distributed.f90.
  !!
  !! A CUDA kernel there runs one thread per pencil line, (i, b) = (threadIdx,
  !! blockIdx), each sweeping its line serially. Here every kernel offloads
  !! the same (b, i) iteration space with `teams distribute parallel do` and
  !! hands each line to a `declare target` routine ending in `_line`. Keeping
  !! the sweep in its own routine makes all of its temporaries per-line locals
  !! with no private clause to keep in step, and the explicit `distribute
  !! parallel do` keeps NVHPC from serialising each team, as it does when
  !! `teams loop` meets a procedure call.
  !!
  !! Arrays are explicit-shape or assumed-size dummies named in
  !! `OMP_DEVICE_ADDR` (see m_omptgt_common). Fields are (SZ, n, n_groups),
  !! halo buffers (SZ, 4, n_groups) and the 2x2-system buffers
  !! (SZ, 1, n_groups). The boundary stencils are 9 wide (n_halo = 4).
  use m_common, only: dp
  use m_omptgt_common, only: SZ

  implicit none

  private
  public :: der_univ_dist, der_univ_subs, transeq_3fused_dist, &
            transeq_3fused_subs, der_penta_full, der_penta_periodic

contains

  ! --------------------------------------------------------------------------
  ! Boundary stencil helpers. Row j of the first four uses coeffs_s(:, j) over
  ! points j-4..j+4, those before 1 coming from the start halo; row r of the
  ! last four (j = n_rhs - 4 + r) uses coeffs_e(:, r), those past n_rhs from
  ! the end halo. The sum runs k = 1..9 like the unrolled CUDA expressions.
  ! --------------------------------------------------------------------------

  pure real(dp) function stencil_s(cs, j, i, b, n, u, u_s) result(r)
    !$omp declare target
    real(dp), intent(in) :: cs(9, *)
    integer, intent(in) :: j, i, b, n
    real(dp), intent(in) :: u(SZ, n, *), u_s(SZ, 4, *)

    integer :: k, m

    r = 0._dp
    do k = 1, 9
      m = j - 5 + k
      if (m >= 1) then
        r = r + cs(k, j)*u(i, m, b)
      else
        r = r + cs(k, j)*u_s(i, 4 + m, b)
      end if
    end do
  end function stencil_s

  pure real(dp) function stencil_e(ce, row, n_rhs, i, b, n, u, u_e) result(r)
    !$omp declare target
    real(dp), intent(in) :: ce(9, *)
    integer, intent(in) :: row, n_rhs, i, b, n
    real(dp), intent(in) :: u(SZ, n, *), u_e(SZ, 4, *)

    integer :: k, m

    r = 0._dp
    do k = 1, 9
      m = n_rhs - 4 + row - 5 + k
      if (m <= n_rhs) then
        r = r + ce(k, row)*u(i, m, b)
      else
        r = r + ce(k, row)*u_e(i, m - n_rhs, b)
      end if
    end do
  end function stencil_e

  pure real(dp) function stencil_s_prod(cs, j, i, b, n, u, u_s, v, v_s) &
    result(r)
    !$omp declare target
    !! As stencil_s, applied to the product u*v.
    real(dp), intent(in) :: cs(9, *)
    integer, intent(in) :: j, i, b, n
    real(dp), intent(in) :: u(SZ, n, *), u_s(SZ, 4, *)
    real(dp), intent(in) :: v(SZ, n, *), v_s(SZ, 4, *)

    integer :: k, m

    r = 0._dp
    do k = 1, 9
      m = j - 5 + k
      if (m >= 1) then
        r = r + cs(k, j)*u(i, m, b)*v(i, m, b)
      else
        r = r + cs(k, j)*u_s(i, 4 + m, b)*v_s(i, 4 + m, b)
      end if
    end do
  end function stencil_s_prod

  pure real(dp) function stencil_e_prod(ce, row, n_rhs, i, b, n, u, u_e, &
                                        v, v_e) result(r)
    !$omp declare target
    !! As stencil_e, applied to the product u*v.
    real(dp), intent(in) :: ce(9, *)
    integer, intent(in) :: row, n_rhs, i, b, n
    real(dp), intent(in) :: u(SZ, n, *), u_e(SZ, 4, *)
    real(dp), intent(in) :: v(SZ, n, *), v_e(SZ, 4, *)

    integer :: k, m

    r = 0._dp
    do k = 1, 9
      m = n_rhs - 4 + row - 5 + k
      if (m <= n_rhs) then
        r = r + ce(k, row)*u(i, m, b)*v(i, m, b)
      else
        r = r + ce(k, row)*u_e(i, m - n_rhs, b)*v_e(i, m - n_rhs, b)
      end if
    end do
  end function stencil_e_prod

  ! --------------------------------------------------------------------------
  ! der_univ_dist / der_univ_subs
  ! --------------------------------------------------------------------------

  subroutine der_univ_dist( &
    du, send_u_s, send_u_e, u, u_s, u_e, n_tds, n_rhs, &
    coeffs_s, coeffs_e, coeffs, ffr, fbc, faf, n, n_groups &
    )
    !! First phase of the distributed solve: builds the RHS and runs the
    !! forward and backward hybrid sweeps, leaving the 2x2 boundary values in
    !! the send buffers.
    integer, intent(in) :: n_tds, n_rhs, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(inout) :: send_u_s(SZ, 1, n_groups), &
                               send_u_e(SZ, 1, n_groups)
    real(dp), intent(in) :: u(SZ, n, n_groups), u_s(SZ, 4, n_groups), &
                            u_e(SZ, 4, n_groups)
    real(dp), intent(in) :: coeffs_s(9, 4), coeffs_e(9, 4), coeffs(9)
    real(dp), intent(in) :: ffr(n_rhs), fbc(n_rhs), faf(n_rhs)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, send_u_s, send_u_e, u, u_s, u_e, coeffs_s, &
    !$omp                 coeffs_e, coeffs, ffr, fbc, faf)
    do b = 1, n_groups
      do i = 1, SZ
        call der_univ_dist_line(i, b, n, du, send_u_s, send_u_e, u, u_s, &
                                u_e, n_tds, n_rhs, coeffs_s, coeffs_e, &
                                coeffs, ffr, fbc, faf)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_univ_dist

  subroutine der_univ_dist_line(i, b, n, du, send_u_s, send_u_e, u, u_s, &
                                u_e, n_tds, n_rhs, coeffs_s, coeffs_e, &
                                coeffs, ffr, fbc, faf)
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds, n_rhs
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(inout) :: send_u_s(SZ, 1, *), send_u_e(SZ, 1, *)
    real(dp), intent(in) :: u(SZ, n, *), u_s(SZ, 4, *), u_e(SZ, 4, *)
    real(dp), intent(in) :: coeffs_s(9, *), coeffs_e(9, *), coeffs(*)
    real(dp), intent(in) :: ffr(*), fbc(*), faf(*)

    integer :: j, row
    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, &
                temp_du, alpha, last_r

    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)
    last_r = ffr(1)

    du(i, 1, b) = stencil_s(coeffs_s, 1, i, b, n, u, u_s)*faf(1)
    du(i, 2, b) = stencil_s(coeffs_s, 2, i, b, n, u, u_s)*faf(2)
    do j = 3, 4
      du(i, j, b) = ffr(j)*(stencil_s(coeffs_s, j, i, b, n, u, u_s) &
                            - faf(j)*du(i, j - 1, b))
    end do

    alpha = faf(5)
    do j = 5, n_rhs - 4
      temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
      du(i, j, b) = ffr(j)*(temp_du - alpha*du(i, j - 1, b))
    end do

    do row = 1, 4
      j = n_rhs - 4 + row
      du(i, j, b) = ffr(j)*(stencil_e(coeffs_e, row, n_rhs, i, b, n, u, u_e) &
                            - faf(j)*du(i, j - 1, b))
    end do

    send_u_e(i, 1, b) = du(i, n_tds, b)

    ! Backward pass of the hybrid algorithm
    do j = n_tds - 2, 2, -1
      du(i, j, b) = du(i, j, b) - fbc(j)*du(i, j + 1, b)
    end do
    du(i, 1, b) = last_r*(du(i, 1, b) - fbc(1)*du(i, 2, b))
    send_u_s(i, 1, b) = du(i, 1, b)

  end subroutine der_univ_dist_line

  subroutine der_univ_subs(du, recv_u_s, recv_u_e, n_tds, dist_sa, dist_sc, &
                           strch, n, n_groups)
    !! Second phase of the distributed solve: resolves the 2x2 boundary
    !! systems with the neighbours' values and substitutes back.
    integer, intent(in) :: n_tds, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(in) :: recv_u_s(SZ, 1, n_groups), &
                            recv_u_e(SZ, 1, n_groups)
    real(dp), intent(in) :: dist_sa(n_tds), dist_sc(n_tds), strch(n_tds)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, recv_u_s, recv_u_e, dist_sa, dist_sc, strch)
    do b = 1, n_groups
      do i = 1, SZ
        call der_univ_subs_line(i, b, n, du, recv_u_s, recv_u_e, n_tds, &
                                dist_sa, dist_sc, strch)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_univ_subs

  subroutine der_univ_subs_line(i, b, n, du, recv_u_s, recv_u_e, n_tds, &
                                dist_sa, dist_sc, strch)
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(in) :: recv_u_s(SZ, 1, *), recv_u_e(SZ, 1, *)
    real(dp), intent(in) :: dist_sa(*), dist_sc(*), strch(*)

    integer :: j
    real(dp) :: ur, bl, recp, du_s, du_e

    ! The 2x2 trick is valid for symmetric Toeplitz matrices, which ours are
    ! in the (5:n-4) region; see cuda/kernels/distributed.f90.
    ! At the start we have the 'bl', and assume 'ur'
    bl = dist_sa(1)
    ur = dist_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    du_s = recp*(du(i, 1, b) - bl*recv_u_s(i, 1, b))

    ! At the end we have the 'ur', and assume 'bl'
    bl = dist_sc(n_tds)
    ur = dist_sc(n_tds)
    recp = 1._dp/(1._dp - ur*bl)
    du_e = recp*(du(i, n_tds, b) - ur*recv_u_e(i, 1, b))

    du(i, 1, b) = du_s*strch(1)
    do j = 2, n_tds - 1
      du(i, j, b) = (du(i, j, b) - dist_sa(j)*du_s - dist_sc(j)*du_e)*strch(j)
    end do
    du(i, n_tds, b) = du_e*strch(n_tds)

  end subroutine der_univ_subs_line

  ! --------------------------------------------------------------------------
  ! Fused transport equation kernels
  ! --------------------------------------------------------------------------

  subroutine transeq_3fused_dist( &
    du, dud, d2u, send_du_s, send_du_e, send_dud_s, send_dud_e, &
    send_d2u_s, send_d2u_e, u, u_s, u_e, v, v_s, v_e, n_tds, n_rhs, &
    du_coeffs_s, du_coeffs_e, du_coeffs, du_fw, du_bw, du_af, &
    dud_coeffs_s, dud_coeffs_e, dud_coeffs, dud_fw, dud_bw, dud_af, &
    d2u_coeffs_s, d2u_coeffs_e, d2u_coeffs, d2u_fw, d2u_bw, d2u_af, &
    n, n_groups &
    )
    !! First phase for du = d(u), dud = d(u*v) and d2u = d2(u) in one sweep.
    integer, intent(in) :: n_tds, n_rhs, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups), dud(SZ, n, n_groups), &
                               d2u(SZ, n, n_groups)
    real(dp), intent(inout) :: send_du_s(SZ, 1, n_groups), &
                               send_du_e(SZ, 1, n_groups), &
                               send_dud_s(SZ, 1, *), send_dud_e(SZ, 1, *), &
                               send_d2u_s(SZ, 1, *), send_d2u_e(SZ, 1, *)
    real(dp), intent(in) :: u(SZ, n, n_groups), u_s(SZ, 4, n_groups), &
                            u_e(SZ, 4, n_groups)
    real(dp), intent(in) :: v(SZ, n, n_groups), v_s(SZ, 4, n_groups), &
                            v_e(SZ, 4, n_groups)
    real(dp), intent(in) :: du_coeffs_s(9, 4), du_coeffs_e(9, 4), &
                            du_coeffs(9), du_fw(n_rhs), du_bw(n_rhs), &
                            du_af(n_rhs)
    real(dp), intent(in) :: dud_coeffs_s(9, 4), dud_coeffs_e(9, 4), &
                            dud_coeffs(9), dud_fw(n_rhs), dud_bw(n_rhs), &
                            dud_af(n_rhs)
    real(dp), intent(in) :: d2u_coeffs_s(9, 4), d2u_coeffs_e(9, 4), &
                            d2u_coeffs(9), d2u_fw(n_rhs), d2u_bw(n_rhs), &
                            d2u_af(n_rhs)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, dud, d2u, send_du_s, send_du_e, send_dud_s, &
    !$omp                 send_dud_e, send_d2u_s, send_d2u_e, u, u_s, u_e, &
    !$omp                 v, v_s, v_e, du_coeffs_s, du_coeffs_e, du_coeffs, &
    !$omp                 du_fw, du_bw, du_af, dud_coeffs_s, dud_coeffs_e, &
    !$omp                 dud_coeffs, dud_fw, dud_bw, dud_af, d2u_coeffs_s, &
    !$omp                 d2u_coeffs_e, d2u_coeffs, d2u_fw, d2u_bw, d2u_af)
    do b = 1, n_groups
      do i = 1, SZ
        call transeq_3fused_dist_line( &
          i, b, n, du, dud, d2u, send_du_s, send_du_e, send_dud_s, &
          send_dud_e, send_d2u_s, send_d2u_e, u, u_s, u_e, v, v_s, v_e, &
          n_tds, n_rhs, &
          du_coeffs_s, du_coeffs_e, du_coeffs, du_fw, du_bw, du_af, &
          dud_coeffs_s, dud_coeffs_e, dud_coeffs, dud_fw, dud_bw, dud_af, &
          d2u_coeffs_s, d2u_coeffs_e, d2u_coeffs, d2u_fw, d2u_bw, d2u_af)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine transeq_3fused_dist

  subroutine transeq_3fused_dist_line( &
    i, b, n, du, dud, d2u, send_du_s, send_du_e, send_dud_s, send_dud_e, &
    send_d2u_s, send_d2u_e, u, u_s, u_e, v, v_s, v_e, n_tds, n_rhs, &
    du_coeffs_s, du_coeffs_e, du_coeffs, du_fw, du_bw, du_af, &
    dud_coeffs_s, dud_coeffs_e, dud_coeffs, dud_fw, dud_bw, dud_af, &
    d2u_coeffs_s, d2u_coeffs_e, d2u_coeffs, d2u_fw, d2u_bw, d2u_af &
    )
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds, n_rhs
    real(dp), intent(inout) :: du(SZ, n, *), dud(SZ, n, *), d2u(SZ, n, *)
    real(dp), intent(inout) :: send_du_s(SZ, 1, *), send_du_e(SZ, 1, *), &
                               send_dud_s(SZ, 1, *), send_dud_e(SZ, 1, *), &
                               send_d2u_s(SZ, 1, *), send_d2u_e(SZ, 1, *)
    real(dp), intent(in) :: u(SZ, n, *), u_s(SZ, 4, *), u_e(SZ, 4, *)
    real(dp), intent(in) :: v(SZ, n, *), v_s(SZ, 4, *), v_e(SZ, 4, *)
    real(dp), intent(in) :: du_coeffs_s(9, *), du_coeffs_e(9, *), &
                            du_coeffs(*), du_fw(*), du_bw(*), du_af(*)
    real(dp), intent(in) :: dud_coeffs_s(9, *), dud_coeffs_e(9, *), &
                            dud_coeffs(*), dud_fw(*), dud_bw(*), dud_af(*)
    real(dp), intent(in) :: d2u_coeffs_s(9, *), d2u_coeffs_e(9, *), &
                            d2u_coeffs(*), d2u_fw(*), d2u_bw(*), d2u_af(*)

    integer :: j, row
    real(dp) :: du_c_m4, du_c_m3, du_c_m2, du_c_m1, du_c_j, &
                du_c_p1, du_c_p2, du_c_p3, du_c_p4, du_alpha
    real(dp) :: dud_c_m4, dud_c_m3, dud_c_m2, dud_c_m1, dud_c_j, &
                dud_c_p1, dud_c_p2, dud_c_p3, dud_c_p4, dud_alpha
    real(dp) :: d2u_c_m4, d2u_c_m3, d2u_c_m2, d2u_c_m1, d2u_c_j, &
                d2u_c_p1, d2u_c_p2, d2u_c_p3, d2u_c_p4, d2u_alpha
    real(dp) :: temp_du, temp_dud, temp_d2u
    real(dp) :: u_m4, u_m3, u_m2, u_m1, u_j, u_p1, u_p2, u_p3, u_p4
    real(dp) :: v_m4, v_m3, v_m2, v_m1, v_j, v_p1, v_p2, v_p3, v_p4
    real(dp) :: old_du, old_dud, old_d2u

    ! j = 1, 2
    do j = 1, 2
      du(i, j, b) = stencil_s(du_coeffs_s, j, i, b, n, u, u_s)*du_af(j)
      dud(i, j, b) = stencil_s_prod(dud_coeffs_s, j, i, b, n, u, u_s, &
                                    v, v_s)*dud_af(j)
      d2u(i, j, b) = stencil_s(d2u_coeffs_s, j, i, b, n, u, u_s)*d2u_af(j)
    end do
    ! j = 3, 4
    do j = 3, 4
      du(i, j, b) = du_fw(j)*(stencil_s(du_coeffs_s, j, i, b, n, u, u_s) &
                              - du_af(j)*du(i, j - 1, b))
      dud(i, j, b) = dud_fw(j)*(stencil_s_prod(dud_coeffs_s, j, i, b, n, &
                                               u, u_s, v, v_s) &
                                - dud_af(j)*dud(i, j - 1, b))
      d2u(i, j, b) = d2u_fw(j)*(stencil_s(d2u_coeffs_s, j, i, b, n, u, u_s) &
                                - d2u_af(j)*d2u(i, j - 1, b))
    end do

    du_alpha = du_af(5)
    dud_alpha = dud_af(5)
    d2u_alpha = d2u_af(5)

    du_c_m4 = du_coeffs(1); du_c_m3 = du_coeffs(2)
    du_c_m2 = du_coeffs(3); du_c_m1 = du_coeffs(4)
    du_c_j = du_coeffs(5)
    du_c_p1 = du_coeffs(6); du_c_p2 = du_coeffs(7)
    du_c_p3 = du_coeffs(8); du_c_p4 = du_coeffs(9)

    dud_c_m4 = dud_coeffs(1); dud_c_m3 = dud_coeffs(2)
    dud_c_m2 = dud_coeffs(3); dud_c_m1 = dud_coeffs(4)
    dud_c_j = dud_coeffs(5)
    dud_c_p1 = dud_coeffs(6); dud_c_p2 = dud_coeffs(7)
    dud_c_p3 = dud_coeffs(8); dud_c_p4 = dud_coeffs(9)

    d2u_c_m4 = d2u_coeffs(1); d2u_c_m3 = d2u_coeffs(2)
    d2u_c_m2 = d2u_coeffs(3); d2u_c_m1 = d2u_coeffs(4)
    d2u_c_j = d2u_coeffs(5)
    d2u_c_p1 = d2u_coeffs(6); d2u_c_p2 = d2u_coeffs(7)
    d2u_c_p3 = d2u_coeffs(8); d2u_c_p4 = d2u_coeffs(9)

    ! Carry d?(i, j - 1, b) and the stencils in registers through the bulk
    old_du = du(i, 4, b)
    old_dud = dud(i, 4, b)
    old_d2u = d2u(i, 4, b)

    u_m4 = u(i, 1, b); u_m3 = u(i, 2, b)
    u_m2 = u(i, 3, b); u_m1 = u(i, 4, b)
    u_j = u(i, 5, b); u_p1 = u(i, 6, b)
    u_p2 = u(i, 7, b); u_p3 = u(i, 8, b)
    v_m4 = v(i, 1, b); v_m3 = v(i, 2, b)
    v_m2 = v(i, 3, b); v_m1 = v(i, 4, b)
    v_j = v(i, 5, b); v_p1 = v(i, 6, b)
    v_p2 = v(i, 7, b); v_p3 = v(i, 8, b)

    do j = 5, n_rhs - 4
      u_p4 = u(i, j + 4, b); v_p4 = v(i, j + 4, b)

      temp_du = du_c_m4*u_m4 + du_c_m3*u_m3 + du_c_m2*u_m2 + du_c_m1*u_m1 &
                + du_c_j*u_j &
                + du_c_p1*u_p1 + du_c_p2*u_p2 + du_c_p3*u_p3 + du_c_p4*u_p4
      du(i, j, b) = du_fw(j)*(temp_du - du_alpha*old_du)
      old_du = du(i, j, b)

      temp_dud = dud_c_m4*u_m4*v_m4 + dud_c_m3*u_m3*v_m3 &
                 + dud_c_m2*u_m2*v_m2 + dud_c_m1*u_m1*v_m1 &
                 + dud_c_j*u_j*v_j &
                 + dud_c_p1*u_p1*v_p1 + dud_c_p2*u_p2*v_p2 &
                 + dud_c_p3*u_p3*v_p3 + dud_c_p4*u_p4*v_p4
      dud(i, j, b) = dud_fw(j)*(temp_dud - dud_alpha*old_dud)
      old_dud = dud(i, j, b)

      temp_d2u = d2u_c_m4*u_m4 + d2u_c_m3*u_m3 &
                 + d2u_c_m2*u_m2 + d2u_c_m1*u_m1 &
                 + d2u_c_j*u_j &
                 + d2u_c_p1*u_p1 + d2u_c_p2*u_p2 &
                 + d2u_c_p3*u_p3 + d2u_c_p4*u_p4
      d2u(i, j, b) = d2u_fw(j)*(temp_d2u - d2u_alpha*old_d2u)
      old_d2u = d2u(i, j, b)

      u_m4 = u_m3; u_m3 = u_m2; u_m2 = u_m1; u_m1 = u_j
      u_j = u_p1; u_p1 = u_p2; u_p2 = u_p3; u_p3 = u_p4
      v_m4 = v_m3; v_m3 = v_m2; v_m2 = v_m1; v_m1 = v_j
      v_j = v_p1; v_p1 = v_p2; v_p2 = v_p3; v_p3 = v_p4
    end do

    do row = 1, 4
      j = n_rhs - 4 + row
      du(i, j, b) = du_fw(j)*(stencil_e(du_coeffs_e, row, n_rhs, i, b, n, &
                                        u, u_e) &
                              - du_af(j)*du(i, j - 1, b))
      dud(i, j, b) = dud_fw(j)*(stencil_e_prod(dud_coeffs_e, row, n_rhs, &
                                               i, b, n, u, u_e, v, v_e) &
                                - dud_af(j)*dud(i, j - 1, b))
      d2u(i, j, b) = d2u_fw(j)*(stencil_e(d2u_coeffs_e, row, n_rhs, i, b, &
                                          n, u, u_e) &
                                - d2u_af(j)*d2u(i, j - 1, b))
    end do

    send_du_e(i, 1, b) = du(i, n_tds, b)
    send_dud_e(i, 1, b) = dud(i, n_tds, b)
    send_d2u_e(i, 1, b) = d2u(i, n_tds, b)

    ! Backward pass of the hybrid algorithm
    do j = n_tds - 2, 2, -1
      du(i, j, b) = du(i, j, b) - du_bw(j)*du(i, j + 1, b)
      dud(i, j, b) = dud(i, j, b) - dud_bw(j)*dud(i, j + 1, b)
      d2u(i, j, b) = d2u(i, j, b) - d2u_bw(j)*d2u(i, j + 1, b)
    end do
    du(i, 1, b) = du_fw(1)*(du(i, 1, b) - du_bw(1)*du(i, 2, b))
    dud(i, 1, b) = dud_fw(1)*(dud(i, 1, b) - dud_bw(1)*dud(i, 2, b))
    d2u(i, 1, b) = d2u_fw(1)*(d2u(i, 1, b) - d2u_bw(1)*d2u(i, 2, b))

    send_du_s(i, 1, b) = du(i, 1, b)
    send_dud_s(i, 1, b) = dud(i, 1, b)
    send_d2u_s(i, 1, b) = d2u(i, 1, b)

  end subroutine transeq_3fused_dist_line

  subroutine transeq_3fused_subs( &
    r_du, conv, dud, d2u, &
    recv_du_s, recv_du_e, recv_dud_s, recv_dud_e, recv_d2u_s, recv_d2u_e, &
    n_tds, nu, du_sa, du_sc, du_strch, dud_sa, dud_sc, dud_strch, &
    d2u_sa, d2u_sc, d2u_strch, d2u_strch_cor, n, n_groups &
    )
    !! Second phase: substitutes back and combines
    !! r_du = -1/2*(conv*du + dud) + nu*d2u.
    integer, intent(in) :: n_tds, n, n_groups
    real(dp), intent(in) :: nu
    !> Holds du on entry and the result on exit
    real(dp), intent(inout) :: r_du(SZ, n, n_groups)
    real(dp), intent(in) :: conv(SZ, n, n_groups), dud(SZ, n, n_groups), &
                            d2u(SZ, n, n_groups)
    real(dp), intent(in) :: recv_du_s(SZ, 1, n_groups), &
                            recv_du_e(SZ, 1, n_groups), recv_dud_s(SZ, 1, *), &
                            recv_dud_e(SZ, 1, *), recv_d2u_s(SZ, 1, *), &
                            recv_d2u_e(SZ, 1, *)
    real(dp), intent(in) :: du_sa(n_tds), du_sc(n_tds), du_strch(n_tds), &
                            dud_sa(n_tds), dud_sc(n_tds), dud_strch(n_tds), &
                            d2u_sa(n_tds), d2u_sc(n_tds), d2u_strch(n_tds), &
                            d2u_strch_cor(n_tds)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(r_du, conv, dud, d2u, recv_du_s, recv_du_e, &
    !$omp                 recv_dud_s, recv_dud_e, recv_d2u_s, recv_d2u_e, &
    !$omp                 du_sa, du_sc, du_strch, dud_sa, dud_sc, dud_strch, &
    !$omp                 d2u_sa, d2u_sc, d2u_strch, d2u_strch_cor)
    do b = 1, n_groups
      do i = 1, SZ
        call transeq_3fused_subs_line( &
          i, b, n, r_du, conv, dud, d2u, recv_du_s, recv_du_e, recv_dud_s, &
          recv_dud_e, recv_d2u_s, recv_d2u_e, n_tds, nu, &
          du_sa, du_sc, du_strch, dud_sa, dud_sc, dud_strch, &
          d2u_sa, d2u_sc, d2u_strch, d2u_strch_cor)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine transeq_3fused_subs

  subroutine transeq_3fused_subs_line( &
    i, b, n, r_du, conv, dud, d2u, recv_du_s, recv_du_e, recv_dud_s, &
    recv_dud_e, recv_d2u_s, recv_d2u_e, n_tds, nu, &
    du_sa, du_sc, du_strch, dud_sa, dud_sc, dud_strch, &
    d2u_sa, d2u_sc, d2u_strch, d2u_strch_cor &
    )
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds
    real(dp), intent(in) :: nu
    real(dp), intent(inout) :: r_du(SZ, n, *)
    real(dp), intent(in) :: conv(SZ, n, *), dud(SZ, n, *), d2u(SZ, n, *)
    real(dp), intent(in) :: recv_du_s(SZ, 1, *), recv_du_e(SZ, 1, *), &
                            recv_dud_s(SZ, 1, *), recv_dud_e(SZ, 1, *), &
                            recv_d2u_s(SZ, 1, *), recv_d2u_e(SZ, 1, *)
    real(dp), intent(in) :: du_sa(*), du_sc(*), du_strch(*), &
                            dud_sa(*), dud_sc(*), dud_strch(*), &
                            d2u_sa(*), d2u_sc(*), d2u_strch(*), &
                            d2u_strch_cor(*)

    integer :: j
    real(dp) :: ur, bl, recp
    real(dp) :: du_temp, dud_temp, d2u_temp
    real(dp) :: du_s, du_e, dud_s, dud_e, d2u_s, d2u_e

    ! Start: we have the 'bl', and assume 'ur'
    bl = du_sa(1)
    ur = du_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    du_s = recp*(r_du(i, 1, b) - bl*recv_du_s(i, 1, b))

    bl = dud_sa(1)
    ur = dud_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    dud_s = recp*(dud(i, 1, b) - bl*recv_dud_s(i, 1, b))

    bl = d2u_sa(1)
    ur = d2u_sa(1)
    recp = 1._dp/(1._dp - ur*bl)
    d2u_s = recp*(d2u(i, 1, b) - bl*recv_d2u_s(i, 1, b))

    ! End: we have the 'ur', and assume 'bl'
    bl = du_sc(n_tds)
    ur = du_sc(n_tds)
    recp = 1._dp/(1._dp - ur*bl)
    du_e = recp*(r_du(i, n_tds, b) - ur*recv_du_e(i, 1, b))

    bl = dud_sc(n_tds)
    ur = dud_sc(n_tds)
    recp = 1._dp/(1._dp - ur*bl)
    dud_e = recp*(dud(i, n_tds, b) - ur*recv_dud_e(i, 1, b))

    bl = d2u_sc(n_tds)
    ur = d2u_sc(n_tds)
    recp = 1._dp/(1._dp - ur*bl)
    d2u_e = recp*(d2u(i, n_tds, b) - ur*recv_d2u_e(i, 1, b))

    ! Final substitution
    r_du(i, 1, b) = -0.5_dp*(conv(i, 1, b)*du_s*du_strch(1) &
                             + dud_s*dud_strch(1)) &
                    + nu*(d2u_s*d2u_strch(1) &
                          + du_s*du_strch(1)*d2u_strch_cor(1))
    do j = 2, n_tds - 1
      du_temp = (r_du(i, j, b) - du_sa(j)*du_s - du_sc(j)*du_e)*du_strch(j)
      dud_temp = (dud(i, j, b) - dud_sa(j)*dud_s - dud_sc(j)*dud_e) &
                 *dud_strch(j)
      d2u_temp = (d2u(i, j, b) - d2u_sa(j)*d2u_s - d2u_sc(j)*d2u_e) &
                 *d2u_strch(j) + du_temp*d2u_strch_cor(j)
      r_du(i, j, b) = -0.5_dp*(conv(i, j, b)*du_temp + dud_temp) + nu*d2u_temp
    end do
    r_du(i, n_tds, b) = -0.5_dp*(conv(i, n_tds, b)*du_e*du_strch(n_tds) &
                                 + dud_e*dud_strch(n_tds)) &
                        + nu*(d2u_e*d2u_strch(n_tds) &
                              + du_e*du_strch(n_tds)*d2u_strch_cor(n_tds))

  end subroutine transeq_3fused_subs_line

  ! --------------------------------------------------------------------------
  ! Pentadiagonal kernels
  ! --------------------------------------------------------------------------

  subroutine penta_rhs_and_solve_line(i, b, n, du, u, u_s, u_e, n_tds, &
                                      n_rhs, coeffs_s, coeffs_e, coeffs, &
                                      ffr, faf, fsa, fbw, beta_lhs, &
                                      beta_lhs_s)
    !$omp declare target
    !! Builds the RHS of one line and solves the non-periodic pentadiagonal
    !! system in place: forward L-solve (unit diagonal, faf = l1, fsa = l2),
    !! then the U-solve dividing by d_j = 1/ffr(j). Shared by both
    !! pentadiagonal kernels, whose first stages are identical.
    integer, intent(in) :: i, b, n, n_tds, n_rhs
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(in) :: u(SZ, n, *), u_s(SZ, 4, *), u_e(SZ, 4, *)
    real(dp), intent(in) :: coeffs_s(9, *), coeffs_e(9, *), coeffs(*)
    real(dp), intent(in) :: ffr(*), faf(*), fsa(*), fbw(*)
    real(dp), intent(in) :: beta_lhs, beta_lhs_s

    integer :: j, row
    real(dp) :: c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3

    ! Bulk stencil uses indices 2..8 of the 9-element array
    c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8)

    do j = 1, 4
      du(i, j, b) = stencil_s(coeffs_s, j, i, b, n, u, u_s)
    end do
    do j = 5, n_rhs - 4
      du(i, j, b) = c_m3*u(i, j - 3, b) + c_m2*u(i, j - 2, b) &
                    + c_m1*u(i, j - 1, b) + c_j*u(i, j, b) &
                    + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                    + c_p3*u(i, j + 3, b)
    end do
    do row = 1, 4
      j = n_rhs - 4 + row
      du(i, j, b) = stencil_e(coeffs_e, row, n_rhs, i, b, n, u, u_e)
    end do

    ! Forward substitution: y_j = r_j - l1_j*y_{j-1} - l2_j*y_{j-2}
    du(i, 2, b) = du(i, 2, b) - faf(2)*du(i, 1, b)
    do j = 3, n_rhs
      du(i, j, b) = du(i, j, b) &
                    - faf(j)*du(i, j - 1, b) &
                    - fsa(j)*du(i, j - 2, b)
    end do

    ! Backward substitution
    du(i, n_tds, b) = du(i, n_tds, b)*ffr(n_tds)
    du(i, n_tds - 1, b) = (du(i, n_tds - 1, b) &
                           - fbw(n_tds - 1)*du(i, n_tds, b))*ffr(n_tds - 1)
    do j = n_tds - 2, 2, -1
      du(i, j, b) = (du(i, j, b) &
                     - fbw(j)*du(i, j + 1, b) &
                     - beta_lhs*du(i, j + 2, b))*ffr(j)
    end do
    ! j = 1 uses beta_lhs_s (0 for BC_NEUMANN sym, 2*beta antisym, else beta)
    du(i, 1, b) = (du(i, 1, b) &
                   - fbw(1)*du(i, 2, b) &
                   - beta_lhs_s*du(i, 3, b))*ffr(1)

  end subroutine penta_rhs_and_solve_line

  subroutine der_penta_full( &
    du, u, u_s, u_e, n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, &
    ffr, faf, fsa, fbw, beta_lhs, beta_lhs_s, n, n_groups &
    )
    !! Full non-periodic pentadiagonal solve on a single subdomain.
    integer, intent(in) :: n_tds, n_rhs, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(in) :: u(SZ, n, n_groups), u_s(SZ, 4, n_groups), &
                            u_e(SZ, 4, n_groups)
    real(dp), intent(in) :: coeffs_s(9, 4), coeffs_e(9, 4), coeffs(9)
    real(dp), intent(in) :: ffr(n_rhs), faf(n_rhs), fsa(n_rhs), fbw(n_rhs)
    real(dp), intent(in) :: beta_lhs, beta_lhs_s

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, u, u_s, u_e, coeffs_s, coeffs_e, coeffs, &
    !$omp                 ffr, faf, fsa, fbw)
    do b = 1, n_groups
      do i = 1, SZ
        call penta_rhs_and_solve_line(i, b, n, du, u, u_s, u_e, n_tds, &
                                      n_rhs, coeffs_s, coeffs_e, coeffs, &
                                      ffr, faf, fsa, fbw, beta_lhs, &
                                      beta_lhs_s)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_penta_full

  subroutine der_penta_periodic( &
    du, u, u_s, u_e, n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, &
    ffr, faf, fsa, fbw, beta_lhs, beta_lhs_s, alp, bet, z, n, n_groups &
    )
    !! Cyclic pentadiagonal solve via Sherman-Morrison-Woodbury (rank 4):
    !! A_cyc = A_np + U*C with U = [e_1|e_2|e_{n-1}|e_n]. `z` holds
    !! A_np^{-1}*U (see omptgt_tdsops_t%penta_z_dev); each line forms
    !! M = I + C*Z, solves M*c = C*y and applies y -= Z*c.
    integer, intent(in) :: n_tds, n_rhs, n, n_groups
    real(dp), intent(inout) :: du(SZ, n, n_groups)
    real(dp), intent(in) :: u(SZ, n, n_groups), u_s(SZ, 4, n_groups), &
                            u_e(SZ, 4, n_groups)
    real(dp), intent(in) :: coeffs_s(9, 4), coeffs_e(9, 4), coeffs(9)
    real(dp), intent(in) :: ffr(n_rhs), faf(n_rhs), fsa(n_rhs), fbw(n_rhs)
    real(dp), intent(in) :: beta_lhs, beta_lhs_s, alp, bet
    real(dp), intent(in) :: z(n_tds, 4)

    integer :: i, b

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(du, u, u_s, u_e, coeffs_s, coeffs_e, coeffs, &
    !$omp                 ffr, faf, fsa, fbw, z)
    do b = 1, n_groups
      do i = 1, SZ
        call penta_rhs_and_solve_line(i, b, n, du, u, u_s, u_e, n_tds, &
                                      n_rhs, coeffs_s, coeffs_e, coeffs, &
                                      ffr, faf, fsa, fbw, beta_lhs, &
                                      beta_lhs_s)
        call penta_periodic_correct_line(i, b, n, du, n_tds, alp, bet, z)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine der_penta_periodic

  subroutine penta_periodic_correct_line(i, b, n, du, n_tds, alp, bet, z)
    !$omp declare target
    integer, intent(in) :: i, b, n, n_tds
    real(dp), intent(inout) :: du(SZ, n, *)
    real(dp), intent(in) :: alp, bet
    real(dp), intent(in) :: z(n_tds, 4)

    integer :: j
    real(dp) :: m11, m12, m13, m14, m21, m22, m23, m24, &
                m31, m32, m33, m34, m41, m42, m43, m44
    real(dp) :: tmp1, tmp2, tmp3, tmp4, c4_1, c4_2, c4_3, c4_4

    ! M = I + C*Z (4x4), held in scalars so it stays in registers
    m11 = 1._dp + bet*z(n_tds - 1, 1) + alp*z(n_tds, 1)
    m12 = bet*z(n_tds - 1, 2) + alp*z(n_tds, 2)
    m13 = bet*z(n_tds - 1, 3) + alp*z(n_tds, 3)
    m14 = bet*z(n_tds - 1, 4) + alp*z(n_tds, 4)
    m21 = bet*z(n_tds, 1)
    m22 = 1._dp + bet*z(n_tds, 2)
    m23 = bet*z(n_tds, 3)
    m24 = bet*z(n_tds, 4)
    m31 = bet*z(1, 1)
    m32 = bet*z(1, 2)
    m33 = 1._dp + bet*z(1, 3)
    m34 = bet*z(1, 4)
    m41 = alp*z(1, 1) + bet*z(2, 1)
    m42 = alp*z(1, 2) + bet*z(2, 2)
    m43 = alp*z(1, 3) + bet*z(2, 3)
    m44 = 1._dp + alp*z(1, 4) + bet*z(2, 4)

    ! Doolittle LU factorisation of M, no pivoting
    m21 = m21/m11
    m22 = m22 - m21*m12
    m23 = m23 - m21*m13
    m24 = m24 - m21*m14
    m31 = m31/m11
    m32 = m32 - m31*m12
    m33 = m33 - m31*m13
    m34 = m34 - m31*m14
    m32 = m32/m22
    m33 = m33 - m32*m23
    m34 = m34 - m32*m24
    m41 = m41/m11
    m42 = m42 - m41*m12
    m43 = m43 - m41*m13
    m44 = m44 - m41*m14
    m42 = m42/m22
    m43 = m43 - m42*m23
    m44 = m44 - m42*m24
    m43 = m43/m33
    m44 = m44 - m43*m34

    ! rhs4 = C*y, then solve M*c4 = rhs4
    tmp1 = bet*du(i, n_tds - 1, b) + alp*du(i, n_tds, b)
    tmp2 = bet*du(i, n_tds, b)
    tmp3 = bet*du(i, 1, b)
    tmp4 = alp*du(i, 1, b) + bet*du(i, 2, b)

    tmp2 = tmp2 - m21*tmp1
    tmp3 = tmp3 - m31*tmp1 - m32*tmp2
    tmp4 = tmp4 - m41*tmp1 - m42*tmp2 - m43*tmp3

    c4_4 = tmp4/m44
    c4_3 = (tmp3 - m34*c4_4)/m33
    c4_2 = (tmp2 - m23*c4_3 - m24*c4_4)/m22
    c4_1 = (tmp1 - m12*c4_2 - m13*c4_3 - m14*c4_4)/m11

    ! Apply the correction du -= Z*c4
    do j = 1, n_tds
      du(i, j, b) = du(i, j, b) &
                    - z(j, 1)*c4_1 - z(j, 2)*c4_2 &
                    - z(j, 3)*c4_3 - z(j, 4)*c4_4
    end do

  end subroutine penta_periodic_correct_line

end module m_omptgt_kernels_dist
