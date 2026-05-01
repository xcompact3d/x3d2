program test_omp_penta
  !! Verification test for the compact10_penta pentadiagonal first-derivative
  !! on the CPU/OMP backend.
  !!
  !! Runs four grid-refinement convergence studies:
  !!   1. BC_DIRICHLET: f = sin^3(pi*x),    require rate >= 4.
  !!   2. BC_NEUMANN sym=.true.:  f = cos(10*pi*x), require rate >= 9.
  !!   3. BC_NEUMANN sym=.false.: f = sin(10*pi*x), require rate >= 9.
  !!   4. BC_PERIODIC: f = sin(2*pi*x)+0.3*cos(4*pi*x), require rate >= 9.
  !!
  !! BC_DIRICHLET uses compact one-sided closures (4th-order, same alpha/beta as
  !! interior) at boundary rows; halos unused.  BC_NEUMANN uses mirror-ghost
  !! extension; halos filled with exact symmetric/antisymmetric extension.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_PERIODIC, BC_DIRICHLET, BC_NEUMANN
  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_penta_compact, exec_dist_penta_periodic
  use m_tdsops, only: tdsops_t, tdsops_init

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc, ierr

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  call run_dirichlet_test()
  call run_neumann_sym_true()
  call run_neumann_sym_false()
  call run_periodic_test()

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if
  call MPI_Finalize(ierr)

contains

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_DIRICHLET: f = sin(pi*x), f' = pi*cos(pi*x).
  ! Uses compact one-sided closures (4th-order) at rows 1-2 and N-1..N; halos unused.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_dirichlet_test()
    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    integer, parameter :: min_rate_tol = 4
    integer, parameter :: n_skip = 2   ! boundary rows with 4th-order closures
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev, l2_int, l2_int_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp; l2_int_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_DIRICHLET: f = sin^3(pi*x) on (0,1)  [global L2]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob + 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_interior(u, n, n_block, dx, 0, 'sin3')
      u_s(:, :, :) = 0._dp
      u_e(:, :, :) = 0._dp
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm(du, n, n_block, dx, 0, nproc, '3pi_sin2cos')
      call report_rate(l2_err, l2_prev, n_glob, isize, real(min_rate_tol, dp), &
                       'BC_DIRICHLET')
      l2_prev = l2_err
      l2_int = l2_norm_interior(du, n, n_block, dx, 0, nproc, '3pi_sin2cos', n_skip)
      l2_int_prev = l2_int   ! stored for interior table below
      deallocate (u, du, u_s, u_e)
    end do

    ! ---- Interior-only table (skip n_skip rows at each end, diagnostic only) ----
    if (nrank == 0) then
      print '(a)', ''
      print '(a,i0,a)', 'BC_DIRICHLET: interior only (skip ', n_skip, &
                         ' rows each end)  [diagnostic, no pass/fail]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if
    l2_int_prev = 0._dp
    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob + 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_interior(u, n, n_block, dx, 0, 'sin3')
      u_s(:, :, :) = 0._dp
      u_e(:, :, :) = 0._dp
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_int = l2_norm_interior(du, n, n_block, dx, 0, nproc, '3pi_sin2cos', n_skip)
      if (nrank == 0) then
        if (isize == 1) then
          print '(i6, es16.4, a10)', n_glob, l2_int, '   ---'
        else
          print '(i6, es16.4, f10.2)', n_glob, l2_int, &
            log(l2_int_prev/l2_int)/log(2.0_dp)
        end if
      end if
      l2_int_prev = l2_int
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_dirichlet_test

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.true.: f = cos(10*pi*x), f' = -10*pi*sin(10*pi*x).
  ! Ghost: even extension about x=0 and x=1 (f(-kh)=f(kh)).
  ! cos(10*pi*x) has natural even symmetry at both walls so extension is exact.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_true()
    integer, parameter :: n_sizes = 4
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256]
    real(dp), parameter :: min_rate_tol = 9.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.true.: f = cos(10*pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'cos10')
      ! Start ghost: even extension about x=0 (f(-kh)=f(kh), x_1=0).
      u_s(:, 4, :) = u(:, 2, :)     ! f(-h)  = f(h)  = u_2
      u_s(:, 3, :) = u(:, 3, :)     ! f(-2h) = f(2h) = u_3
      u_s(:, 2, :) = u(:, 4, :)     ! f(-3h) = f(3h) = u_4
      u_s(:, 1, :) = u(:, 5, :)     ! (coeff=0, fill anyway)
      ! End ghost: even extension about x=1 (f(1+kh)=f(1-kh), x_N=1).
      u_e(:, 1, :) = u(:, n - 1, :) ! f(1+h)  = f(1-h)  = u_{N-1}
      u_e(:, 2, :) = u(:, n - 2, :) ! f(1+2h) = f(1-2h) = u_{N-2}
      u_e(:, 3, :) = u(:, n - 3, :) ! f(1+3h) = f(1-3h) = u_{N-3}
      u_e(:, 4, :) = u(:, n - 4, :) ! (coeff=0, fill anyway)
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                           sym=.true.)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, 'neg10pi_sin10')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, &
                       'BC_NEUMANN sym=T')
      l2_prev = l2_err
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_neumann_sym_true

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.false.: f = sin(10*pi*x), f' = 10*pi*cos(10*pi*x).
  ! Ghost: odd extension about x=0 and x=1 (f(-kh)=-f(kh)).
  ! sin(10*pi*x) has natural odd symmetry at both walls so extension is exact.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_false()
    integer, parameter :: n_sizes = 4
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256]
    real(dp), parameter :: min_rate_tol = 9.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.false.: f = sin(10*pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'sin10')
      ! Start ghost: odd extension about x=0 (f(-kh)=-f(kh), x_1=0).
      u_s(:, 4, :) = -u(:, 2, :)     ! f(-h)  = -f(h)  = -u_2
      u_s(:, 3, :) = -u(:, 3, :)     ! f(-2h) = -f(2h) = -u_3
      u_s(:, 2, :) = -u(:, 4, :)     ! f(-3h) = -f(3h) = -u_4
      u_s(:, 1, :) = -u(:, 5, :)     ! (coeff=0, fill anyway)
      ! End ghost: odd extension about x=1 (f(1+kh)=-f(1-kh), x_N=1).
      u_e(:, 1, :) = -u(:, n - 1, :) ! f(1+h)  = -f(1-h)  = -u_{N-1}
      u_e(:, 2, :) = -u(:, n - 2, :) ! f(1+2h) = -f(1-2h) = -u_{N-2}
      u_e(:, 3, :) = -u(:, n - 3, :) ! f(1+3h) = -f(1-3h) = -u_{N-3}
      u_e(:, 4, :) = -u(:, n - 4, :) ! (coeff=0, fill anyway)
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                           sym=.false.)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, '10pi_cos10')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, &
                       'BC_NEUMANN sym=F')
      l2_prev = l2_err
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_neumann_sym_false

  ! ─────────────────────────────────────────────────────────────────────────────
  ! Field-fill helpers
  ! ─────────────────────────────────────────────────────────────────────────────

  subroutine fill_interior(u, n, n_block, dx, nrank, func)
    real(dp), intent(out) :: u(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank
    real(dp), intent(in) :: dx
    character(*), intent(in) :: func
    integer :: i, j, k
    real(dp) :: x
    do k = 1, n_block
      do j = 1, n
        x = real(nrank*n + j, dp)*dx
        do i = 1, SZ
          select case (func)
          case ('sin');       u(i, j, k) = sin(pi*x)
          case ('cos');       u(i, j, k) = cos(pi*x)
          case ('sin3');      u(i, j, k) = sin(pi*x)**3
          case ('sin_per');   u(i, j, k) = sin(2._dp*pi*x) + 0.3_dp*cos(4._dp*pi*x)
          end select
        end do
      end do
    end do
  end subroutine fill_interior

  subroutine fill_wall(u, n, n_block, dx, nrank, func)
    real(dp), intent(out) :: u(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank
    real(dp), intent(in) :: dx
    character(*), intent(in) :: func
    integer :: i, j, k
    real(dp) :: x
    do k = 1, n_block
      do j = 1, n
        x = real(nrank*n + j - 1, dp)*dx
        do i = 1, SZ
          select case (func)
          case ('sin');       u(i, j, k) = sin(pi*x)
          case ('cos');       u(i, j, k) = cos(pi*x)
          case ('cos10');      u(i, j, k) = cos(10._dp*pi*x)
          case ('sin10');      u(i, j, k) = sin(10._dp*pi*x)
          end select
        end do
      end do
    end do
  end subroutine fill_wall

  ! ─────────────────────────────────────────────────────────────────────────────
  ! L2-error helpers
  ! ─────────────────────────────────────────────────────────────────────────────

  real(dp) function l2_norm(du, n, n_block, dx, nrank, nproc, deriv)
    real(dp), intent(in) :: du(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank, nproc
    real(dp), intent(in) :: dx
    character(*), intent(in) :: deriv
    integer :: j, k
    real(dp) :: x, err, sum_sq, global_sum
    sum_sq = 0._dp
    do k = 1, n_block
      do j = 1, n
        x = real(nrank*n + j, dp)*dx
        err = du(1, j, k) - exact_deriv(x, deriv)
        sum_sq = sum_sq + err*err
      end do
    end do
    call MPI_Allreduce(sum_sq, global_sum, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    l2_norm = sqrt(global_sum/real(n*nproc*n_block, dp))
  end function l2_norm

  real(dp) function l2_norm_interior(du, n, n_block, dx, nrank, nproc, deriv, n_skip)
    !! L2 error on interior grid, skipping n_skip points at each end.
    real(dp), intent(in) :: du(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank, nproc, n_skip
    real(dp), intent(in) :: dx
    character(*), intent(in) :: deriv
    integer :: j, k, n_pts
    real(dp) :: x, err, sum_sq, global_sum
    sum_sq = 0._dp
    n_pts = 0
    do k = 1, n_block
      do j = n_skip + 1, n - n_skip
        x = real(nrank*n + j, dp)*dx
        err = du(1, j, k) - exact_deriv(x, deriv)
        sum_sq = sum_sq + err*err
        n_pts = n_pts + 1
      end do
    end do
    call MPI_Allreduce(sum_sq, global_sum, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    l2_norm_interior = sqrt(global_sum/real(n_pts*nproc, dp))
  end function l2_norm_interior

  real(dp) function l2_norm_wall(du, n, n_block, dx, nrank, nproc, deriv)
    real(dp), intent(in) :: du(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank, nproc
    real(dp), intent(in) :: dx
    character(*), intent(in) :: deriv
    integer :: j, k
    real(dp) :: x, err, sum_sq, global_sum
    sum_sq = 0._dp
    do k = 1, n_block
      do j = 1, n
        x = real(nrank*n + j - 1, dp)*dx
        err = du(1, j, k) - exact_deriv(x, deriv)
        sum_sq = sum_sq + err*err
      end do
    end do
    call MPI_Allreduce(sum_sq, global_sum, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    l2_norm_wall = sqrt(global_sum/real(n*nproc*n_block, dp))
  end function l2_norm_wall

  real(dp) pure function exact_deriv(x, which)
    real(dp), intent(in) :: x
    character(*), intent(in) :: which
    select case (which)
    case ('pi_cos');          exact_deriv = pi*cos(pi*x)
    case ('neg_pi_sin');      exact_deriv = -pi*sin(pi*x)
    case ('3pi_sin2cos');     exact_deriv = 3._dp*pi*sin(pi*x)**2*cos(pi*x)
    case ('neg10pi_sin10');   exact_deriv = -10._dp*pi*sin(10._dp*pi*x)
    case ('10pi_cos10');      exact_deriv = 10._dp*pi*cos(10._dp*pi*x)
    case ('per_deriv');       exact_deriv = 2._dp*pi*cos(2._dp*pi*x) &
                                            - 1.2_dp*pi*sin(4._dp*pi*x)
    end select
  end function exact_deriv

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_PERIODIC: f = sin(2πx)+0.3cos(4πx). Grid: x_j=j/n, j=1..n, h=1/n.
  ! SMW rank-4 correction restores full 10th-order accuracy.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_periodic_test()
    integer, parameter :: n_sizes = 4
    integer, parameter :: n_glob_arr(n_sizes) = [8, 16, 32, 64]
    real(dp), parameter :: min_rate_tol = 9.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev, rate
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_PERIODIC: f = sin(2πx)+0.3cos(4πx) on [0,1)'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_interior(u, n, n_block, dx, 0, 'sin_per')
      ! Periodic halos: u_s wraps the end of u back to the start.
      u_s(:, 4, :) = u(:, n, :)       ! x_{j-1} for j=1: x_{n}  -> u_n
      u_s(:, 3, :) = u(:, n - 1, :)   ! x_{j-2} for j=1: x_{n-1}
      u_s(:, 2, :) = u(:, n - 2, :)   ! x_{j-3} for j=1: x_{n-2}
      u_s(:, 1, :) = u(:, n - 3, :)   ! x_{j-4} for j=1: x_{n-3}
      u_e(:, 1, :) = u(:, 1, :)       ! x_{j+1} for j=n: x_1 (wrap)
      u_e(:, 2, :) = u(:, 2, :)       ! x_{j+2} for j=n: x_2
      u_e(:, 3, :) = u(:, 3, :)       ! x_{j+3} for j=n: x_3
      u_e(:, 4, :) = u(:, 4, :)       ! x_{j+4} for j=n: x_4
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
      call exec_dist_penta_periodic(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm(du, n, n_block, dx, 0, nproc, 'per_deriv')
      if (nrank == 0) then
        if (isize == 1) then
          print '(i6, es16.4, a10)', n_glob, l2_err, '   ---'
        else
          rate = log(l2_prev/l2_err)/log(2.0_dp)
          print '(i6, es16.4, f10.2)', n_glob, l2_err, rate
          if (rate < min_rate_tol) then
            allpass = .false.
            write (stderr, '(a,f5.2,a,f4.1,a)') &
              'BC_PERIODIC convergence check... failed (rate = ', rate, &
              ', expected >= ', min_rate_tol, ')'
          end if
        end if
      end if
      l2_prev = l2_err
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_periodic_test

  ! ─────────────────────────────────────────────────────────────────────────────
  ! Print one convergence-table row and check the rate.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine report_rate(l2_err, l2_prev, n_glob, isize, min_rate, label)
    real(dp), intent(in) :: l2_err, l2_prev, min_rate
    integer, intent(in) :: n_glob, isize
    character(*), intent(in) :: label
    real(dp) :: rate
    if (nrank /= 0) return
    if (isize == 1) then
      print '(i6, es16.4, a10)', n_glob, l2_err, '   ---'
    else if (l2_err < 1e-12_dp) then
      print '(i6, es16.4, a10)', n_glob, l2_err, '  <eps'
    else
      rate = log(l2_prev/l2_err)/log(2.0_dp)
      print '(i6, es16.4, f10.2)', n_glob, l2_err, rate
      if (rate < min_rate) then
        allpass = .false.
        write (stderr, '(a, f5.2, a, f4.1, a)') &
          trim(label)//' convergence check... failed (rate = ', rate, &
          ', expected >= ', min_rate, ')'
      end if
    end if
  end subroutine report_rate

end program test_omp_penta
