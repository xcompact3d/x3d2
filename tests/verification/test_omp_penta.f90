program test_omp_penta
  !! Verification test for the compact10_penta pentadiagonal first-derivative
  !! on the CPU/OMP backend.
  !!
  !! Runs three grid-refinement convergence studies:
  !!   1. BC_DIRICHLET: f = sin(pi*x), x in (0,1), require rate >= 4.
  !!   2. BC_NEUMANN sym=.true.:  f = cos(pi*x), require rate >= 4.
  !!   3. BC_NEUMANN sym=.false.: f = sin(pi*x), require rate >= 4.
  !!
  !! Grid conventions match test_cuda_penta.f90 exactly.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_DIRICHLET, BC_NEUMANN
  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_penta_compact
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

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if
  call MPI_Finalize(ierr)

contains

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_DIRICHLET: f = sin(pi*x), f' = pi*cos(pi*x), zero halos.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_dirichlet_test()
    integer, parameter :: n_sizes = 3
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_DIRICHLET: f = sin(pi*x) on (0,1)'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob + 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_interior(u, n, n_block, dx, 0, 'sin')
      u_s = 0._dp; u_e = 0._dp
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm(du, n, n_block, dx, 0, nproc, 'pi_cos')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, 'BC_DIRICHLET')
      l2_prev = l2_err
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_dirichlet_test

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.true.: f = cos(pi*x), f' = -pi*sin(pi*x).
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_true()
    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.true.: f = cos(pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'cos')
      u_s = 0._dp; u_e = 0._dp
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                           sym=.true.)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, 'neg_pi_sin')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, &
                       'BC_NEUMANN sym=T')
      l2_prev = l2_err
      deallocate (u, du, u_s, u_e)
    end do
  end subroutine run_neumann_sym_true

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.false.: f = sin(pi*x), f' = pi*cos(pi*x).
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_false()
    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du, u_s, u_e
    type(tdsops_t) :: tdsops

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.false.: f = sin(pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_s(SZ, n_halo, n_block), u_e(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'sin')
      u_s = 0._dp; u_e = 0._dp
      tdsops = tdsops_init(n, dx, operation='first-deriv', &
                           scheme='compact10_penta', &
                           bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                           sym=.false.)
      call exec_dist_penta_compact(du, u, u_s, u_e, tdsops, n_block)
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, 'pi_cos')
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
          case ('sin'); u(i, j, k) = sin(pi*x)
          case ('cos'); u(i, j, k) = cos(pi*x)
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
          case ('sin'); u(i, j, k) = sin(pi*x)
          case ('cos'); u(i, j, k) = cos(pi*x)
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
    case ('pi_cos');     exact_deriv = pi*cos(pi*x)
    case ('neg_pi_sin'); exact_deriv = -pi*sin(pi*x)
    end select
  end function exact_deriv

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
    else
      rate = log(l2_prev/l2_err)/log(2.0_dp)
      print '(i6, es16.4, f10.2)', n_glob, l2_err, rate
      if (rate < min_rate) then
        allpass = .false.
        write (stderr, '(a, f5.2, a)') &
          trim(label)//' convergence check... failed (rate = ', rate, &
          ', expected >= 4)'
      end if
    end if
  end subroutine report_rate

end program test_omp_penta
