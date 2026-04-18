program test_cuda_penta
  !! Verification test for the compact10_penta pentadiagonal first-derivative.
  !!
  !! Runs three grid-refinement convergence studies:
  !!   1. BC_DIRICHLET: f = sin(pi*x), x in (0,1), require rate >= 4.
  !!   2. BC_NEUMANN sym=.true.:  f = cos(pi*x), require rate >= 4.
  !!   3. BC_NEUMANN sym=.false.: f = sin(pi*x), require rate >= 4.
  !!
  !! Grid conventions:
  !!   BC_DIRICHLET: N interior points, x_j = j*h,   h = 1/(N+1), j = 1..N.
  !!                 Walls at x=0 and x=1 are NOT grid points; halos = 0.
  !!   BC_NEUMANN:   N total points,   x_j = (j-1)*h, h = 1/(N-1), j = 1..N.
  !!                 Walls ARE grid points x_1=0 and x_N=1.  Halos = 0 (stencils
  !!                 at wall rows are designed not to read halo positions).
  !!                 Placing the wall exactly at x_1 makes the even/odd ghost-cell
  !!                 extension exact for cos and sin, yielding interior-class
  !!                 convergence rates everywhere.
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_DIRICHLET, BC_NEUMANN
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_penta_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc, ierr, ndevs, devnum

  call initialise_mpi()
  call select_device()
  call run_dirichlet_test()
  call run_neumann_sym_true()
  call run_neumann_sym_false()
  call finalise()

contains

  subroutine initialise_mpi()
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'
  end subroutine initialise_mpi

  subroutine select_device()
    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(devnum)
  end subroutine select_device

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_DIRICHLET: f = sin(pi*x), f' = pi*cos(pi*x), zero halos.
  ! Grid: N interior points, x_j = j*h, h = 1/(N+1).
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_dirichlet_test()
    integer, parameter :: n_sizes = 3
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du
    real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
    real(dp), device, allocatable, dimension(:, :, :) :: u_s_dev, u_e_dev
    type(cuda_tdsops_t) :: tdsops
    type(dim3) :: blocks, threads

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_DIRICHLET: f = sin(pi*x) on (0,1)'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob + 1, dp)   ! interior grid: x_j = j*h
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
      allocate (u_s_dev(SZ, n_halo, n_block), u_e_dev(SZ, n_halo, n_block))
      call fill_interior(u, n, n_block, dx, 0, 'sin')
      u_dev = u; u_s_dev = 0._dp; u_e_dev = 0._dp
      tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                                scheme='compact10_penta', &
                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
      blocks = dim3(n_block, 1, 1); threads = dim3(SZ, 1, 1)
      call exec_dist_penta_compact(du_dev, u_dev, u_s_dev, u_e_dev, &
                                   tdsops, blocks, threads)
      du = du_dev
      l2_err = l2_norm(du, n, n_block, dx, 0, nproc, 'pi_cos')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, 'BC_DIRICHLET')
      l2_prev = l2_err
      deallocate (u, du, u_dev, du_dev, u_s_dev, u_e_dev)
    end do
  end subroutine run_dirichlet_test

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.true.: f = cos(pi*x), f' = -pi*sin(pi*x).
  ! Grid: N total points, x_j = (j-1)*h, h = 1/(N-1).
  ! Walls at x_1=0 and x_N=1; f'(0)=f'(1)=0 exactly enforced.
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_true()
    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du
    real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
    real(dp), device, allocatable, dimension(:, :, :) :: u_s_dev, u_e_dev
    type(cuda_tdsops_t) :: tdsops
    type(dim3) :: blocks, threads

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.true.: f = cos(pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)   ! wall grid: x_j = (j-1)*h, x_1=0
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
      allocate (u_s_dev(SZ, n_halo, n_block), u_e_dev(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'cos')
      ! Wall stencils do not reference halo positions; set to zero.
      u_dev = u; u_s_dev = 0._dp; u_e_dev = 0._dp
      tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                                scheme='compact10_penta', &
                                bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                                sym=.true.)
      blocks = dim3(n_block, 1, 1); threads = dim3(SZ, 1, 1)
      call exec_dist_penta_compact(du_dev, u_dev, u_s_dev, u_e_dev, &
                                   tdsops, blocks, threads)
      du = du_dev
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, 'neg_pi_sin')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, &
                       'BC_NEUMANN sym=T')
      l2_prev = l2_err
      deallocate (u, du, u_dev, du_dev, u_s_dev, u_e_dev)
    end do
  end subroutine run_neumann_sym_true

  ! ─────────────────────────────────────────────────────────────────────────────
  ! BC_NEUMANN sym=.false.: f = sin(pi*x), f' = pi*cos(pi*x).
  ! Grid: N total points, x_j = (j-1)*h, h = 1/(N-1).
  ! Walls at x_1=0 and x_N=1; f(0)=f(1)=0 (odd extension exact for sin).
  ! ─────────────────────────────────────────────────────────────────────────────
  subroutine run_neumann_sym_false()
    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    real(dp), parameter :: min_rate_tol = 4.0_dp
    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err, l2_prev
    real(dp), allocatable, dimension(:, :, :) :: u, du
    real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
    real(dp), device, allocatable, dimension(:, :, :) :: u_s_dev, u_e_dev
    type(cuda_tdsops_t) :: tdsops
    type(dim3) :: blocks, threads

    n_block = 1; n_halo = 4; l2_prev = 0._dp

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'BC_NEUMANN sym=.false.: f = sin(pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      n = n_glob
      dx = 1._dp/real(n_glob - 1, dp)   ! wall grid: x_j = (j-1)*h, x_1=0
      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
      allocate (u_s_dev(SZ, n_halo, n_block), u_e_dev(SZ, n_halo, n_block))
      call fill_wall(u, n, n_block, dx, 0, 'sin')
      ! Wall stencils do not reference halo positions; set to zero.
      u_dev = u; u_s_dev = 0._dp; u_e_dev = 0._dp
      tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                                scheme='compact10_penta', &
                                bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
                                sym=.false.)
      blocks = dim3(n_block, 1, 1); threads = dim3(SZ, 1, 1)
      call exec_dist_penta_compact(du_dev, u_dev, u_s_dev, u_e_dev, &
                                   tdsops, blocks, threads)
      du = du_dev
      l2_err = l2_norm_wall(du, n, n_block, dx, 0, nproc, 'pi_cos')
      call report_rate(l2_err, l2_prev, n_glob, isize, min_rate_tol, &
                       'BC_NEUMANN sym=F')
      l2_prev = l2_err
      deallocate (u, du, u_dev, du_dev, u_s_dev, u_e_dev)
    end do
  end subroutine run_neumann_sym_false

  ! ─────────────────────────────────────────────────────────────────────────────
  ! Field-fill helpers
  ! ─────────────────────────────────────────────────────────────────────────────

  subroutine fill_interior(u, n, n_block, dx, nrank, func)
    !! Fill using interior grid x_j = (nrank*n + j) * dx.
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
    !! Fill using wall grid x_j = (nrank*n + j - 1) * dx  (x_1=0 for rank 0).
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
    !! L2 error on interior grid x_j = (nrank*n + j) * dx.
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
    !! L2 error on wall grid x_j = (nrank*n + j - 1) * dx.
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
    !! Analytical first derivative.
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

  subroutine finalise()
    if (allpass) then
      if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if
    call MPI_Finalize(ierr)
  end subroutine finalise

end program test_cuda_penta
