program test_cuda_penta
  !! Verification test for the compact10_penta pentadiagonal first-derivative.
  !!
  !! Strategy: sweep over a sequence of grid sizes (N=32,64,128,256,512) with
  !! BC_DIRICHLET on both ends.  At each N compare the GPU result against the
  !! analytical derivative of f(x) = sin(pi*x) on [0,1] and compute the L2
  !! error.  The scheme should converge at ~10th order in the interior; we
  !! require the measured rate between consecutive grid sizes to be >= 4.
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_DIRICHLET
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_penta_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc, ierr, ndevs, devnum

  call initialise_mpi()
  call select_device()
  call run_convergence_test()
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

  subroutine run_convergence_test()
    !! Grid-refinement convergence study.
    !!
    !! f(x) = sin(pi*x) on [0,1], homogeneous Dirichlet at both ends.
    !! Analytical derivative: f'(x) = pi*cos(pi*x).
    !! Grid points: x_j = j/(N+1), j = 1..N  (interior, NOT including boundary).
    !! u_recv_s/e are zero (homogeneous Dirichlet halos).

    integer, parameter :: n_sizes = 5
    integer, parameter :: n_glob_arr(n_sizes) = [32, 64, 128, 256, 512]
    real(dp), parameter :: rate_tol = 4.0_dp
    real(dp), parameter :: min_rate_tol = 4.0_dp

    integer :: isize, n_glob, n, n_block, n_halo
    real(dp) :: dx, l2_err
    real(dp) :: l2_prev
    real(dp) :: rate

    real(dp), allocatable, dimension(:, :, :) :: u, du
    real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
    real(dp), device, allocatable, dimension(:, :, :) :: u_s_dev, u_e_dev

    type(cuda_tdsops_t) :: tdsops
    type(dim3) :: blocks, threads

    n_block = 1
    n_halo = 4

    if (nrank == 0) then
      print '(a)', ''
      print '(a)', 'compact10_penta BC_DIRICHLET convergence: f = sin(pi*x) on [0,1]'
      print '(a6, a16, a10)', 'N', 'L2 error', 'Rate'
    end if

    l2_prev = 0._dp

    do isize = 1, n_sizes
      n_glob = n_glob_arr(isize)
      ! Interior points only; boundaries are x=0 and x=1 (Dirichlet, not stored)
      n = n_glob / nproc
      ! Grid spacing: (N+1) intervals over [0,1] → h = 1/(N+1)
      dx = 1._dp / real(n_glob + 1, dp)

      allocate (u(SZ, n, n_block), du(SZ, n, n_block))
      allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
      allocate (u_s_dev(SZ, n_halo, n_block))
      allocate (u_e_dev(SZ, n_halo, n_block))

      call fill_field(u, n, n_block, dx, nrank, nproc)
      u_dev = u
      ! Homogeneous Dirichlet: halos are zero (f(0)=f(1)=0 for sin(pi*x))
      u_s_dev = 0._dp
      u_e_dev = 0._dp

      tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                                scheme='compact10_penta', &
                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

      blocks = dim3(n_block, 1, 1)
      threads = dim3(SZ, 1, 1)

      call exec_dist_penta_compact(du_dev, u_dev, u_s_dev, u_e_dev, &
                                   tdsops, blocks, threads)
      du = du_dev

      l2_err = compute_l2_error(du, n, n_block, dx, nrank, nproc)

      if (nrank == 0) then
        if (isize == 1) then
          print '(i6, es16.4, a10)', n_glob, l2_err, '   ---'
        else
          rate = log(l2_prev / l2_err) / log(2.0_dp)
          print '(i6, es16.4, f10.2)', n_glob, l2_err, rate
          if (rate < min_rate_tol) then
            allpass = .false.
            write (stderr, '(a, f5.2, a)') &
              'Convergence check... failed (rate = ', rate, &
              ', expected >= 4)'
          end if
        end if
      end if

      l2_prev = l2_err

      deallocate (u, du, u_dev, du_dev, u_s_dev, u_e_dev)
    end do

    if (nrank == 0 .and. allpass) then
      write (stderr, '(a)') 'Convergence check... passed'
    end if
  end subroutine run_convergence_test

  subroutine fill_field(u, n, n_block, dx, nrank, nproc)
    !! Fill u with sin(pi*x) at interior grid points.
    !! x_j = (global_j) * dx, global_j = nrank*n + local_j.
    real(dp), intent(out) :: u(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank, nproc
    real(dp), intent(in) :: dx
    integer :: i, j, k
    real(dp) :: x
    do k = 1, n_block
      do j = 1, n
        x = real(nrank * n + j, dp) * dx
        do i = 1, SZ
          u(i, j, k) = sin(pi * x)
        end do
      end do
    end do
  end subroutine fill_field

  real(dp) function compute_l2_error(du, n, n_block, dx, nrank, nproc)
    !! L2 norm of (du_kernel - pi*cos(pi*x)) over all interior points.
    real(dp), intent(in) :: du(SZ, n, n_block)
    integer, intent(in) :: n, n_block, nrank, nproc
    real(dp), intent(in) :: dx
    integer :: i, j, k
    real(dp) :: x, err, sum_sq, global_sum
    sum_sq = 0._dp
    do k = 1, n_block
      do j = 1, n
        x = real(nrank * n + j, dp) * dx
        err = du(1, j, k) - pi * cos(pi * x)
        sum_sq = sum_sq + err * err
      end do
    end do
    call MPI_Allreduce(sum_sq, global_sum, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    compute_l2_error = sqrt(global_sum / real(n * nproc * n_block, dp))
  end function compute_l2_error

  subroutine finalise()
    if (allpass) then
      if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if
    call MPI_Finalize(ierr)
  end subroutine finalise

end program test_cuda_penta
