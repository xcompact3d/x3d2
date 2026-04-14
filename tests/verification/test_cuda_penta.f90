program test_cuda_penta
  !! Verification test for the compact10_penta pentadiagonal first-derivative.
  !!
  !! Strategy: build the same non-periodic pentadiagonal LHS and RHS on the
  !! host, solve via a host-side pentadiagonal Thomas, then compare with the
  !! CUDA kernel output.  This verifies the kernel's RHS assembly, LU
  !! factorisation, and Thomas forward/backward passes are all correct.
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_penta_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, du_kernel
  real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
  real(dp), device, allocatable, dimension(:, :, :) :: &
    u_recv_s_dev, u_recv_e_dev

  type(cuda_tdsops_t) :: tdsops

  integer :: n, n_block, n_halo, n_glob
  integer :: nrank, nproc
  integer :: ierr, ndevs, devnum

  type(dim3) :: blocks, threads
  real(dp) :: dx

  call initialise_mpi()
  call select_device()
  call setup_geometry()
  call allocate_fields()
  call initialise_input()
  call setup_backend()
  call run_kernel()
  call check_kernel_vs_host()
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

  subroutine setup_geometry()
    n_glob = 256
    n = n_glob/nproc
    n_block = 1
    n_halo = 4
    dx = 2*pi/n_glob
  end subroutine setup_geometry

  subroutine allocate_fields()
    allocate (u(SZ, n, n_block), du_kernel(SZ, n, n_block))
    allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))
    allocate (u_recv_s_dev(SZ, n_halo, n_block))
    allocate (u_recv_e_dev(SZ, n_halo, n_block))
  end subroutine allocate_fields

  subroutine initialise_input()
    integer :: i, j, k
    do k = 1, n_block
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = sin((j - 1 + nrank*n)*dx)
        end do
      end do
    end do
    u_dev = u
  end subroutine initialise_input

  subroutine setup_backend()
    tdsops = cuda_tdsops_init(n, dx, operation='first-deriv', &
                              scheme='compact10_penta', &
                              bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
    blocks = dim3(n_block, 1, 1)
    threads = dim3(SZ, 1, 1)
  end subroutine setup_backend

  subroutine run_kernel()
    u_recv_s_dev(:, :, :) = u_dev(:, n - n_halo + 1:n, :)
    u_recv_e_dev(:, :, :) = u_dev(:, 1:n_halo, :)
    call exec_dist_penta_compact(du_dev, u_dev, u_recv_s_dev, u_recv_e_dev, &
                                 tdsops, blocks, threads)
    du_kernel = du_dev
  end subroutine run_kernel

  subroutine check_kernel_vs_host()
    !! Build the same RHS and solve the same non-periodic pentadiagonal system
    !! on the host, then compare with the CUDA kernel result.
    real(dp), parameter :: kernel_tol = 1.0e-10_dp

    real(dp) :: alp, bet
    real(dp) :: afi, bfi, cfi
    real(dp), allocatable :: u_ext(:)
    real(dp), allocatable :: rhs(:), du_ref(:)
    ! LU factors for host Thomas
    real(dp), allocatable :: d_inv(:), l1(:), l2(:), u1(:)
    integer :: j
    real(dp) :: d_j, norm_diff, norm_ref, maxdiff

    alp = 0.5_dp
    bet = 0.05_dp
    afi = (17._dp/12._dp)/(2._dp*dx)
    bfi = (101._dp/150._dp)/(4._dp*dx)
    cfi = (1._dp/100._dp)/(6._dp*dx)

    ! Extended u with periodic halos
    allocate (u_ext(1 - n_halo:n + n_halo))
    do j = 1, n
      u_ext(j) = u(1, j, 1)
    end do
    do j = 1, n_halo
      u_ext(1 - j) = u_ext(n - j + 1)
      u_ext(n + j) = u_ext(j)
    end do

    ! Build RHS
    allocate (rhs(n), du_ref(n))
    do j = 1, n
      rhs(j) = -cfi*u_ext(j - 3) - bfi*u_ext(j - 2) - afi*u_ext(j - 1) &
               + afi*u_ext(j + 1) + bfi*u_ext(j + 2) + cfi*u_ext(j + 3)
    end do

    ! --- Host pentadiagonal Thomas (non-periodic) ---
    ! LHS: bet*x_{i-2} + alp*x_{i-1} + 1*x_i + alp*x_{i+1} + bet*x_{i+2}
    allocate (d_inv(n), l1(n), l2(n), u1(n))

    ! LU factorisation (same as preprocess_penta_dist)
    l2(1) = 0._dp; l1(1) = 0._dp
    d_inv(1) = 1._dp; u1(1) = alp

    l2(2) = 0._dp; l1(2) = alp
    d_j = 1._dp - alp*alp
    d_inv(2) = 1._dp/d_j
    u1(2) = alp*(1._dp - bet)

    do j = 3, n
      l2(j) = bet*d_inv(j - 2)
      l1(j) = (alp - l2(j)*u1(j - 2))*d_inv(j - 1)
      d_j = (1._dp - l2(j)*bet) - l1(j)*u1(j - 1)
      d_inv(j) = 1._dp/d_j
      u1(j) = alp - l1(j)*bet
    end do

    ! Forward substitution
    du_ref(:) = rhs(:)
    du_ref(1) = du_ref(1)*d_inv(1)
    du_ref(2) = (du_ref(2) - l1(2)*du_ref(1))*d_inv(2)
    do j = 3, n
      du_ref(j) = (du_ref(j) - l1(j)*du_ref(j - 1) &
                   - l2(j)*du_ref(j - 2))*d_inv(j)
    end do

    ! Backward substitution
    du_ref(n - 1) = du_ref(n - 1) - u1(n - 1)*du_ref(n)
    do j = n - 2, 1, -1
      du_ref(j) = du_ref(j) - u1(j)*du_ref(j + 1) - bet*du_ref(j + 2)
    end do

    ! Compare
    norm_diff = 0._dp; norm_ref = 0._dp; maxdiff = 0._dp
    do j = 1, n
      norm_diff = norm_diff + (du_kernel(1, j, 1) - du_ref(j))**2
      norm_ref = norm_ref + du_ref(j)**2
      maxdiff = max(maxdiff, abs(du_kernel(1, j, 1) - du_ref(j)))
    end do
    norm_diff = sqrt(norm_diff/n)
    norm_ref = sqrt(norm_ref/n)

    if (nrank == 0) then
      print '(a, es12.4)', ' kernel vs host L2 norm: ', norm_diff
      print '(a, es12.4)', ' host reference L2 norm: ', norm_ref
      print '(a, es12.4)', ' max pointwise diff:     ', maxdiff
      print '(a, es12.4)', ' relative L2 error:      ', norm_diff/norm_ref

      if (norm_diff > kernel_tol) then
        allpass = .false.
        write (stderr, '(a)') 'Check penta kernel vs host... failed'
      else
        write (stderr, '(a)') 'Check penta kernel vs host... passed'
      end if
    end if

    deallocate (u_ext, rhs, du_ref, d_inv, l1, l2, u1)
  end subroutine check_kernel_vs_host

  subroutine finalise()
    if (allpass) then
      if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if
    call MPI_Finalize(ierr)
  end subroutine finalise

end program test_cuda_penta
