program test_thom

#ifdef CUDA
  use cudafor
#endif
  use m_common, only: dp, pi, BC_PERIODIC, BC_DIRICHLET
#ifdef CUDA
  use m_cuda_common, only: SZ
  use m_cuda_exec_thom, only: exec_thom_tds_compact
  use m_cuda_tdsops, only: tdsops_t => cuda_tdsops_t, tdsops_init => cuda_tdsops_init
#else
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_exec_thom, only: exec_thom_tds_compact
#endif
  use m_test_utils, only: checkerr
  implicit none

  real(dp), parameter :: residual_tol = 1.0e-8_dp
  logical :: allpass = .true.

#ifdef CUDA
  type(dim3) :: blocks, threads
#endif
  integer :: i, j, k
  integer :: n_glob, n, n_groups

  real(dp), allocatable, dimension(:, :, :) :: u, du
#ifdef CUDA
  real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
#endif
  real(dp) :: dx, dx_per

  type(tdsops_t) :: tdsops

#ifdef CUDA
  n_glob = 1024
  n_groups = 128*128/SZ
#else
  n_glob = 1024
  n_groups = 64*64/SZ
#endif
  n = n_glob

  allocate (u(SZ, n, n_groups), du(SZ, n, n_groups))
#ifdef CUDA
  allocate (u_dev(SZ, n, n_groups), du_dev(SZ, n, n_groups))
#endif

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  !! Periodic case
  print *, "=== Testing periodic case ==="

  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx_per)
      end do
    end do
  end do

#ifdef CUDA
  ! move data to device
  u_dev = u
#endif

  ! preprocess the operator and coefficient arrays
  tdsops = tdsops_init(n, dx_per, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

#ifdef CUDA
  blocks = dim3(n_groups, 1, 1)
  threads = dim3(SZ, 1, 1)
#endif
#ifdef CUDA
  call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
#else
  call exec_thom_tds_compact(du, u, tdsops, n_groups)
#endif

#ifdef CUDA
  ! move data to host
  du = du_dev
#endif

  call checkerr(u, du, residual_tol, 'thom_periodic', allpass)

  !! Dirichlet case
  print *, "=== Testing Dirichlet case ==="

  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx)
      end do
    end do
  end do

#ifdef CUDA
  ! move data to device
  u_dev = u
#endif

  ! preprocess the operator and coefficient arrays
  tdsops = tdsops_init(n, dx, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

#ifdef CUDA
  call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
#else
  call exec_thom_tds_compact(du, u, tdsops, n_groups)
#endif

#ifdef CUDA
  ! move data to host
  du = du_dev
#endif

  call checkerr(u, du, residual_tol, 'thom_dirichlet', allpass)

  if (allpass) then
      print *, 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

end program test_thom
