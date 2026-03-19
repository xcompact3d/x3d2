program test_thom

  use iso_fortran_env, only: stderr => error_unit
#ifdef CUDA
  use cudafor
#else
  use omp_lib
<<<<<<<< HEAD:tests/unit/test_thom.f90
#endif
  use m_common, only: dp, nbytes, pi, BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, BC_HALO
#ifdef CUDA
  use m_cuda_common, only: SZ
  use m_cuda_exec_thom, only: exec_thom_tds_compact
  use m_cuda_tdsops, only: tdsops_t => cuda_tdsops_t, tdsops_init => cuda_tdsops_init
#else
========
  use m_common, only: dp, pi, BC_PERIODIC, BC_DIRICHLET
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_exec_thom, only: exec_thom_tds_compact
#endif
  implicit none

  logical :: allpass = .true.

#ifdef CUDA
  type(dim3) :: blocks, threads
#endif
  integer :: i, j, k
  integer :: n_glob, n, n_groups
<<<<<<<< HEAD:tests/unit/test_thom.f90
  integer :: n_iters
  integer :: ndof
  integer :: ierr

  real(kind(0.0d0)) :: tstart, tend
  real(dp), allocatable, dimension(:, :, :) :: u, du
#ifdef CUDA
  real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
#endif
  real(dp) :: dx, dx_per, norm_du
========

  real(dp), dimension(:, :, :), allocatable :: u, du
  real(dp) :: dx, dx_per
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90

  type(tdsops_t) :: tdsops

  !! Verification test
#ifdef CUDA
  n_glob = 1024
  n_groups = 512*512/SZ
  n_iters = 100
#else
  n_glob = 1024
  n_groups = 64*64/SZ
<<<<<<<< HEAD:tests/unit/test_thom.f90
  n_iters = 1
#endif
========
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
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

<<<<<<<< HEAD:tests/unit/test_thom.f90
#ifdef CUDA
  ! move data to device
  u_dev = u
#endif

  ! preprocess the operator and coefficient arrays
========
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  tdsops = tdsops_init(n, dx_per, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

<<<<<<<< HEAD:tests/unit/test_thom.f90
#ifdef CUDA
  blocks = dim3(n_groups, 1, 1)
  threads = dim3(SZ, 1, 1)
#endif

#ifdef CUDA
  call cpu_time(tstart)
#else
  tstart = omp_get_wtime()
#endif
  do i = 1, n_iters
#ifdef CUDA
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
#else
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
#endif
  end do
#ifdef CUDA
  call cpu_time(tend)
#else
  tend = omp_get_wtime()
#endif
  print *, "Total time", tend - tstart

#ifdef CUDA
  ! move data to host
  du = du_dev
#endif

#ifdef CUDA
  ! 2 in fw pass, 2 in bw pass, 2 in final periodic pass: 6 in total
  call checkperf(tend - tstart, n_iters, ndof, 6._dp)
#else
  call checkperf(tend - tstart, n_iters, ndof, 3.0_dp)
#endif
========
  call exec_thom_tds_compact(du, u, tdsops, n_groups)
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  call checkerr(u, du, 1.0e-8_dp)

  !! Dirichlet case
  print *, "=== Testing Dirichlet case ==="

  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx)
      end do
    end do
  end do

<<<<<<<< HEAD:tests/unit/test_thom.f90
#ifdef CUDA
  ! move data to device
  u_dev = u
#endif

  ! preprocess the operator and coefficient arrays
========
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  tdsops = tdsops_init(n, dx, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

<<<<<<<< HEAD:tests/unit/test_thom.f90
#ifdef CUDA
  call cpu_time(tstart)
#else
  tstart = omp_get_wtime()
#endif
  do i = 1, n_iters
#ifdef CUDA
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
#else
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
#endif
  end do
#ifdef CUDA
  call cpu_time(tend)
#else
  tend = omp_get_wtime()
#endif
  print *, "Total time", tend - tstart

#ifdef CUDA
  ! move data to host
  du = du_dev
#endif

#ifdef CUDA
  ! 2 in fw pass, 2 in bw pass: 4 in total
  call checkperf(tend - tstart, n_iters, ndof, 4._dp)
#else
  call checkperf(tend - tstart, n_iters, ndof, 3.0_dp)
#endif
========
  call exec_thom_tds_compact(du, u, tdsops, n_groups)
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  call checkerr(u, du, 1.0e-8_dp)

  if (allpass) then
    print *, 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

contains

<<<<<<<< HEAD:tests/unit/test_thom.f90
  subroutine checkperf(trun, n_iters, ndof, consumed_bw)
    implicit none

    real(kind(0.0d0)), intent(in) :: trun
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw

    integer :: ierr
    integer :: memClockRt, memBusWidth
    real(dp) :: achievedBW, deviceBW

#ifdef CUDA
    ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
    ierr = cudaDeviceGetAttribute(memBusWidth, &
                                  cudaDevAttrGlobalMemoryBusWidth, 0)
#else
    memClockRt = 3200000
    memBusWidth = 64
#endif

    ! BW utilisation and performance checks
    achievedBW = consumed_bw*n_iters*ndof*nbytes/trun
    deviceBW = 2.0_dp*memBusWidth/nbytes*memClockRt*(10**3)

    print *, "Check performance:"
    print'(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'
    print'(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
    print'(a, f5.2)', 'Effective BW util: %', achievedBW/deviceBW*100

  end subroutine checkperf

========
>>>>>>>> 0d54766 (split into performance and verification tests):tests/verification/test_omp_thom.f90
  subroutine checkerr(u, du, tol)

    real(dp), dimension(:, :, :), intent(in) :: u, du
    real(dp), intent(in) :: tol

    real(dp) :: norm_du

    norm_du = sum((u + du)**2)/n_glob/n_groups/SZ
    norm_du = sqrt(norm_du)

    print *, "Check error:"
    print *, "min:", minval(u + du), "max: ", maxval(u + du)
    print *, "error norm", norm_du

    if (norm_du > tol) then
      print *, "Check second derivatives... FAILED"
      allpass = .false.
    else
      print *, "Check second derivatives... PASSED"
    end if

  end subroutine checkerr

end program test_thom
