program test_thom
  use iso_fortran_env, only: stderr => error_unit
  use cudafor

  use m_common, only: dp, pi, BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, BC_HALO
  use m_cuda_common, only: SZ
  use m_cuda_exec_thom, only: exec_thom_tds_compact
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, du
  real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev

  type(cuda_tdsops_t) :: tdsops

  integer :: n, n_block, i, j, k, n_iters, ndof
  integer :: n_glob
  integer :: ierr, ndevs, devnum, memClockRt, memBusWidth

  type(dim3) :: blocks, threads
  real(dp) :: dx, dx_per, norm_du, tol = 1d-8, tstart, tend
  real(dp) :: achievedBW, deviceBW

  n_glob = 512*2
  n = n_glob
  n_block = 512*512/SZ
  n_iters = 100
  ndof = n_glob*n_block*SZ

  allocate (u(SZ, n, n_block), du(SZ, n, n_block))
  allocate (u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx_per)
      end do
    end do
  end do

  ! move data to device
  u_dev = u

  ! preprocess the operator and coefficient arrays
  tdsops = cuda_tdsops_init(n, dx_per, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  blocks = dim3(n_block, 1, 1)
  threads = dim3(SZ, 1, 1)

  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  call cpu_time(tend)
  print *, 'Total time', tend - tstart

  ! 2 in fw pass, 2 in bw pass, 2 in final periodic pass: 6 in total
  call checkperf(tend - tstart, n_iters, ndof, 6._dp)

  ! check error
  du = du_dev
  norm_du = norm2(u + du)
  norm_du = norm_du*norm_du/n_glob/n_block/SZ
  norm_du = sqrt(norm_du)

  print *, 'error norm', norm_du

  if (norm_du > tol) then
    allpass = .false.
    write (stderr, '(a)') 'Check periodic second derivatives... failed'
  else
    write (stderr, '(a)') 'Check periodic second derivatives... passed'
  end if

  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx)
      end do
    end do
  end do

  ! move data to device
  u_dev = u

  ! preprocess the operator and coefficient arrays
  tdsops = cuda_tdsops_init(n, dx, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du_dev, u_dev, tdsops, blocks, threads)
  end do
  call cpu_time(tend)
  print *, 'Total time', tend - tstart

  ! 2 in fw pass, 2 in bw pass: 4 in total
  call checkperf(tend - tstart, n_iters, ndof, 4._dp)

  ! check error
  du = du_dev
  norm_du = norm2(u + du)
  norm_du = norm_du*norm_du/n_glob/n_block/SZ
  norm_du = sqrt(norm_du)

  print *, 'error norm', norm_du

  if (norm_du > tol) then
    allpass = .false.
    write (stderr, '(a)') 'Check dirichlet second derivatives... failed'
  else
    write (stderr, '(a)') 'Check dirichlet second derivatives... passed'
  end if
  if (allpass) then
    write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

contains

  subroutine checkperf(t_tot, n_iters, ndof, consumed_bw)
    implicit none

    real(dp), intent(in) :: t_tot, consumed_bw
    integer, intent(in) :: n_iters, ndof

    real(dp) :: achievedBW, devBW
    integer :: ierr, memClockRt, memBusWidth

    ! BW utilisation and performance checks
    achievedBW = consumed_bw*n_iters*ndof*dp/t_tot

    print'(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'

    ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
    ierr = cudaDeviceGetAttribute(memBusWidth, &
                                  cudaDevAttrGlobalMemoryBusWidth, 0)
    devBW = 2*memBusWidth/8._dp*memClockRt*1000

    print'(a, f8.3, a)', 'Device BW:   ', devBW/2**30, ' GiB/s'
    print'(a, f5.2)', 'Effective BW util: %', achievedBW/devBW*100

  end subroutine checkperf

end program test_thom

