program test_thom

  use m_common, only: dp, pi
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_exec_thom, only: exec_thom_tds_compact
  
  implicit none
  
  integer :: n_glob, n, n_groups
  integer :: n_iters
  integer :: ndof

  real(dp), dimension(:, :, :), allocatable :: u, du
  real(dp) :: dx, dx_per

  integer :: i, j, k

  type(tdsops_t) :: tdsops

  real(dp) :: tstart, tend

  logical :: allpass

  allpass = .true.
  
  n_glob = 1024
  n = n_glob
  n_groups = 64 * 64 / SZ
  n_iters = 1
  ndof = n_glob * n_groups * SZ

  allocate(u(SZ, n, n_groups), du(SZ, n, n_groups))

  dx_per = 2 * pi / n_glob
  dx = 2 * pi / (n_glob - 1)

  !! Periodic case
  
  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1) * dx_per)
      end do
    end do
  end do

  tdsops = tdsops_init(n, dx_per, &
                       operation = "second-deriv", scheme = "compact6", &
                       bc_start = "periodic", bc_end = "periodic")

  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du, u, tdsops)
  end do
  call cpu_time(tend)
  print *, "Total time", tend - tstart

  call checkperf(tend - tstart, n_iters, ndof, 5.0_dp)
  call checkerr(u, du, 1.0e-8_dp)

  !! Dirichlet case
  
  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1) * dx)
      end do
    end do
  end do

  tdsops = tdsops_init(n, dx_per, &
                       operation = "second-deriv", scheme = "compact6", &
                       bc_start = "dirichlet", bc_end = "dirichlet")

  call cpu_time(tstart)
  do i = 1, n_iters
    call exec_thom_tds_compact(du, u, tdsops)
  end do
  call cpu_time(tend)
  print *, "Total time", tend - tstart

  call checkperf(tend - tstart, n_iters, ndof, 5.0_dp)
  call checkerr(u, du, 1.0e-8_dp)

  
  error stop "STILL UNDER IMPLEMENTATION"

contains

  subroutine checkperf(trun, n_iters, ndof, consumed_bw)

    real(dp), intent(in) :: trun
    integer, intent(in) :: n_iters
    integer, intent(in) :: ndof
    real(dp), intent(in) :: consumed_bw

    real(dp) :: nbytes
    real(dp) :: achievedBW

    real(dp) :: memClockRt, memBusWidth, deviceBW

    if (dp == kind(0.0d0)) then
      nbytes = 8_dp
    else
      nbytes = 4_dp
    end if
    achievedBW = consumed_bw * n_iters * ndof * nbytes / trun

    print *, "Achieved BW: ", achievedBW / (2**30), " GiB/s"

    memClockRt = 3200
    memBusWidth = 64
    deviceBW = 2 * memBusWidth / nbytes * memClockRt * (10**6)

    print *, "Available BW: ", deviceBW / 2**30, " GiB/s / NUMA zone on ARCHER2"
    print *, "Utilised BW: ", 100 * (achievedBW / deviceBW), " %"
    
  end subroutine checkperf

  subroutine checkerr(u, du, tol)

    real(dp), dimension(:, :, :), intent(in) :: u, du
    real(dp), intent(in) :: tol
    
    real(dp) :: norm_du

    norm_du = sum((u + du)**2) / n_glob / n_groups / SZ
    norm_du = sqrt(norm_du)

    print *, "error norm", norm_du

    if (norm_du > tol) then
      print *, "Check Dirichlet second derivatives... FAILED"
      allpass = .false.
    else
      print *,  "Check Dirichlet second derivatives... PASSED"
    end if
    
  end subroutine checkerr
  
end program test_thom
