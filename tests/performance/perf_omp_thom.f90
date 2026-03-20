program perf_omp_thom

  use omp_lib
  use m_common, only: dp, nbytes, pi, BC_PERIODIC, BC_DIRICHLET
  use m_omp_common, only: SZ
  use m_tdsops, only: tdsops_t, tdsops_init
  use m_exec_thom, only: exec_thom_tds_compact
  use m_test_utils, only: report_perf

  implicit none

  integer :: i, j, k
  integer :: n_glob, n, n_groups
  integer :: n_iters, n_warmup
  integer :: ndof

  real(dp) :: tstart, tend
  real(dp), dimension(:, :, :), allocatable :: u, du
  real(dp) :: dx_per, dx

  type(tdsops_t) :: tdsops

  ! Performance configuration
  n_glob = 512
  n_groups = 128*128/SZ
  n_iters = 500
  n_warmup = 10
  n = n_glob
  ndof = n_glob*n_groups*SZ

  allocate (u(SZ, n, n_groups), du(SZ, n, n_groups))

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  ! Periodic case benchmark
  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx_per)
      end do
    end do
  end do

  tdsops = tdsops_init(n, dx_per, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  ! Warmup
  do i = 1, n_warmup
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
  end do

  tstart = omp_get_wtime()
  do i = 1, n_iters
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
  end do
  tend = omp_get_wtime()

  call report_perf('omp_thom_periodic', tend - tstart, n_iters, ndof, 3.0_dp)

  ! Dirichlet case benchmark
  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1)*dx)
      end do
    end do
  end do

  tdsops = tdsops_init(n, dx, &
                       operation="second-deriv", scheme="compact6", &
                       bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

  ! Warmup
  do i = 1, n_warmup
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
  end do

  tstart = omp_get_wtime()
  do i = 1, n_iters
    call exec_thom_tds_compact(du, u, tdsops, n_groups)
  end do
  tend = omp_get_wtime()

  call report_perf('omp_thom_dirichlet', tend - tstart, n_iters, ndof, 3.0_dp)

end program perf_omp_thom
