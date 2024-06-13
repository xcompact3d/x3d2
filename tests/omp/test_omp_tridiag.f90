program test_omp_tridiag
  use iso_fortran_env, only: stderr => error_unit
  use mpi
  use omp_lib

  use m_common, only: dp, pi, DIR_X, VERT
  use m_omp_common, only: SZ
  use m_omp_sendrecv, only: sendrecv_fields
  use m_omp_exec_dist, only: exec_dist_tds_compact
  use m_allocator, only: allocator_t
  use m_field, only: field_t

  use m_tdsops, only: tdsops_t, tdsops_init
  use m_mesh, only: mesh_t

  implicit none

  logical :: allpass = .true.

  real(dp), allocatable, dimension(:, :, :) :: u, du
  class(field_t), pointer :: u_field, du_field
  real(dp), allocatable, dimension(:, :, :) :: u_recv_s, u_recv_e, &
                                               u_send_s, u_send_e

  real(dp), allocatable, dimension(:, :, :) :: send_s, send_e, &
                                               recv_s, recv_e

  real(dp), allocatable, dimension(:) :: sin_0_2pi_per, cos_0_2pi_per, &
                                         sin_0_2pi, cos_0_2pi, &
                                         sin_stag, cos_stag

  type(tdsops_t) :: tdsops

  character(len=20) :: bc_start, bc_end

  integer :: n, n_groups, j, n_halo, n_iters, n_loc
  integer :: n_glob
  integer :: nrank, nproc, pprev, pnext
  integer :: ierr, memClockRt, memBusWidth
  real(dp), dimension(3) :: L_global
  integer, dimension(3) :: dims_global, nproc_dir
  class(mesh_t), allocatable :: mesh
  class(allocator_t), allocatable :: allocator

  real(dp) :: dx, dx_per, norm_du, tol = 1d-8, tstart, tend
  real(dp) :: achievedBW, deviceBW, achievedBWmax, achievedBWmin

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  ! Global number of cells in each direction
  dims_global = [1024, 64, 64]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [nproc, 1, 1]

  mesh = mesh_t(dims_global, nproc_dir, L_global)
  allocator = allocator_t(mesh, SZ)
  
  pnext = mesh%par%pnext(DIR_X)
  pprev = mesh%par%pprev(DIR_X)

  n_glob = dims_global(1)
  n = mesh%get_n(DIR_X, VERT)
  n_groups = mesh%get_n_groups(DIR_X)
  n_iters = 1

  u_field => allocator%get_block(DIR_X, VERT)
  du_field => allocator%get_block(DIR_X, VERT)
  u = u_field%data
  du = du_field%data

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  allocate (sin_0_2pi_per(n), cos_0_2pi_per(n))
  allocate (sin_0_2pi(n), cos_0_2pi(n))
  allocate (sin_stag(n), cos_stag(n))
  do j = 1, n
    sin_0_2pi_per(j) = sin(((j - 1) + nrank*n)*dx_per)
    cos_0_2pi_per(j) = cos(((j - 1) + nrank*n)*dx_per)
    sin_0_2pi(j) = sin(((j - 1) + nrank*n)*dx)
    cos_0_2pi(j) = cos(((j - 1) + nrank*n)*dx)
    sin_stag(j) = sin(((j - 1) + nrank*n)*dx + dx/2._dp)
    cos_stag(j) = cos(((j - 1) + nrank*n)*dx + dx/2._dp)
  end do

  n_halo = 4

  ! arrays for exchanging data between ranks
  allocate (u_send_s(SZ, n_halo, n_groups))
  allocate (u_send_e(SZ, n_halo, n_groups))
  allocate (u_recv_s(SZ, n_halo, n_groups))
  allocate (u_recv_e(SZ, n_halo, n_groups))

  allocate (send_s(SZ, 1, n_groups), send_e(SZ, 1, n_groups))
  allocate (recv_s(SZ, 1, n_groups), recv_e(SZ, 1, n_groups))

  ! =========================================================================
  ! second derivative with periodic BC
  tdsops = tdsops_init(n, dx_per, operation='second-deriv', scheme='compact6')

  call set_u(u, sin_0_2pi_per, n, n_groups)

  tstart = omp_get_wtime()

  call run_kernel(n_iters, n_groups, u, du, tdsops, n, &
                  u_recv_s, u_recv_e, u_send_s, u_send_e, &
                  recv_s, recv_e, send_s, send_e, &
                  nproc, pprev, pnext &
                  )

  tend = omp_get_wtime()
  if (nrank == 0) print *, 'Total time', tend - tstart

  call check_error_norm(du, sin_0_2pi_per, n, n_glob, n_groups, 1, norm_du)
  if (nrank == 0) print *, 'error norm second-deriv periodic', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check 2nd derivatives, periodic BCs... failed'
    else
      write (stderr, '(a)') 'Check 2nd derivatives, periodic BCs... passed'
    end if
  end if

  ! =========================================================================
  ! first derivative with periodic BC
  tdsops = tdsops_init(n, dx_per, operation='first-deriv', scheme='compact6')

  call set_u(u, sin_0_2pi_per, n, n_groups)

  call run_kernel(n_iters, n_groups, u, du, tdsops, n, &
                  u_recv_s, u_recv_e, u_send_s, u_send_e, &
                  recv_s, recv_e, send_s, send_e, &
                  nproc, pprev, pnext &
                  )

  call check_error_norm(du, cos_0_2pi_per, n, n_glob, n_groups, -1, norm_du)
  if (nrank == 0) print *, 'error norm first-deriv periodic', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check 1st derivatives, periodic BCs... failed'
    else
      write (stderr, '(a)') 'Check 1st derivatives, periodic BCs... passed'
    end if
  end if

  ! =========================================================================
  ! first derivative with dirichlet and neumann
  if (nrank == 0) then
    bc_start = 'dirichlet'!'neumann'!'dirichlet'
  else
    bc_start = 'null'
  end if
  if (nrank == nproc - 1) then
    bc_end = 'neumann'!'dirichlet'!'neumann'
  else
    bc_end = 'null'
  end if

  tdsops = tdsops_init(n, dx, operation='first-deriv', scheme='compact6', &
                       bc_start=trim(bc_start), bc_end=trim(bc_end), &
                       sym=.false.)

  call set_u(u, sin_0_2pi, n, n_groups)

  call run_kernel(n_iters, n_groups, u, du, tdsops, n, &
                  u_recv_s, u_recv_e, u_send_s, u_send_e, &
                  recv_s, recv_e, send_s, send_e, &
                  nproc, pprev, pnext &
                  )

  call check_error_norm(du, cos_0_2pi, n, n_glob, n_groups, -1, norm_du)
  if (nrank == 0) print *, 'error norm first deriv dir-neu', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check 1st derivatives, dir-neu... failed'
    else
      write (stderr, '(a)') 'Check 1st derivatives, dir-neu... passed'
    end if
  end if

  ! =========================================================================
  ! stag interpolate with neumann sym
  n_loc = n
  if (nrank == nproc - 1) n_loc = n - 1
  tdsops = tdsops_init(n_loc, dx, operation='interpolate', scheme='classic', &
                       bc_start=trim(bc_start), bc_end=trim(bc_end), &
                       from_to='v2p')

  call set_u(u, cos_0_2pi, n, n_groups)

  call run_kernel(n_iters, n_groups, u, du, tdsops, n_loc, &
                  u_recv_s, u_recv_e, u_send_s, u_send_e, &
                  recv_s, recv_e, send_s, send_e, &
                  nproc, pprev, pnext &
                  )

  call check_error_norm(du, cos_stag, n_loc, n_glob, n_groups, -1, norm_du)
  if (nrank == 0) print *, 'error norm interpolate', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check interpolation... failed'
    else
      write (stderr, '(a)') 'Check interpolation... passed'
    end if
  end if

  ! =========================================================================
  ! second derivative and hyperviscousity on with dirichlet and neumann
  ! c_nu = 0.22 and nu0_nu = 63 results in alpha = 0.40869111947709036
  tdsops = tdsops_init(n, dx, operation='second-deriv', &
                       scheme='compact6-hyperviscous', &
                       bc_start=trim(bc_start), bc_end=trim(bc_end), &
                       sym=.false., c_nu=0.22_dp, nu0_nu=63._dp)

  call set_u(u, sin_0_2pi, n, n_groups)

  call run_kernel(n_iters, n_groups, u, du, tdsops, n, &
                  u_recv_s, u_recv_e, u_send_s, u_send_e, &
                  recv_s, recv_e, send_s, send_e, &
                  nproc, pprev, pnext &
                  )

  call check_error_norm(du, sin_0_2pi, n, n_glob, n_groups, 1, norm_du)
  if (nrank == 0) print *, 'error norm hyperviscous', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check 2nd ders, hyperviscous, dir-neu... failed'
    else
      write (stderr, '(a)') 'Check 2nd ders, hyperviscous, dir-neu... passed'
    end if
  end if

  ! =========================================================================
  ! BW utilisation and performance checks
  ! 3 in the first phase, 2 in the second phase, so 5 in total
  achievedBW = 5._dp*n_iters*n*n_groups*SZ*dp/(tend - tstart)
  call MPI_Allreduce(achievedBW, achievedBWmax, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(achievedBW, achievedBWmin, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MIN, MPI_COMM_WORLD, ierr)
  if (nrank == 0) then
    print'(a, f8.3, a)', 'Achieved BW min: ', achievedBWmin/2**30, ' GiB/s'
    print'(a, f8.3, a)', 'Achieved BW max: ', achievedBWmax/2**30, ' GiB/s'
  end if

  memClockRt = 3200
  memBusWidth = 64
  deviceBW = 2*memBusWidth/8._dp*memClockRt*1000000

  if (nrank == 0) then
    print'(a, f8.3, a)', 'Available BW:   ', deviceBW/2**30, &
      ' GiB/s (per NUMA zone on ARCHER2)'
    print'(a, f5.2)', 'Effective BW util min: %', achievedBWmin/deviceBW*100
    print'(a, f5.2)', 'Effective BW util max: %', achievedBWmax/deviceBW*100
  end if

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine run_kernel(n_iters, n_groups, u, du, tdsops, n, &
                        u_recv_s, u_recv_e, u_send_s, u_send_e, &
                        recv_s, recv_e, send_s, send_e, &
                        nproc, pprev, pnext)
    implicit none

    integer, intent(in) :: n_iters, n_groups
    real(dp), intent(in), dimension(:, :, :) :: u
    real(dp), intent(out), dimension(:, :, :) :: du
    type(tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n
    real(dp), intent(inout), dimension(:, :, :) :: u_recv_s, u_recv_e, &
                                                   u_send_s, u_send_e
    real(dp), intent(inout), dimension(:, :, :) :: recv_s, recv_e, &
                                                   send_s, send_e
    integer, intent(in) :: nproc, pprev, pnext

    integer :: iters, i, j, k

    do iters = 1, n_iters
      ! first copy halo data into buffers
      !$omp parallel do
      do k = 1, n_groups
        do j = 1, 4
          !$omp simd
          do i = 1, SZ
            u_send_s(i, j, k) = u(i, j, k)
            u_send_e(i, j, k) = u(i, n - n_halo + j, k)
          end do
          !$omp end simd
        end do
      end do
      !$omp end parallel do

      ! halo exchange
      call sendrecv_fields(u_recv_s, u_recv_e, u_send_s, u_send_e, &
                           SZ*n_halo*n_groups, nproc, pprev, pnext)

     call exec_dist_tds_compact(du, u, u_recv_s, u_recv_e, &
                                send_s, send_e, recv_s, recv_e, &
                                tdsops, nproc, pprev, pnext, n_groups)

    end do
  end subroutine run_kernel

  subroutine set_u(u, line, n, n_groups)
    implicit none

    real(dp), intent(out), dimension(:, :, :) :: u
    real(dp), intent(in), dimension(:) :: line
    integer, intent(in) :: n, n_groups

    integer :: i, j, k

    do k = 1, n_groups
      do j = 1, n
        do i = 1, SZ
          u(i, j, k) = line(j)
        end do
      end do
    end do

  end subroutine set_u

  subroutine check_error_norm(du, line, n, n_glob, n_groups, c, norm)
    implicit none

    real(dp), intent(inout), dimension(:, :, :) :: du
    real(dp), intent(in), dimension(:) :: line
    integer, intent(in) :: n, n_glob, n_groups, c
    real(dp), intent(out) :: norm

    integer :: i, j, k

    do k = 1, n_groups
      do j = 1, n
        do i = 1, SZ
          du(i, j, k) = du(i, j, k) + c*line(j)
        end do
      end do
    end do

    norm = norm2(du(:, 1:n, :))
    norm = norm*norm/n_glob/n_groups/SZ
    call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    norm = sqrt(norm)

  end subroutine check_error_norm

end program test_omp_tridiag

