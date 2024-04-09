program test_omp_transeq
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, pi, globs_t, set_pprev_pnext, DIR_X, DIR_Y, DIR_Z
  use m_omp_common, only: SZ
  use m_omp_sendrecv, only: sendrecv_fields
  use m_omp_backend, only: omp_backend_t, transeq_x_omp, base_backend_t
  use m_tdsops, only: dirps_t, tdsops_t
  use m_solver, only: allocate_tdsops

  implicit none

  logical :: allpass = .true.
  class(field_t), pointer :: u, v, w
  class(field_t), pointer :: du, dv, dw
  real(dp), dimension(:, :, :), allocatable :: r_u

  integer :: n, n_block, i, j, k
  integer :: n_glob
  integer :: nrank, nproc
  integer :: ierr

  real(dp) :: dx, dx_per, nu, norm_du, tol = 1d-8, tstart, tend

  type(globs_t) :: globs
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator

  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
  type(dirps_t) :: xdirps, ydirps, zdirps

  ! Initialise variables and arrays
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  globs%nx = 96
  globs%ny = 96
  globs%nz = 96

  globs%nx_loc = globs%nx/nproc
  globs%ny_loc = globs%ny/nproc
  globs%nz_loc = globs%nz/nproc

  globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
  globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
  globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

  xdirps%nproc = nproc
  ydirps%nproc = nproc
  zdirps%nproc = nproc

  call set_pprev_pnext( &
    xdirps%pprev, xdirps%pnext, &
    ydirps%pprev, ydirps%pnext, &
    zdirps%pprev, zdirps%pnext, &
    xdirps%nproc, ydirps%nproc, zdirps%nproc, nrank &
    )

  xdirps%dir = DIR_X
  ydirps%dir = DIR_Y
  zdirps%dir = DIR_Z

  omp_allocator = allocator_t(globs%nx_loc, globs%ny_loc, globs%nz_loc, SZ)
  allocator => omp_allocator
  print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(allocator)
  backend => omp_backend
  print *, 'OpenMP backend instantiated'

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  n_glob = globs%nx
  n = n_glob/nproc
  n_block = globs%n_groups_x

  nu = 1._dp
  omp_backend%nu = nu

  u => allocator%get_block(DIR_X)
  v => allocator%get_block(DIR_X)
  w => allocator%get_block(DIR_X)

  du => allocator%get_block(DIR_X)
  dv => allocator%get_block(DIR_X)
  dw => allocator%get_block(DIR_X)

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)
  globs%dx = dx

  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u%data(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
        v%data(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
      end do
    end do
  end do
  w%data(:, :, :) = 0.d0

  call allocate_tdsops(xdirps, globs%nx_loc, dx_per, omp_backend)

  call cpu_time(tstart)
  call transeq_x_omp(omp_backend, du, dv, dw, u, v, w, xdirps)
  call cpu_time(tend)

  if (nrank == 0) print *, 'Total time', tend - tstart

  allocate (r_u(SZ, n, n_block))

  ! check error
  ! dv = -1/2*(u*dv/dx + d(u*v)/dx) + nu*d2v/dx2
  ! u is sin, v is cos;
  ! dv = -1/2*(u*(-u) + v*v + u*(-u)) + nu*(-v)
  !    = u*u - 1/2*v*v - nu*v
  r_u = dv%data - (u%data*u%data - 0.5_dp*v%data*v%data - nu*v%data)
  norm_du = norm2(r_u)
  norm_du = norm_du*norm_du/n_glob/n_block/SZ
  call MPI_Allreduce(MPI_IN_PLACE, norm_du, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  norm_du = sqrt(norm_du)

  if (nrank == 0) print *, 'error norm', norm_du
  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check second derivatives... failed'
    else
      write (stderr, '(a)') 'Check second derivatives... passed'
    end if
  end if

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

end program test_omp_transeq

