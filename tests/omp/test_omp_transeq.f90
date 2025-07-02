program test_omp_transeq
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, pi, DIR_X, DIR_Y, DIR_Z, VERT
  use m_omp_common, only: SZ
  use m_omp_sendrecv, only: sendrecv_fields
  use m_omp_backend, only: omp_backend_t, transeq_x_omp, base_backend_t
  use m_tdsops, only: dirps_t, tdsops_t
  use m_solver, only: allocate_tdsops
  use m_mesh, only: mesh_t

  implicit none

  logical :: allpass = .true.
  class(field_t), pointer :: u, v, w
  class(field_t), pointer :: du, dv, dw
  real(dp), dimension(:, :, :), allocatable :: r_u
  class(mesh_t), allocatable :: mesh
  integer, dimension(3) :: dims_global, nproc_dir
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)

  integer :: n, n_groups, i, j, k
  integer :: nrank, nproc
  integer :: ierr

  real(dp) :: dx_per, nu, norm_du, tol = 1d-8, tstart, tend

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator

  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
  type(dirps_t) :: xdirps, ydirps, zdirps

  ! Initialise variables and arrays
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! Global number of cells in each direction
  dims_global = [96, 96, 96]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, nproc]

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)

  xdirps%dir = DIR_X; ydirps%dir = DIR_Y; zdirps%dir = DIR_Z

  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator
  print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  print *, 'OpenMP backend instantiated'

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  n = mesh%get_n(DIR_X, VERT)
  n_groups = mesh%get_n_groups(DIR_X)

  nu = 1._dp

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)

  du => allocator%get_block(DIR_X, VERT)
  dv => allocator%get_block(DIR_X, VERT)
  dw => allocator%get_block(DIR_X, VERT)

  dx_per = mesh%geo%d(DIR_X)

  do k = 1, n_groups
    do j = 1, n
      do i = 1, SZ
        u%data(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
        v%data(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
      end do
    end do
  end do
  w%data(:, :, :) = 0.d0

  call allocate_tdsops(xdirps, omp_backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')

  call cpu_time(tstart)
  call transeq_x_omp(omp_backend, du, dv, dw, u, v, w, nu, xdirps)
  call cpu_time(tend)

  if (nrank == 0) print *, 'Total time', tend - tstart

  allocate (r_u(SZ, n, n_groups))

  ! check error
  ! dv = -1/2*(u*dv/dx + d(u*v)/dx) + nu*d2v/dx2
  ! u is sin, v is cos;
  ! dv = -1/2*(u*(-u) + v*v + u*(-u)) + nu*(-v)
  !    = u*u - 1/2*v*v - nu*v
  r_u = dv%data - (u%data*u%data - 0.5_dp*v%data*v%data - nu*v%data)
  norm_du = norm2(r_u)
  norm_du = norm_du*norm_du/dims_global(DIR_X)/n_groups/SZ
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

