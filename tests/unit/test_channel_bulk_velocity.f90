program test_channel_bulk_velocity
  !! Verifies the channel bulk velocity on a y stretched mesh.
  !!
  !! The profile u = y**2*(2 - y) on 0 <= y <= 2 is zero at both walls and has
  !! the exact bulk velocity 2/3. It is not symmetric, so a wrong reduction
  !! over the y decomposition changes the result. The y direction is split
  !! over all ranks, so every rank uses its own slice of vert_ds.
  !!
  !! Check 1 compares channel_bulk_velocity with a host reference that has its
  !! own MPI_Allreduce and does not use the backend. Check 2 compares it with
  !! the exact value, within the quadrature error of the trapezoid rule.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_INTEGER, MPI_LAND, &
                   MPI_LOGICAL, MPI_SUM, MPI_Allreduce, MPI_Comm_rank, &
                   MPI_Comm_size, MPI_Finalize, MPI_Init

  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_case_channel, only: channel_bulk_velocity
  use m_common, only: dp, MPI_X3D2_DP, DIR_X, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  class(field_t), pointer :: u

  ! ny divides evenly over 1 and 2 ranks, so both walls are present.
  integer, parameter :: ny = 64
  ! beta of the channel smoke case
  real(dp), parameter :: beta = 0.259065151_dp
  real(dp), parameter :: ub_exact = 2._dp/3._dp
  ! Check 1: the two sums differ only in the order of the additions. The
  ! round-off of a sum of n terms is at most n*epsilon, here n = 4*ny.
  real(dp), parameter :: tol_ref = 4._dp*real(ny, dp)*epsilon(1._dp)
  ! Check 2: trapezoid error of a profile that vanishes at the walls is
  ! (h**2/12)*|u'(L)|*y'(L)**2/L, with h = 1/(ny - 1), u'(L) = -4 and
  ! y'(L) at most L = 2. This is 8*h**2/12, and 2*h**2 bounds it.
  real(dp), parameter :: tol_quad = 2._dp/real(ny - 1, dp)**2
  character(len=10), parameter :: uni = 'uniform', tb = 'top-bottom'
  character(len=8), parameter :: periodic(2) = ['periodic', 'periodic']
  character(len=9), parameter :: walls(2) = ['dirichlet', 'dirichlet']

  real(dp), allocatable :: host_u(:, :, :)
  real(dp) :: ub, ub_ref, y, err_ref, err_quad
  integer :: dims(3), nproc_dir(3)
  integer :: nvert, j, ierr, nrank, nproc
  logical :: test_pass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  nproc_dir = [1, nproc, 1]
  mesh = mesh_t([16, ny, 16], nproc_dir, [4._dp, 2._dp, 2._dp], &
                periodic, walls, periodic, &
                stretching=[uni, tb, uni], beta=[0._dp, beta, 0._dp])

  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend

  ! Exactly one rank must own each vertex plane.
  dims = mesh%get_dims(VERT)
  nvert = dims(2)
  call MPI_Allreduce(MPI_IN_PLACE, nvert, 1, MPI_INTEGER, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)
  test_pass = nvert == ny

  allocate (host_u(dims(1), dims(2), dims(3)))
  do j = 1, dims(2)
    y = mesh%geo%vert_coords(j, 2)
    host_u(:, j, :) = y**2*(2._dp - y)
  end do

  u => allocator%get_block(DIR_X, VERT)
  call backend%set_field_data(u, host_u)

  ub = channel_bulk_velocity(backend, mesh, u)

  ! Host reference: a plane mean of a field that is constant in x and z is
  ! its value, so the weighted sum needs no plane sums.
  ub_ref = 0._dp
  do j = 1, dims(2)
    ub_ref = ub_ref + host_u(1, j, 1)/mesh%geo%vert_ds(j, 2)
  end do
  ub_ref = ub_ref*mesh%geo%d(2)/mesh%geo%L(2)
  call MPI_Allreduce(MPI_IN_PLACE, ub_ref, 1, MPI_X3D2_DP, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)

  err_ref = abs(ub - ub_ref)
  err_quad = abs(ub - ub_exact)
  test_pass = test_pass .and. ieee_is_finite(ub) .and. .not. ieee_is_nan(ub) &
              .and. ieee_is_finite(ub_ref) .and. .not. ieee_is_nan(ub_ref)
  test_pass = test_pass .and. err_ref <= tol_ref
  test_pass = test_pass .and. err_quad < tol_quad
  call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, MPI_LOGICAL, MPI_LAND, &
                     MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    print '(a,i0,a,i0)', 'ranks: ', nproc, ', vertex planes: ', nvert
    print '(a,es25.16e3)', 'ub:           ', ub
    print '(a,es25.16e3)', 'host ref:     ', ub_ref
    print '(a,es25.16e3,a,es25.16e3)', '|ub - ref|:   ', err_ref, &
      '  tol: ', tol_ref
    print '(a,es25.16e3,a,es25.16e3)', '|ub - 2/3|:   ', err_quad, &
      '  tol: ', tol_quad
    if (test_pass) then
      print *, 'PASS'
    else
      print *, 'FAIL'
    end if
  end if

  call allocator%release_block(u)
  call MPI_Finalize(ierr)

  if (.not. test_pass) error stop 'Channel bulk velocity test failed.'

end program test_channel_bulk_velocity
