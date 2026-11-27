program test_monitoring_volume_mean
  !! Verifies the volume mean of a squared vector norm that the monitoring
  !! uses for the kinetic energy and the enstrophy.
  !!
  !! The fields on 0 <= y <= 2 are
  !!   a = y*(2 - y)*(1 + cos(2*pi*x/Lx)),
  !!   b = y**2*(2 - y),
  !!   c = y*(2 - y)*sin(2*pi*z/Lz).
  !! b is not symmetric in y, so a wrong reduction over the y decomposition
  !! changes the result. The exact mean of a**2 + b**2 + c**2 is 176/105.
  !! The y direction is split over all ranks, so every rank uses its own slice
  !! of vert_ds.
  !!
  !! The mesh has y stretching. Check 1 compares volume_mean_norm_squared
  !! with a host reference that has its own MPI_Allreduce and does not use the
  !! backend. Check 2 compares it with the exact value. The uniform mesh path
  !! is vector_norm_squared over the point count, which test_kinetic_energy
  !! covers.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_INTEGER, MPI_LAND, &
                   MPI_LOGICAL, MPI_SUM, MPI_Allreduce, MPI_Comm_rank, &
                   MPI_Comm_size, MPI_Finalize, MPI_Init

  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_C, DIR_X, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_monitoring, only: volume_mean_norm_squared

  implicit none

  ! ny divides evenly over 1 and 2 ranks, so both walls are present.
  integer, parameter :: nx = 16, ny = 64, nz = 16
  ! beta of the channel smoke case
  real(dp), parameter :: beta = 0.259065151_dp
  real(dp), parameter :: mean_exact = 176._dp/105._dp
  ! Check 1: the two sums differ only in the order of the additions.
  ! The sums have nx*ny*nz terms of size at most 2.
  real(dp), parameter :: tol_ref = 2._dp*real(nx*ny*nz, dp)*epsilon(1._dp)
  ! Check 2: the integrand and its y derivative vanish at both walls, so the
  ! end corrections of the trapezoid rule vanish to order h**2. 2*h**2, with
  ! h = 1/(ny - 1), is a loose bound.
  real(dp), parameter :: tol_quad = 2._dp/real(ny - 1, dp)**2
  character(len=10), parameter :: uni = 'uniform', tb = 'top-bottom'
  character(len=8), parameter :: periodic(2) = ['periodic', 'periodic']
  character(len=9), parameter :: walls(2) = ['dirichlet', 'dirichlet']

  real(dp) :: mean, mean_ref, err_ref, err_quad
  integer :: nvert, ierr, nrank, nproc
  logical :: test_pass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  call run_case(mean, mean_ref, nvert)
  err_ref = abs(mean - mean_ref)
  err_quad = abs(mean - mean_exact)
  test_pass = nvert == ny
  test_pass = test_pass .and. finite(mean) .and. finite(mean_ref)
  test_pass = test_pass .and. err_ref <= tol_ref*mean_exact
  test_pass = test_pass .and. err_quad < tol_quad

  if (nrank == 0) then
    print '(a,i0,a,i0)', 'ranks: ', nproc, ', vertex planes: ', nvert
    print '(a,es25.16e3)', 'mean:            ', mean
    print '(a,es25.16e3)', 'host ref:        ', mean_ref
    print '(a,es25.16e3,a,es25.16e3)', '|mean - ref|:    ', err_ref, &
      '  tol: ', tol_ref*mean_exact
    print '(a,es25.16e3,a,es25.16e3)', '|mean - exact|:  ', err_quad, &
      '  tol: ', tol_quad
  end if

  call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, MPI_LOGICAL, MPI_LAND, &
                     MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    if (test_pass) then
      print *, 'PASS'
    else
      print *, 'FAIL'
    end if
  end if

  call MPI_Finalize(ierr)

  if (.not. test_pass) error stop 'Monitoring volume mean test failed.'

contains

  logical function finite(x)
    real(dp), intent(in) :: x

    finite = ieee_is_finite(x) .and. .not. ieee_is_nan(x)
  end function finite

  subroutine run_case(mean, mean_ref, nvert)
    real(dp), intent(out) :: mean, mean_ref
    integer, intent(out) :: nvert

    type(mesh_t), target :: mesh
    type(backend_runtime_t), target :: runtime
    class(allocator_t), pointer :: allocator
    class(base_backend_t), pointer :: backend
    class(field_t), pointer :: a, b, c
    real(dp), allocatable :: host_a(:, :, :), host_b(:, :, :), &
                             host_c(:, :, :)
    real(dp) :: x, y, z, plane
    integer :: dims(3), dims_padded(3), i, j, k

    mesh = mesh_t([nx, ny, nz], [1, nproc, 1], [4._dp, 2._dp, 2._dp], &
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

    ! The field blocks are padded, so the host arrays are padded too. The
    ! padding stays zero and the sums only read the unpadded range.
    dims_padded = allocator%get_padded_dims(DIR_C)
    allocate (host_a(dims_padded(1), dims_padded(2), dims_padded(3)))
    allocate (host_b, host_c, mold=host_a)
    host_a = 0._dp
    host_b = 0._dp
    host_c = 0._dp
    do k = 1, dims(3)
      z = mesh%geo%vert_coords(k, 3)
      do j = 1, dims(2)
        y = mesh%geo%vert_coords(j, 2)
        do i = 1, dims(1)
          x = mesh%geo%vert_coords(i, 1)
          host_a(i, j, k) = y*(2._dp - y) &
                            *(1._dp + cos(2._dp*pi*x/mesh%geo%L(1)))
          host_b(i, j, k) = y**2*(2._dp - y)
          host_c(i, j, k) = y*(2._dp - y)*sin(2._dp*pi*z/mesh%geo%L(3))
        end do
      end do
    end do

    a => allocator%get_block(DIR_X, VERT)
    b => allocator%get_block(DIR_X, VERT)
    c => allocator%get_block(DIR_X, VERT)
    call backend%set_field_data(a, host_a)
    call backend%set_field_data(b, host_b)
    call backend%set_field_data(c, host_c)

    mean = volume_mean_norm_squared(backend, mesh, a, b, c)

    ! Host reference: each local plane mean has the weight
    ! d(2)/(vert_ds*L(2)).
    mean_ref = 0._dp
    do j = 1, dims(2)
      plane = sum(host_a(1:dims(1), j, 1:dims(3))**2 &
                  + host_b(1:dims(1), j, 1:dims(3))**2 &
                  + host_c(1:dims(1), j, 1:dims(3))**2) &
              /real(nx*nz, dp)
      mean_ref = mean_ref &
                 + plane*mesh%geo%d(2)/(mesh%geo%vert_ds(j, 2)*mesh%geo%L(2))
    end do
    call MPI_Allreduce(MPI_IN_PLACE, mean_ref, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

    call allocator%release_block(a)
    call allocator%release_block(b)
    call allocator%release_block(c)

  end subroutine run_case

end program test_monitoring_volume_mean
