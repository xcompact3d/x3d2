program test_omp_transeq_species
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_X, DIR_Y, DIR_Z, VERT
  use m_omp_common, only: SZ
  use m_omp_backend, only: omp_backend_t
  use m_tdsops, only: dirps_t
  use m_solver, only: allocate_tdsops
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  logical :: allpass = .true.
  class(field_t), pointer :: u, v, w, spec
  class(field_t), pointer :: dspec
  real(dp), dimension(:, :, :), allocatable :: r_u
  class(mesh_t), allocatable :: mesh

  integer, dimension(3) :: dims_global
  integer :: n, n_groups
  integer :: nrank, nproc
  integer :: ierr

  real(dp) :: dx_per, nu, norm_du
  ! Roundoff floor of the second-derivative term is ~eps/dx^2, which at
  ! 96^3 cells is ~1e-4 in single precision.
#ifdef SINGLE_PREC
  real(dp) :: tol = 5e-3
#else
  real(dp) :: tol = 1d-8
#endif

  class(allocator_t), pointer :: allocator

  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
  type(dirps_t) :: xdirps, ydirps, zdirps

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'
  call setup_mesh()
  call setup_backend()
  call initialise_input()
  call run_kernel()
  call check_result()
  call finalise_test(allpass, nrank)

contains


  subroutine setup_mesh()
    integer, dimension(3) :: nproc_dir
    real(dp), dimension(3) :: L_global
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)

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
  end subroutine setup_mesh

  subroutine setup_backend()
    xdirps%dir = DIR_X; ydirps%dir = DIR_Y; zdirps%dir = DIR_Z

    omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
    allocator => omp_allocator
    print *, 'OpenMP allocator instantiated'

    omp_backend = omp_backend_t(mesh, allocator)
    print *, 'OpenMP backend instantiated'

    n = mesh%get_n(DIR_X, VERT)
    n_groups = allocator%get_n_groups(DIR_X)
    dx_per = mesh%geo%d(DIR_X)
    nu = 1._dp

    call allocate_tdsops(xdirps, omp_backend, mesh, &
                         'compact6', 'compact6', 'classic', 'compact6')
  end subroutine setup_backend

  subroutine initialise_input()
    integer :: i, j, k

    u => allocator%get_block(DIR_X, VERT)
    v => allocator%get_block(DIR_X, VERT)
    w => allocator%get_block(DIR_X, VERT)

    spec => allocator%get_block(DIR_X, VERT)
    dspec => allocator%get_block(DIR_X, VERT)

    do k = 1, n_groups
      do j = 1, n
        do i = 1, SZ
          u%data(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
          v%data(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
          spec%data(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
        end do
      end do
    end do
    w%data(:, :, :) = 0.d0
  end subroutine initialise_input

  subroutine run_kernel()
    ! Compute species convection-diffusion
    call omp_backend%transeq_species(dspec, u, spec, &
                                     nu, xdirps, .true.)
  end subroutine run_kernel

  subroutine check_result()
    allocate (r_u(SZ, n, n_groups))

    ! check error
    ! dspec = -1/2*(u*dspec/dx + d(u*spec)/dx) + nu*d2spec/dx2
    ! u is sin(x), v is cos(x), spec is cos(x);
    ! dspec = -1/2*(u*(-u) + spec*spec + u*(-u)) + nu*(-spec)
    !       = u*u - 1/2*spec*spec - nu*spec
    r_u = dspec%data &
          - (u%data*u%data - 0.5_dp*spec%data*spec%data - nu*spec%data)
    norm_du = norm2(r_u)
    norm_du = norm_du*norm_du/dims_global(DIR_X)/n_groups/SZ
    call MPI_Allreduce(MPI_IN_PLACE, norm_du, 1, MPI_X3D2_DP, &
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

  end subroutine check_result

end program test_omp_transeq_species
