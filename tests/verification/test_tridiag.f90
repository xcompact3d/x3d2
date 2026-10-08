program test_tridiag
  !! Verification of the distributed tridiagonal solver through the backend
  !! interface (backend%tds_solve) for the compact first/second derivative,
  !! staggered interpolation/derivative and hyperviscous operators.  The
  !! domain is decomposed in x, so with more than one rank the boundary
  !! conditions at interior subdomain edges are BC_HALO.
  use iso_fortran_env, only: stderr => error_unit
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_sz
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_X, VERT, &
                      BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, BC_HALO
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc, ierr

  integer, parameter :: n_glob = 1024
  ! Single precision roundoff floors at n_glob=1024: ~eps/dx (~2e-5) for
  ! first derivatives/interpolation, ~eps/dx^2 (~3e-3) for second
  ! derivatives, and ~63x more for the hyperviscous operator (nu0_nu=63).
#ifdef SINGLE_PREC
  real(dp), parameter :: tol = 1e-4_dp, tol_2nd = 2e-2_dp, tol_hyper = 2.0_dp
#else
  real(dp), parameter :: tol = 1e-8_dp, tol_2nd = 1e-8_dp, &
                         tol_hyper = 1e-8_dp
#endif

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  ! One operator per case, never deallocated: NVHPC (25.3) crashes in
  ! pgf90_dealloc_poly03 when a class(tdsops_t) holding a cuda_tdsops_t is
  ! deallocated.
  type :: tdsops_slot_t
    class(tdsops_t), allocatable :: op
  end type tdsops_slot_t
  type(tdsops_slot_t) :: ops(8)
  integer :: icase = 0
  class(field_t), pointer :: u_field, du_field

  real(dp), allocatable, dimension(:, :, :) :: u, du
  real(dp), allocatable, dimension(:) :: sin_0_2pi_per, cos_0_2pi_per, &
                                         sin_0_2pi, cos_0_2pi, &
                                         cos_0_pi, cos_0_pi_stag, &
                                         sin_0_pi
  real(dp) :: dx, dx_per, dx_pi
  integer :: n, n_loc, n_groups, bc_start, bc_end, dims(3), j, offset

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  n_groups = 64*64/backend_sz
  mesh = mesh_t([n_glob, backend_sz, n_groups], [nproc, 1, 1], &
                [2*pi, 1._dp, 1._dp], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])

  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator
  if (nrank == 0) print *, trim(runtime%backend_name), ' backend instantiated'

  n = mesh%get_n(DIR_X, VERT)
  offset = mesh%par%nrank_dir(DIR_X)*n

  u_field => allocator%get_block(DIR_X)
  du_field => allocator%get_block(DIR_X)
  dims = u_field%get_shape()
  allocate (u(dims(1), dims(2), dims(3)), du(dims(1), dims(2), dims(3)))
  u = 0._dp

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)
  dx_pi = pi/(n_glob - 1)

  allocate (sin_0_2pi_per(n), cos_0_2pi_per(n), sin_0_2pi(n), cos_0_2pi(n))
  allocate (cos_0_pi(n), cos_0_pi_stag(n), sin_0_pi(n))
  do j = 1, n
    sin_0_2pi_per(j) = sin((j - 1 + offset)*dx_per)
    cos_0_2pi_per(j) = cos((j - 1 + offset)*dx_per)
    sin_0_2pi(j) = sin((j - 1 + offset)*dx)
    cos_0_2pi(j) = cos((j - 1 + offset)*dx)
    cos_0_pi(j) = cos((j - 1 + offset)*dx_pi)
    cos_0_pi_stag(j) = cos((j - 1 + offset)*dx_pi + dx_pi/2._dp)
    sin_0_pi(j) = sin((j - 1 + offset)*dx_pi)
  end do

  ! ===========================================================================
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx_per, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  call run_case(sin_0_2pi_per, n, n, sin_0_2pi_per, 1, tol_2nd, &
                '2nd derivatives, periodic BCs')

  ! ===========================================================================
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx_per, operation='first-deriv', &
                            scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  call run_case(sin_0_2pi_per, n, n, cos_0_2pi_per, -1, tol, &
                '1st derivatives, periodic BCs')

  ! ===========================================================================
  ! Non-periodic cases: physical BCs at the ends of the global domain only
  bc_start = BC_HALO
  bc_end = BC_HALO
  if (nrank == 0) bc_start = BC_DIRICHLET
  if (nrank == nproc - 1) bc_end = BC_NEUMANN

  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx, operation='first-deriv', &
                            scheme='compact6', &
                            bc_start=bc_start, bc_end=bc_end, sym=.false.)
  call run_case(sin_0_2pi, n, n, cos_0_2pi, -1, tol, &
                '1st derivatives, dir-neu')

  ! ===========================================================================
  ! Staggered operators: the last subdomain holds one fewer cell centre
  if (nrank == 0) bc_start = BC_NEUMANN
  n_loc = n
  if (nrank == nproc - 1) n_loc = n - 1

  ! stag-interpolate v2p requires an even, cos-type function
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n_loc, dx_pi, operation='interpolate', &
                            scheme='classic', &
                            bc_start=bc_start, bc_end=bc_end, from_to='v2p')
  call run_case(cos_0_pi, n, n_loc, cos_0_pi_stag, -1, tol, &
                'interpolation "v2p"')

  ! stag-interpolate p2v requires an even, cos-type function
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx_pi, operation='interpolate', &
                            scheme='classic', &
                            bc_start=bc_start, bc_end=bc_end, from_to='p2v')
  call run_case(cos_0_pi_stag, n_loc, n, cos_0_pi, -1, tol, &
                'interpolation "p2v"')

  ! stag-derivative v2p requires an odd, sin-type function
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n_loc, dx_pi, operation='stag-deriv', &
                            scheme='compact6', &
                            bc_start=bc_start, bc_end=bc_end, from_to='v2p')
  call run_case(sin_0_pi, n, n_loc, cos_0_pi_stag, -1, tol, &
                'stag derivative "v2p"')

  ! stag-derivative p2v requires an even, cos-type function
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx_pi, operation='stag-deriv', &
                            scheme='compact6', &
                            bc_start=bc_start, bc_end=bc_end, from_to='p2v')
  call run_case(cos_0_pi_stag, n_loc, n, sin_0_pi, 1, tol, &
                'stag derivative "p2v"')

  ! ===========================================================================
  ! c_nu = 0.22 and nu0_nu = 63 results in alpha = 0.40869111947709036
  icase = icase + 1
  call backend%alloc_tdsops(ops(icase)%op, n, dx, operation='second-deriv', &
                            scheme='compact6-hyperviscous', &
                            bc_start=bc_start, bc_end=bc_end, &
                            sym=.false., c_nu=0.22_dp, nu0_nu=63._dp)
  call run_case(sin_0_2pi, n, n, sin_0_2pi, 1, tol_hyper, &
                '2nd ders, hyperviscous, dir-neu')

  call allocator%release_block(u_field)
  call allocator%release_block(du_field)

  call finalise_test(allpass, nrank)

contains

  subroutine run_case(input, n_in, n_out, expected, c, tolerance, label)
    !! Solve with the latest operator and check du + c*expected ~ 0 over the
    !! first n_out points.  Only the first n_in points of u are overwritten,
    !! which leaves the rest from the previous case as the original
    !! kernel-level test did.
    real(dp), intent(in) :: input(:), expected(:)
    integer, intent(in) :: n_in, n_out, c
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label

    real(dp) :: norm
    integer :: jj

    do jj = 1, n_in
      u(:, jj, 1:n_groups) = input(jj)
    end do
    call backend%set_field_data(u_field, u, DIR_X)

    call backend%tds_solve(du_field, u_field, ops(icase)%op)
    call backend%get_field_data(du, du_field, DIR_X)

    do jj = 1, n_out
      du(:, jj, 1:n_groups) = du(:, jj, 1:n_groups) + c*expected(jj)
    end do

    norm = sum(du(:, 1:n_out, 1:n_groups)**2)
    norm = norm/real(n_glob, dp)/real(n_groups*backend_sz, dp)
    call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    norm = sqrt(norm)

    if (nrank == 0) then
      print *, 'error norm ', label, norm
      if (norm > tolerance) then
        allpass = .false.
        write (stderr, '(a)') 'Check '//label//'... failed'
      else
        write (stderr, '(a)') 'Check '//label//'... passed'
      end if
    end if
  end subroutine run_case

end program test_tridiag
