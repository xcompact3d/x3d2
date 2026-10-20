program test_tridiag
  !! Verification of the distributed tridiagonal solver through the backend
  !! interface (backend%tds_solve) for the compact first/second derivative,
  !! staggered interpolation/derivative and hyperviscous operators.  The
  !! domain is decomposed in x, so with more than one rank the boundary
  !! conditions at interior subdomain edges are BC_HALO.
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_sz, tdsops_slot_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, VERT, &
                      BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, BC_HALO
  use m_logger, only: abort_run
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc

  ! The generic decomposition splits x evenly, so n_glob is rounded down to a
  ! multiple of the rank count (1024 on a power-of-two count).
  integer :: n_glob
  ! The hyperviscous case runs on a coarser grid of (256/nproc)*nproc points,
  ! see the tolerance comment.
  integer :: n_hyper_glob
  ! The tolerances bound the relative L2 error. Single precision roundoff
  ! floors at n_glob~1024: ~eps/dx (~2e-5) for first derivatives/interpolation,
  ! ~eps/dx^2 (~3e-3) for second derivatives, and ~63x more for the
  ! hyperviscous operator (nu0_nu=63). The hyperviscous case therefore runs on
  ! a coarser grid of (256/nproc)*nproc points (255 on 3 ranks). Its single
  ! precision error is roundoff dominated and grows ~4x per grid doubling
  ! (7.2e-3 at 128 points, 2.7e-2 at 256, 0.117 at 512), which is why the
  ! coarse grid is used.
  ! Largest values measured over 1, 2, 3, 4 and 8 ranks (gfortran and
  ! nvfortran, Release and Debug), with the rank count of the maximum
  ! (tol_hyper at the coarse grid):
  !   double precision, OpenMP backend:
  !     tol 1.2734e-11 (3 ranks),
  !     tol_2nd 1.9512e-11 (1 or 4 ranks, earlier runs),
  !     tol_hyper 4.2499e-10 (8 ranks, 32 points per rank)
  !   double precision, CUDA backend (1 and 4 ranks):
  !     tol 1.2678e-11, tol_2nd 1.9422e-11, tol_hyper 6.47e-11
  !   single precision, OpenMP backend:
  !     tol 4.0779e-05 (8 ranks, Debug),
  !     tol_2nd 1.1822e-02 (3 ranks, n_glob = 1023),
  !     tol_hyper 2.7724e-02 (2 ranks, Debug)
  !   single precision, CUDA backend (1 and 4 ranks):
  !     tol 3.8401e-05, tol_2nd 4.5822e-03, tol_hyper 2.707e-02
  ! The double precision tol_hyper error grows with the rank count: ~6.5e-11
  ! on 2 ranks, 9.2e-11 on 3, 4.2499e-10 on 8.
  ! Each tolerance is the smallest {1,2,5}x10^k value at least 10x the
  ! largest of these, a margin of 10 or more in every case. The single
  ! precision tol_hyper is below 1, so it rejects zero output (error 1) and a
  ! sign flip (error 2).
#ifdef SINGLE_PREC
  real(dp), parameter :: tol = 5e-4_dp, tol_2nd = 2e-1_dp, tol_hyper = 5e-1_dp
#else
  real(dp), parameter :: tol = 2e-10_dp, tol_2nd = 2e-10_dp, &
                         tol_hyper = 5e-9_dp
#endif

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  ! One operator per case, never deallocated (see tdsops_slot_t)
  type(tdsops_slot_t) :: ops(8)
  integer :: icase = 0
  ! Optional first command-line argument: 'zero_output' or 'sign_flip' damages
  ! the hyperviscous result, so ctest can check that tol_hyper rejects it.
  character(len=16) :: mode = ''
  class(field_t), pointer :: u_field, du_field

  real(dp), allocatable, dimension(:, :, :) :: u, du, ref
  real(dp), allocatable, dimension(:) :: sin_0_2pi_per, cos_0_2pi_per, &
                                         sin_0_2pi, cos_0_2pi, &
                                         cos_0_pi, cos_0_pi_stag, &
                                         sin_0_pi, sin_hyper
  real(dp) :: dx, dx_per, dx_pi, dx_hyper
  integer :: n, n_loc, n_hyper, n_groups, bc_start, bc_end, dims(3), j, &
             offset, offset_hyper

  call initialise_mpi(nrank, nproc)
  if (command_argument_count() >= 1) call get_command_argument(1, mode)
  if (len_trim(mode) > 0 .and. trim(mode) /= 'zero_output' &
      .and. trim(mode) /= 'sign_flip') then
    call abort_run('test_tridiag: unknown mode '//trim(mode), &
                   collective=.true.)
  end if
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  n_glob = (1024/nproc)*nproc
  n_hyper_glob = (256/nproc)*nproc
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
  n_hyper = n_hyper_glob/nproc
  offset_hyper = mesh%par%nrank_dir(DIR_X)*n_hyper

  u_field => allocator%get_block(DIR_X)
  du_field => allocator%get_block(DIR_X)
  dims = u_field%get_shape()
  allocate (u(dims(1), dims(2), dims(3)), du(dims(1), dims(2), dims(3)), &
            ref(dims(1), dims(2), dims(3)))
  u = 0._dp

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)
  dx_pi = pi/(n_glob - 1)
  dx_hyper = 2*pi/(n_hyper_glob - 1)

  allocate (sin_0_2pi_per(n), cos_0_2pi_per(n), sin_0_2pi(n), cos_0_2pi(n))
  allocate (cos_0_pi(n), cos_0_pi_stag(n), sin_0_pi(n), sin_hyper(n_hyper))
  do j = 1, n
    sin_0_2pi_per(j) = sin((j - 1 + offset)*dx_per)
    cos_0_2pi_per(j) = cos((j - 1 + offset)*dx_per)
    sin_0_2pi(j) = sin((j - 1 + offset)*dx)
    cos_0_2pi(j) = cos((j - 1 + offset)*dx)
    cos_0_pi(j) = cos((j - 1 + offset)*dx_pi)
    cos_0_pi_stag(j) = cos((j - 1 + offset)*dx_pi + dx_pi/2._dp)
    sin_0_pi(j) = sin((j - 1 + offset)*dx_pi)
  end do
  do j = 1, n_hyper
    sin_hyper(j) = sin((j - 1 + offset_hyper)*dx_hyper)
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
  call backend%alloc_tdsops(ops(icase)%op, n_hyper, dx_hyper, &
                            operation='second-deriv', &
                            scheme='compact6-hyperviscous', &
                            bc_start=bc_start, bc_end=bc_end, &
                            sym=.false., c_nu=0.22_dp, nu0_nu=63._dp)
  call run_case(sin_hyper, n_hyper, n_hyper, sin_hyper, 1, tol_hyper, &
                '2nd ders, hyperviscous, dir-neu', damage=.true.)

  call allocator%release_block(u_field)
  call allocator%release_block(du_field)

  call finalise_test(allpass, nrank)

contains

  subroutine run_case(input, n_in, n_out, expected, c, tolerance, label, damage)
    !! Solve with the latest operator and check du + c*expected ~ 0 over the
    !! first n_out points.  Only the first n_in points of u are overwritten,
    !! which leaves the rest from the previous case as the original
    !! kernel-level test did.
    real(dp), intent(in) :: input(:), expected(:)
    integer, intent(in) :: n_in, n_out, c
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label
    logical, optional, intent(in) :: damage

    real(dp) :: norm
    integer :: jj

    do jj = 1, n_in
      u(:, jj, 1:n_groups) = input(jj)
    end do
    call backend%set_field_data(u_field, u, DIR_X)

    call backend%tds_solve(du_field, u_field, ops(icase)%op)
    call backend%get_field_data(du, du_field, DIR_X)

    if (present(damage)) then
      if (damage .and. trim(mode) == 'zero_output') du = 0._dp
      if (damage .and. trim(mode) == 'sign_flip') du = -du
    end if

    do jj = 1, n_out
      du(:, jj, 1:n_groups) = du(:, jj, 1:n_groups) + c*expected(jj)
      ref(:, jj, 1:n_groups) = expected(jj)
    end do

    norm = relative_l2_error(du(:, 1:n_out, 1:n_groups), &
                             ref(:, 1:n_out, 1:n_groups))

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
