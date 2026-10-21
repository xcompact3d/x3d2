program test_transeq
  !! Verification of the transport-equation kernels through the backend
  !! interface: transeq_x for the velocity components and transeq_species for
  !! a passive scalar.  The domain is decomposed in x so that, with more than
  !! one rank, the distributed (halo-exchanging) path is exercised.
  !!
  !! With u = sin(x), v = cos(x), w = 0 and phi = cos(x), the skew-symmetric
  !! convection-diffusion term dphi = -1/2*(u*dphi/dx + d(u*phi)/dx)
  !! + nu*d2phi/dx2 gives
  !!   du   = -3/2*u*v - nu*u
  !!   dv   = u*u - 1/2*v*v - nu*v
  !!   dw   = 0
  !!   dphi = u*u - 1/2*phi*phi - nu*phi
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, DIR_Y, DIR_Z, VERT
  use m_mesh, only: mesh_t
  use m_tdsops_setup, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc

  ! The tolerance bounds the relative L2 error. The dw residual is measured
  ! against the du analytical field, since the dw reference is zero. Largest
  ! value measured over du, dv, dw and the species check:
  !   double precision, fixed tol 1e-8 (truncation dominated, the error
  !   falls on the finer grid):
  !     96 points: 9.6863e-10 (margin 10.3)
  !     120 points, 5 ranks: 2.6e-10 (margin 38)
  !     192 points, 8 ranks: 6.2e-11 (margin 160)
  !     384 points, 16 ranks: 4.5e-11
  !   single precision, OpenMP backend (gfortran), tol = K*eps/dx^2:
  !     96 points: 8.3407e-05, 192 points: 2.4383e-04
  ! The single precision error is roundoff dominated, ~eps/dx^2 from the
  ! second-derivative term, so it grows with the grid. With dx = 2*pi/n_glob
  ! the measured error is 3.0 eps/dx^2 at 96 points and 2.2 at 192. K is the
  ! smallest {1,2,5}x10^k value at least 10x the largest ratio, so K = 50.
  ! This gives tol = 1.39e-3 (margin 16.7) at 96 points and tol = 5.57e-3
  ! (margin 22.8) at 192 points.
#ifdef SINGLE_PREC
  real(dp), parameter :: tol_k = 50._dp
  real(dp) :: tol
#else
  real(dp), parameter :: tol = 1e-8_dp
#endif
  real(dp), parameter :: nu = 1._dp
  ! The global grid on which the tolerances for 1 to 4 ranks were measured.
  ! It divides by 1, 2, 3 and 4 ranks.
  integer, parameter :: n_glob_min = 96
  ! The distributed compact solve needs enough points per rank. No constant
  ! in src sets this limit: tdsops_init only warns when the entry that the
  ! distributed solver disregards exceeds a tolerance. The value is measured
  ! at 8 ranks in double precision (tol 1e-8, du error): 12 points per rank
  ! gives 8.9e-6 and 16 gives 1.6e-7, both fail. 20 gives 3.1e-9, which
  ! passes with a margin of only 3.2. 24 gives 6.2e-11 (margin 160). From 5
  ! ranks up the global grid is finer than 96, so the dp error is below the
  ! 96 point value (see above). The sp error grows with the roundoff term
  ! eps/dx^2, which the tolerance follows.
  integer, parameter :: min_points_per_rank = 24
  integer :: n_glob

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(dirps_t) :: xdirps

  class(field_t), pointer :: u, v, w, du, dv, dw
  real(dp), allocatable, dimension(:, :, :) :: u_data, v_data, du_data
  integer :: n, n_groups, dims(3)

  call initialise_mpi(nrank, nproc)
  n_glob = max(n_glob_min, min_points_per_rank*nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  mesh = mesh_t([n_glob, n_glob, n_glob], [nproc, 1, 1], &
                [2*pi, 2*pi, 2*pi], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])
#ifdef SINGLE_PREC
  tol = tol_k*epsilon(1.0_dp)/mesh%geo%d(DIR_X)**2
#endif

  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator
  if (nrank == 0) print *, trim(runtime%backend_name), ' backend instantiated'

  xdirps%dir = DIR_X
  call allocate_tdsops(xdirps, backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')

  n = mesh%get_n(DIR_X, VERT)
  n_groups = allocator%get_n_groups(DIR_X)

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  du => allocator%get_block(DIR_X, VERT)
  dv => allocator%get_block(DIR_X, VERT)
  dw => allocator%get_block(DIR_X, VERT)

  dims = u%get_shape()
  allocate (u_data(dims(1), dims(2), dims(3)))
  allocate (v_data(dims(1), dims(2), dims(3)))
  allocate (du_data(dims(1), dims(2), dims(3)))

  call fill_line(u_data, 'sin')
  call fill_line(v_data, 'cos')
  call backend%set_field_data(u, u_data, DIR_X)
  call backend%set_field_data(v, v_data, DIR_X)
  call w%fill(0._dp)

  ! Velocity components
  call backend%transeq_x(du, dv, dw, u, v, w, nu, xdirps)

  call backend%get_field_data(du_data, du, DIR_X)
  call check('transeq_x du', du_data, &
             -1.5_dp*u_data*v_data - nu*u_data)
  call backend%get_field_data(du_data, dv, DIR_X)
  call check('transeq_x dv', du_data, &
             u_data*u_data - 0.5_dp*v_data*v_data - nu*v_data)
  call backend%get_field_data(du_data, dw, DIR_X)
  call check('transeq_x dw', du_data, 0._dp*u_data, &
             -1.5_dp*u_data*v_data - nu*u_data)

  ! Passive scalar phi = cos(x), transported by u; v holds phi
  call backend%transeq_species(dv, u, v, nu, xdirps, .true.)

  call backend%get_field_data(du_data, dv, DIR_X)
  call check('transeq_species', du_data, &
             u_data*u_data - 0.5_dp*v_data*v_data - nu*v_data)

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)
  call allocator%release_block(du)
  call allocator%release_block(dv)
  call allocator%release_block(dw)

  call finalise_test(allpass, nrank)

contains

  subroutine fill_line(data, func)
    !! Fill a DIR_X array with a function of the global x index.
    real(dp), intent(out) :: data(:, :, :)
    character(*), intent(in) :: func

    real(dp) :: dx, x
    integer :: j, offset

    dx = mesh%geo%d(DIR_X)
    offset = mesh%par%nrank_dir(DIR_X)*n
    data = 0._dp
    do j = 1, n
      x = (j - 1 + offset)*dx
      select case (func)
      case ('sin')
        data(:, j, 1:n_groups) = sin(x)
      case ('cos')
        data(:, j, 1:n_groups) = cos(x)
      end select
    end do
  end subroutine fill_line

  subroutine check(label, computed, expected, ref)
    !! Relative L2 error of computed against expected.  The reference norm
    !! is that of expected unless ref is given, which is needed when expected
    !! is zero.
    character(*), intent(in) :: label
    real(dp), intent(in) :: computed(:, :, :), expected(:, :, :)
    real(dp), intent(in), optional :: ref(:, :, :)

    real(dp) :: norm

    if (present(ref)) then
      norm = relative_l2_error(computed(:, 1:n, 1:n_groups) &
                               - expected(:, 1:n, 1:n_groups), &
                               ref(:, 1:n, 1:n_groups))
    else
      norm = relative_l2_error(computed(:, 1:n, 1:n_groups) &
                               - expected(:, 1:n, 1:n_groups), &
                               expected(:, 1:n, 1:n_groups))
    end if

    if (nrank == 0) then
      print *, 'error norm ', label, norm
      if (norm > tol) then
        allpass = .false.
        write (stderr, '(a)') 'Check '//label//'... failed'
      else
        write (stderr, '(a)') 'Check '//label//'... passed'
      end if
    end if
  end subroutine check

end program test_transeq
