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
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_X, DIR_Y, DIR_Z, VERT
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc, ierr

  ! Roundoff floor of the second-derivative term is ~eps/dx^2, which at
  ! 96 points per period is ~1e-4 in single precision.
#ifdef SINGLE_PREC
  real(dp), parameter :: tol = 5e-3_dp
#else
  real(dp), parameter :: tol = 1e-8_dp
#endif
  real(dp), parameter :: nu = 1._dp
  integer, parameter :: n_glob = 96

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(dirps_t) :: xdirps

  class(field_t), pointer :: u, v, w, du, dv, dw
  real(dp), allocatable, dimension(:, :, :) :: u_data, v_data, du_data
  integer :: n, n_groups, dims(3)

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  mesh = mesh_t([n_glob, n_glob, n_glob], [nproc, 1, 1], &
                [2*pi, 2*pi, 2*pi], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])

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
  call check('transeq_x dw', du_data, 0._dp*u_data)

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

  subroutine check(label, computed, expected)
    character(*), intent(in) :: label
    real(dp), intent(in) :: computed(:, :, :), expected(:, :, :)

    real(dp) :: norm

    norm = sum((computed(:, 1:n, 1:n_groups) &
                - expected(:, 1:n, 1:n_groups))**2)
    norm = norm/real(n_glob, dp)/real(n_groups*size(computed, 1), dp)
    call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    norm = sqrt(norm)

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
