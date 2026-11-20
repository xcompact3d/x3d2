program test_laplacian_stretched
  !! Verification of the Laplacian on a mesh stretched in y, through the
  !! backend interface.
  !!
  !! On a stretched grid d2u/dy2 = vert_ds2*d2u/dxi2
  !! + vert_d2s*(vert_ds*du/dxi). A second derivative `tds_solve` applies
  !! only the vert_ds2 factor, so the Laplacian relies on
  !! `tds_solve_der2nd` for the vert_d2s term.
  !!
  !! With u = sin(kx x)*cos(ky y)*sin(kz z) the Laplacian is
  !! -(kx^2 + ky^2 + kz^2)*u. The mesh is the channel one (top-bottom
  !! stretching in y, beta of examples/channel/input.x3d), periodic in x and
  !! z and non-periodic in y. The domain splits over the ranks in x, and for
  !! several ranks again in y, the stretched direction. The test
  !! also evaluates the same field on a uniform mesh, where the correction
  !! must not change the result. Each case runs at two y resolutions and the
  !! error must fall on the finer one in double precision.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: field_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, VERT
  use m_mesh, only: mesh_t
  use m_tdsops_setup, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error
  use m_vector_calculus, only: vector_calculus_t

  implicit none

  logical :: allpass = .true.
  integer :: nrank, nproc

  ! The acceptance bounds are on the relative L2 error of the Laplacian.
  ! In double precision the stretched error at the finer resolution is below
  ! tol_stretched and falls by at least min_ratio from 65 to 129 points. The
  ! uniform error at 129 points is below tol_uniform and below the 65 point
  ! one. Measured on one and four ranks:
  !   stretched: 65 points 1.196e-3, 129 points 1.672e-5
  !   uniform: 65 points 7.499e-6, 129 points 3.350e-7
  ! In single precision the roundoff of the second derivative dominates and
  ! the error grows with the resolution, so there is no convergence check and
  ! both resolutions must be below the bound. Measured (OpenMP, gfortran):
  !   stretched: 65 points 3.07e-3, 129 points 8.25e-3
  !   uniform: 65 points 7.95e-5, 129 points 1.81e-4
  ! The y decomposed checks on four ranks use larger bounds, because the
  ! roundoff grows with ny squared and their finest mesh has 192 points:
  !   stretched: 96 points 7.8e-3, 192 points 3.15e-2
  !   uniform: 96 points 1.68e-4, 192 points 6.8e-4
  ! Without the stretched grid correction the stretched error is of order one.
#ifdef SINGLE_PREC
  real(dp), parameter :: tol_stretched = 2e-2_dp
  real(dp), parameter :: tol_uniform = 5e-4_dp
  real(dp), parameter :: tol_stretched_y = 6e-2_dp
  real(dp), parameter :: tol_uniform_y = 1.5e-3_dp
  logical, parameter :: must_converge = .false.
#else
  real(dp), parameter :: tol_stretched = 1e-3_dp
  real(dp), parameter :: tol_uniform = 1e-6_dp
  real(dp), parameter :: tol_stretched_y = tol_stretched
  real(dp), parameter :: tol_uniform_y = tol_uniform
  logical, parameter :: must_converge = .true.
#endif
  real(dp), parameter :: min_ratio = 32._dp
  ! beta of examples/channel/input.x3d
  real(dp), parameter :: stretch_beta = 0.259065151_dp
  real(dp), parameter :: kx = 1._dp, ky = 1._dp, kz = 1._dp
  integer, parameter :: ny_coarse = 65, ny_fine = 129, nz = 64
  ! The distributed compact solve needs enough points per rank.
  integer, parameter :: n_x_min = 96, min_points_per_rank = 24

  ! The y decomposed cases need ny to divide over the ranks and to give each
  ! rank at least min_points_per_rank points. The fine mesh doubles the
  ! coarse one, so the error falls as on the x decomposed cases.

  real(dp) :: err_stretched(2), err_uniform(2)
  integer, parameter :: ny_list(2) = [ny_coarse, ny_fine]
  integer :: ny_y(2)
  integer :: i

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  do i = 1, 2
    err_stretched(i) = laplacian_error(ny_list(i), .true., [nproc, 1, 1])
    err_uniform(i) = laplacian_error(ny_list(i), .false., [nproc, 1, 1])
  end do

  call check('stretched', ny_list, err_stretched, tol_stretched, .true.)
  call check('uniform', ny_list, err_uniform, tol_uniform, .false.)

  if (nproc > 1) then
    ny_y(1) = ((max(ny_coarse, min_points_per_rank*nproc) + nproc - 1)/nproc) &
              *nproc
    ny_y(2) = 2*ny_y(1)
    do i = 1, 2
      err_stretched(i) = laplacian_error(ny_y(i), .true., [1, nproc, 1])
      err_uniform(i) = laplacian_error(ny_y(i), .false., [1, nproc, 1])
    end do

    call check('stretched, y decomposed', ny_y, err_stretched, &
               tol_stretched_y, .true.)
    call check('uniform, y decomposed', ny_y, err_uniform, tol_uniform_y, &
               .false.)
  end if

  call finalise_test(allpass, nrank)

contains

  function laplacian_error(ny, stretched, decomp) result(err)
    !! Relative L2 error of the Laplacian of u on a mesh with ny points in y,
    !! decomposed over the ranks as decomp.
    integer, intent(in) :: ny
    logical, intent(in) :: stretched
    integer, intent(in) :: decomp(3)
    real(dp) :: err

    ! The objects are heap allocated. NVHPC never frees them, because its
    ! runtime crashes when they are finalised at the end of the function.
    ! Other compilers free them.
    type(backend_runtime_t), pointer :: runtime
    type(mesh_t), pointer :: mesh
    class(base_backend_t), pointer :: backend
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    type(vector_calculus_t), pointer :: vector_calculus
    class(field_t), pointer :: u, lapl_u
    real(dp), allocatable, dimension(:, :, :) :: u_data, lapl_data, &
                                                 lapl_expected
    real(dp) :: x, y, z
    integer :: nx, dims(3), dims_padded(3), ix, iy, iz

    nx = max(n_x_min, min_points_per_rank*decomp(1))

    allocate (runtime, mesh, xdirps, ydirps, zdirps, vector_calculus)

    if (stretched) then
      mesh = mesh_t([nx, ny, nz], decomp, [2*pi, 2*pi, 2*pi], &
                    ['periodic', 'periodic'], ['dirichlet', 'dirichlet'], &
                    ['periodic', 'periodic'], &
                    stretching=['uniform   ', 'top-bottom', 'uniform   '], &
                    beta=[1._dp, stretch_beta, 1._dp])
    else
      mesh = mesh_t([nx, ny, nz], decomp, [2*pi, 2*pi, 2*pi], &
                    ['periodic', 'periodic'], ['dirichlet', 'dirichlet'], &
                    ['periodic', 'periodic'])
    end if

    call runtime%init(mesh)
    backend => runtime%backend

    xdirps%dir = DIR_X
    ydirps%dir = DIR_Y
    zdirps%dir = DIR_Z
    call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')
    call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')
    call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                         'classic', 'compact6')

    vector_calculus = vector_calculus_t(backend)

    ! The field data arrays carry the padding of the allocator in the
    ! Cartesian orientation. Only the first dims points are the mesh.
    dims = mesh%get_dims(VERT)
    dims_padded = backend%allocator%get_padded_dims(DIR_C)
    allocate (u_data(dims_padded(1), dims_padded(2), dims_padded(3)))
    allocate (lapl_data(dims_padded(1), dims_padded(2), dims_padded(3)))
    allocate (lapl_expected(dims(1), dims(2), dims(3)))
    u_data = 0._dp

    do iz = 1, dims(3)
      do iy = 1, dims(2)
        do ix = 1, dims(1)
          x = mesh%geo%vert_coords(ix, 1)
          y = mesh%geo%vert_coords(iy, 2)
          z = mesh%geo%vert_coords(iz, 3)
          u_data(ix, iy, iz) = sin(kx*x)*cos(ky*y)*sin(kz*z)
        end do
      end do
    end do
    lapl_expected = -(kx**2 + ky**2 + kz**2) &
                    *u_data(1:dims(1), 1:dims(2), 1:dims(3))

    u => backend%allocator%get_block(DIR_X, VERT)
    lapl_u => backend%allocator%get_block(DIR_X, VERT)
    call backend%set_field_data(u, u_data)

    call vector_calculus%laplacian(lapl_u, u, xdirps%der1st, ydirps%der1st, &
                                   zdirps%der1st, xdirps%der2nd, &
                                   ydirps%der2nd, zdirps%der2nd)

    call backend%get_field_data(lapl_data, lapl_u)
    err = relative_l2_error(lapl_data(1:dims(1), 1:dims(2), 1:dims(3)) &
                            - lapl_expected, lapl_expected)

    call backend%allocator%release_block(u)
    call backend%allocator%release_block(lapl_u)

#ifndef __NVCOMPILER
    ! Free the field blocks that the allocator caches in its list.
    call runtime%allocator%destroy()
    deallocate (xdirps, ydirps, zdirps, vector_calculus, runtime, mesh)
#endif

  end function laplacian_error

  subroutine check(label, ny, err, tol, use_ratio)
    !! Checks the errors at the two resolutions ny. In single precision both
    !! must be below tol. In double precision the fine one must be below tol
    !! and fall on the finer mesh, by min_ratio if use_ratio is set. A
    !! non-finite error fails.
    character(*), intent(in) :: label
    integer, intent(in) :: ny(2)
    real(dp), intent(in) :: err(2), tol
    logical, intent(in) :: use_ratio

    logical :: ok

    if (nrank == 0) then
      print *, 'error norm ', label, ' ny =', ny(1), err(1)
      print *, 'error norm ', label, ' ny =', ny(2), err(2)
      ok = ieee_is_finite(err(1)) .and. ieee_is_finite(err(2))
      if (ok) then
        if (must_converge) then
          ok = err(2) < tol .and. err(2) < err(1)
          if (use_ratio) ok = ok .and. err(1) >= min_ratio*err(2)
        else
          ok = all(err < tol)
        end if
      end if
      if (ok) then
        write (stderr, '(a)') 'Check '//label//'... passed'
      else
        allpass = .false.
        write (stderr, '(a)') 'Check '//label//'... failed'
      end if
    end if
  end subroutine check

end program test_laplacian_stretched
