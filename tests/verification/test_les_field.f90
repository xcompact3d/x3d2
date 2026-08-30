program test_les_field
  use m_mpi, only: MPI_Finalize, MPI_Init

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT
  use m_config, only: les_config_t
  use m_field, only: field_t, flist_t
  use m_les, only: les_t, filter_width, horizontal_wall_velocity, &
                   neutral_wall_stress, &
                   wall_damped_mixing_length
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, allocate_tdsops, transeq_default
  use m_tdsops, only: dirps_t

#ifdef CUDA
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  type(mesh_t), target :: mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  type(dirps_t), target :: xdirps, ydirps, zdirps
  type(les_config_t) :: config
  type(les_t) :: les, les_wall
  type(solver_t) :: solver
  type(flist_t) :: rhs(3), variables(3)
  class(field_t), pointer :: u, v, w, rhs_u, rhs_v, rhs_w
  class(field_t), pointer :: wall_nut, dudy, dvdx, dwdy, dvdz

  integer, parameter :: dims_global(3) = [32, 33, 32]
  integer, parameter :: nproc_dir(3) = [1, 1, 1]
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  real(dp), parameter :: shear = 2.5_dp
  real(dp), parameter :: curvature = 1.75_dp
  character(len=9), parameter :: bc_periodic(2) = &
    ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_wall(2) = &
    ['dirichlet', 'dirichlet']

  real(dp), allocatable :: u_data(:, :, :), nut_data(:, :, :)
  real(dp) :: delta, expected, length, max_error, spacing(3), tolerance
  real(dp) :: sampling_height, roughness, kappa, drag_coeff
  real(dp) :: tau_x, tau_z
  real(dp) :: wall_u_mean, wall_w_mean
  integer :: dims(3), dims_padded(3), i, j, k, ierr
  logical :: all_pass

  call MPI_Init(ierr)
  all_pass = .true.
  tolerance = 2000._dp*epsilon(1._dp)

  mesh = mesh_t(dims_global, nproc_dir, lengths, &
                bc_periodic, bc_wall, bc_periodic)
  dims = mesh%get_dims(VERT)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(dims, SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(dims, SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  xdirps%dir = DIR_X
  ydirps%dir = DIR_Y
  zdirps%dir = DIR_Z
  call allocate_tdsops(xdirps, backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')
  call allocate_tdsops(ydirps, backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')
  call allocate_tdsops(zdirps, backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  call v%fill(0._dp)
  call w%fill(0._dp)

  dims_padded = allocator%get_padded_dims(DIR_C)
  allocate (u_data(dims_padded(1), dims_padded(2), dims_padded(3)))
  allocate (nut_data(dims_padded(1), dims_padded(2), dims_padded(3)))
  u_data = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        u_data(i, j, k) = shear*mesh%geo%vert_coords(j, 2)
      end do
    end do
  end do
  call backend%set_field_data(u, u_data)

  config%model = 'smagorinsky'
  config%smagorinsky_constant = 0.2_dp
  les = les_t(config)
  call les%compute_nut(backend, mesh, u, v, w, xdirps, ydirps, zdirps)
  call backend%get_field_data(nut_data, les%nut)

  spacing = mesh%geo%d
  delta = filter_width(spacing)
  expected = (config%smagorinsky_constant*delta)**2*abs(shear)
  max_error = maxval(abs(nut_data(1:dims(1), 1:dims(2), 1:dims(3)) - &
                         expected))
  call check_error('constant-shear nut', max_error, tolerance, all_pass)

  config%wall_damping = .true.
  config%roughness_length = 0._dp
  les_wall = les_t(config)
  call les_wall%compute_nut(backend, mesh, u, v, w, &
                            xdirps, ydirps, zdirps)
  call backend%get_field_data(nut_data, les_wall%nut)

  max_error = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      length = wall_damped_mixing_length( &
        delta, mesh%geo%vert_coords(j, 2), &
        config%smagorinsky_constant, config%von_karman_constant, &
        config%wall_damping_n, config%roughness_length)
      expected = length**2*abs(shear)
      do i = 1, dims(1)
        max_error = max(max_error, abs(nut_data(i, j, k) - expected))
      end do
    end do
  end do
  call check_error('wall-damped constant-shear nut', max_error, &
                   tolerance, all_pass)
  call check_error('smooth-wall nut', abs(nut_data(1, 1, 1)), &
                   tolerance, all_pass)

  config%wall_damping = .false.
  solver%backend => backend
  solver%mesh => mesh
  solver%xdirps => xdirps
  solver%ydirps => ydirps
  solver%zdirps => zdirps
  solver%les = les_t(config)
  solver%nu = 0._dp

  ! For u=(a/2)y^2, div(2*nut*S)_x=2*(Cs*Delta)^2*a^2*y.
  u_data = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        u_data(i, j, k) = &
          0.5_dp*curvature*mesh%geo%vert_coords(j, 2)**2
      end do
    end do
  end do
  call backend%set_field_data(u, u_data)

  rhs_u => allocator%get_block(DIR_X, VERT)
  rhs_v => allocator%get_block(DIR_X, VERT)
  rhs_w => allocator%get_block(DIR_X, VERT)
  call rhs_u%fill(0._dp)
  call rhs_v%fill(0._dp)
  call rhs_w%fill(0._dp)
  rhs(1)%ptr => rhs_u
  rhs(2)%ptr => rhs_v
  rhs(3)%ptr => rhs_w
  variables(1)%ptr => u
  variables(2)%ptr => v
  variables(3)%ptr => w
  call transeq_default(solver, rhs, variables)

  call backend%get_field_data(nut_data, rhs_u)
  max_error = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      expected = 2._dp*(config%smagorinsky_constant*delta)**2* &
                 curvature**2*mesh%geo%vert_coords(j, 2)
      do i = 1, dims(1)
        max_error = max(max_error, abs(nut_data(i, j, k) - expected))
      end do
    end do
  end do
  call check_error('solver-coupled SGS stress divergence', max_error, &
                   10._dp*tolerance, all_pass)

  call backend%get_field_data(nut_data, rhs_v)
  call check_error('SGS y-momentum cross-coupling', &
                   maxval(abs(nut_data(1:dims(1), 1:dims(2), 1:dims(3)))), &
                   10._dp*tolerance, all_pass)
  call backend%get_field_data(nut_data, rhs_w)
  call check_error('SGS z-momentum cross-coupling', &
                   maxval(abs(nut_data(1:dims(1), 1:dims(2), 1:dims(3)))), &
                   10._dp*tolerance, all_pass)

  ! ABL rough-wall correction. Use constant fields with prescribed shear
  ! components so the legacy second-order boundary formula is analytical.
  wall_nut => allocator%get_block(DIR_X, VERT)
  dudy => allocator%get_block(DIR_X, VERT)
  dvdx => allocator%get_block(DIR_X, VERT)
  dwdy => allocator%get_block(DIR_X, VERT)
  dvdz => allocator%get_block(DIR_X, VERT)

  ! Incompact3d samples one horizontally averaged velocity at dsampling*dy.
  ! Exercise the shared reduction with nonuniform columns; the field is
  ! constant in y, so the mean does not depend on which plane is named.
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        u_data(i, j, k) = real(i, dp)
        nut_data(i, j, k) = 2._dp*real(i, dp)
      end do
    end do
  end do
  call backend%set_field_data(u, u_data)
  call backend%set_field_data(w, nut_data)
  call horizontal_wall_velocity( &
    backend, mesh, u, w, 4, wall_u_mean, wall_w_mean)
  call check_error('ABL horizontally averaged wall u', &
                   abs(wall_u_mean - 0.5_dp*real(dims(1) + 1, dp)), &
                   tolerance, all_pass)
  call check_error('ABL horizontally averaged wall w', &
                   abs(wall_w_mean - real(dims(1) + 1, dp)), &
                   tolerance, all_pass)

  ! The wall stress now enters the SGS stress field rather than overwriting
  ! the stress divergence, so what is worth pinning here is the drag law:
  ! tau = (kappa/log(h/z0))^2 * u * |u|, signed positive for positive u so
  ! that differentiating it removes momentum.
  sampling_height = 3._dp*mesh%geo%d(2)
  roughness = 0.001_dp
  kappa = 0.4_dp
  drag_coeff = (kappa/log(sampling_height/roughness))**2
  call neutral_wall_stress(4._dp, 3._dp, kappa, roughness, sampling_height, &
                           tau_x, tau_z)
  call check_error('neutral wall stress x', &
                   abs(tau_x - drag_coeff*4._dp*5._dp), tolerance, all_pass)
  call check_error('neutral wall stress z', &
                   abs(tau_z - drag_coeff*3._dp*5._dp), tolerance, all_pass)
  ! A wall at rest exerts no stress.
  call neutral_wall_stress(0._dp, 0._dp, kappa, roughness, sampling_height, &
                           tau_x, tau_z)
  call check_error('no stress at zero velocity', &
                   abs(tau_x) + abs(tau_z), tolerance, all_pass)

  call les%finalise(backend)
  call les_wall%finalise(backend)
  call solver%finalise()
  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)
  call allocator%release_block(rhs_u)
  call allocator%release_block(rhs_v)
  call allocator%release_block(rhs_w)
  call allocator%release_block(wall_nut)
  call allocator%release_block(dudy)
  call allocator%release_block(dvdx)
  call allocator%release_block(dwdy)
  call allocator%release_block(dvdz)
  deallocate (u_data, nut_data)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  subroutine check_error(label, error, tol, pass)
    character(*), intent(in) :: label
    real(dp), intent(in) :: error, tol
    logical, intent(inout) :: pass

    if (error > tol) then
      print *, 'FAIL: ', label, ' error=', error, ' tolerance=', tol
      pass = .false.
    else
      print *, 'PASS: ', label, ' error=', error
    end if
  end subroutine check_error

end program test_les_field
