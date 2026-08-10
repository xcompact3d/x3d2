program test_les_field
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT
  use m_config, only: les_config_t
  use m_field, only: field_t
  use m_les, only: les_t, filter_width, wall_damped_mixing_length
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
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
  type(dirps_t) :: xdirps, ydirps, zdirps
  type(les_config_t) :: config
  type(les_t) :: les, les_wall
  class(field_t), pointer :: u, v, w

  integer, parameter :: dims_global(3) = [32, 33, 32]
  integer, parameter :: nproc_dir(3) = [1, 1, 1]
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  real(dp), parameter :: shear = 2.5_dp
  character(len=9), parameter :: bc_periodic(2) = &
    ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_wall(2) = &
    ['dirichlet', 'dirichlet']

  real(dp), allocatable :: u_data(:, :, :), nut_data(:, :, :)
  real(dp) :: delta, expected, length, max_error, spacing(3), tolerance
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

  call les%finalise(backend)
  call les_wall%finalise(backend)
  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)
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
