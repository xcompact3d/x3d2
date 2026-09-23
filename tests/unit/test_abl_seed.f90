program test_abl_seed
  !! A fixed ABL noise seed must repeat the initial field exactly, and a
  !! different seed must change it. The noise must stay within init_noise.
  use m_mpi, only: MPI_Init, MPI_Finalize

  use m_abl, only: abl_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_C, VERT
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t

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
  type(allocator_t), target :: host_allocator
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(omp_backend_t), target :: omp_backend
#endif
  class(field_t), pointer :: u, v, w

  integer, parameter :: dims_global(3) = [16, 17, 8]
  integer, parameter :: nproc_dir(3) = [1, 1, 1]
  real(dp), parameter :: lengths(3) = [1000._dp, 500._dp, 1000._dp]
  real(dp), parameter :: noise = 0.5_dp
  character(len=9), parameter :: bc_periodic(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_wall(2) = ['dirichlet', 'neumann  ']

  real(dp), allocatable :: first(:, :, :, :), second(:, :, :, :)
  integer :: dims(3), dims_padded(3), ierr
  logical :: all_pass

  call MPI_Init(ierr)
  all_pass = .true.

  mesh = mesh_t(dims_global, nproc_dir, lengths, &
                bc_periodic, bc_wall, bc_periodic)
  dims = mesh%get_dims(VERT)

  host_allocator = allocator_t(dims, SZ)
#ifdef CUDA
  cuda_allocator = cuda_allocator_t(dims, SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  allocator => host_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  dims_padded = allocator%get_padded_dims(DIR_C)
  allocate (first(dims_padded(1), dims_padded(2), dims_padded(3), 3))
  allocate (second, mold=first)

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)

  call initial_field(42, first)
  call initial_field(42, second)
  call check('same seed repeats the field', &
             maxval(abs(active(first) - active(second))) == 0._dp, all_pass)

  call initial_field(43, second)
  call check('a different seed changes the field', &
             maxval(abs(active(first) - active(second))) > 0._dp, all_pass)

  ! v and w are pure noise, bounded by init_noise
  call check('noise stays within init_noise', &
             maxval(abs(active(first(:, :, :, 2:3)))) <= noise, all_pass)

  ! seed = 0 draws from the clock and must still give a bounded field
  call initial_field(0, second)
  call check('clock seed stays within init_noise', &
             maxval(abs(active(second(:, :, :, 2:3)))) <= noise, all_pass)

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  subroutine initial_field(seed, fields)
    !! Initialise the ABL velocity with the given seed and return u, v, w.
    integer, intent(in) :: seed
    real(dp), intent(out) :: fields(:, :, :, :)

    type(abl_config_t) :: cfg
    type(abl_t) :: abl
    character(len=256) :: nml

    write (nml, '(a,i0,a)') '&abl_nml u_star = 0.45, kappa = 0.4, &
      &init_noise = 3*0.5, pressure_gradient = .true., seed = ', seed, ' /'
    call cfg%read(nml_string=nml)

    abl = abl_t(backend, mesh, host_allocator, cfg, 1._dp)
    call abl%initialise(u, v, w)
    call backend%get_field_data(fields(:, :, :, 1), u)
    call backend%get_field_data(fields(:, :, :, 2), v)
    call backend%get_field_data(fields(:, :, :, 3), w)
  end subroutine initial_field

  function active(fields) result(inner)
    !! The unpadded part of the fields
    real(dp), intent(in) :: fields(:, :, :, :)
    real(dp), allocatable :: inner(:, :, :, :)

    inner = fields(1:dims(1), 1:dims(2), 1:dims(3), :)
  end function active

  subroutine check(label, condition, pass)
    character(*), intent(in) :: label
    logical, intent(in) :: condition
    logical, intent(inout) :: pass

    if (condition) then
      print *, 'PASS: ', label
    else
      print *, 'FAIL: ', label
      pass = .false.
    end if
  end subroutine check

end program test_abl_seed
