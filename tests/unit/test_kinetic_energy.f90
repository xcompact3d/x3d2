program test_kinetic_energy
  !! Verifies the one-pass global squared norm used for kinetic energy.

  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, VERT
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
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  class(field_t), pointer :: u, v, w

  integer :: dims_global(3), nproc_dir(3)
  integer :: ierr, nrank, nproc
  real(dp) :: energy, norm_squared, expected_norm_squared
  real(dp), parameter :: expected_energy = 7._dp
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  real(dp), parameter :: tolerance = 100._dp*epsilon(1._dp)
  character(len=8), parameter :: periodic(2) = &
    ['periodic', 'periodic']
  logical :: test_pass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! ny is not a multiple of SZ, so each field contains padding that the
  ! reduction must exclude.  Keep twenty physical z points on every rank.
  dims_global = [17, 18, 20*nproc]
  nproc_dir = [1, 1, nproc]
  mesh = mesh_t(dims_global, nproc_dir, lengths, &
                periodic, periodic, periodic)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  call u%fill(1._dp)
  call v%fill(2._dp)
  call w%fill(-3._dp)

  norm_squared = backend%vector_norm_squared(u, v, w)
  expected_norm_squared = 14._dp*real(product(dims_global), dp)
  energy = 0.5_dp*norm_squared/real(product(dims_global), dp)
  test_pass = abs(norm_squared - expected_norm_squared) &
              <= tolerance*expected_norm_squared
  test_pass = test_pass .and. abs(energy - expected_energy) <= tolerance
  call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, MPI_LOGICAL, MPI_LAND, &
                     MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    print *, 'kinetic energy:', energy
    print *, 'expected:', expected_energy
    if (test_pass) then
      print *, 'PASS'
    else
      print *, 'FAIL'
    end if
  end if

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)
  call MPI_Finalize(ierr)

  if (.not. test_pass) error stop 'Kinetic-energy test failed.'

end program test_kinetic_energy
