program test_setget_field

  use mpi 

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, VERT
#ifdef CUDA
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif
  use m_mesh, only: mesh_t
  
  implicit none
  
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  type(mesh_t) :: mesh
  
  class(field_t), pointer :: fld
  real(dp), dimension(:, :, :), allocatable :: arr

  integer :: nx, ny, nz

  integer :: ierr
  
  call MPI_Init(ierr)

  mesh = mesh_t ([16, 32, 48], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
    ["periodic", "periodic"], &
    ["periodic", "periodic"], &
    ["periodic", "periodic"])

  nx = mesh%get_n(DIR_X, VERT)
  ny = mesh%get_n(DIR_Y, VERT)
  nz = mesh%get_n(DIR_Z, VERT)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh, SZ)
  allocator => cuda_allocator
#else
  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif
  
  fld => backend%allocator%get_block(DIR_X, VERT)
  arr = 1.0_dp
  call backend%set_field_data(fld, arr)
  
  if (fld%data_loc /= VERT) then
    error stop "Field location was changed by set_field_data"
  end if

  call MPI_Finalize(ierr)
  
end program test_setget_field
