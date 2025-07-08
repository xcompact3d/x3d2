program test_setget_field

  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_C, DIR_X, DIR_Y, DIR_Z, VERT
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

  class(field_t), pointer :: fld, fld_c
  real(dp), dimension(:, :, :), allocatable :: arr
  integer, dimension(3) :: shape_c

  integer :: ierr

  call MPI_Init(ierr)

  mesh = mesh_t([16, 32, 48], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"])

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  fld => backend%allocator%get_block(DIR_X, VERT)
  fld_c => backend%allocator%get_block(DIR_C, VERT)
  shape_c = fld_c%get_shape()
  allocate (arr(shape_c(1), shape_c(2), shape_c(3)))
  arr = 1.0_dp
  call backend%set_field_data(fld, arr)

  if (fld%data_loc /= VERT) then
    error stop "Field location was changed by set_field_data"
  end if

  deallocate (arr)
  call backend%allocator%release_block(fld)
  call backend%allocator%release_block(fld_c)

  call MPI_Finalize(ierr)

end program test_setget_field
