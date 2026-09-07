program test_setget_field

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t
  use m_common, only: dp, DIR_C, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  type(backend_runtime_t), target :: runtime
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  type(mesh_t), target :: mesh

  class(field_t), pointer :: fld, fld_c
  real(dp), dimension(:, :, :), allocatable :: arr
  integer, dimension(3) :: shape_c

  integer :: nrank, nproc
  logical :: allpass

  call initialise_mpi(nrank, nproc)

  print *, "Initialised MPI"

  mesh = mesh_t([16, 32, 48], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"])

  print *, "Initialised mesh"

  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend

  print *, "Initialised backend"

  fld => backend%allocator%get_block(DIR_X, VERT)
  fld_c => backend%allocator%get_block(DIR_C, VERT)
  shape_c = fld_c%get_shape()
  print *, shape_c
  allocate (arr(shape_c(1), shape_c(2), shape_c(3)))
  arr = 1.0_dp

  print *, "Initialised data"

  call backend%set_field_data(fld, arr)

  print *, "Set field data"

  allpass = .true.
  if (fld%data_loc /= VERT) then
    print *, "Field location was changed by set_field_data"
    allpass = .false.
  end if

  arr = 0.0_dp
  call backend%get_field_data(arr, fld)
  if (any(arr /= 1.0_dp)) then
    print *, "Getting/setting field data failed"
    allpass = .false.
  end if

  print *, "Get field data"

  deallocate (arr)
  call backend%allocator%release_block(fld)
  call backend%allocator%release_block(fld_c)

  call finalise_test(allpass, nrank)

end program test_setget_field
