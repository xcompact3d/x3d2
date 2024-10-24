program test_reorder
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_tdsops, only: dirps_t

  use m_common, only: dp, pi, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT

  use m_ordering, only: get_index_dir, get_index_ijk
  use m_mesh, only: mesh_t

#ifdef CUDA
  use cudafor

  use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_common, only: SZ
  use m_omp_backend, only: omp_backend_t
#endif

  implicit none

  logical :: allpass = .true.
  class(field_t), pointer :: u_x, u_y, u_z
  class(field_t), pointer :: u_x_original

  real(dp), allocatable, dimension(:, :, :) :: u_array, temp_1, temp_2

  integer :: dims(3)

  integer :: nrank, nproc
  integer :: ierr, i, j, k

  real(dp) :: dx, dx_per

  class(base_backend_t), pointer :: backend
  class(mesh_t), allocatable :: mesh
  class(allocator_t), pointer :: allocator
  type(dirps_t), target :: xdirps, ydirps, zdirps
  integer, dimension(3) :: dims_padded, dims_global, nproc_dir
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  logical :: pass_X, pass_Y, pass_Z

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
#endif

  ! Initialise variables and arrays
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
#endif

  ! Global number of cells in each direction
  dims_global = [32, 64, 96]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, nproc]

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh, SZ)
  allocator => cuda_allocator
  print *, 'CUDA allocator instantiated'

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator
  print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  print *, 'OpenMP backend instantiated'
#endif

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'
  pass_X = .true.
  pass_Y = .true.
  pass_Z = .true.

  dims_padded = mesh%get_padded_dims(DIR_C)

  ! Test indexing only
  do k = 1, mesh%get_n(DIR_Z, VERT)
    do j = 1, mesh%get_n(DIR_Y, VERT)
      do i = 1, mesh%get_n(DIR_X, VERT)
        call test_index_reversing(pass_X, i, j, k, DIR_X, &
                                  mesh%get_sz(), dims_padded(1), &
                                  dims_padded(2), dims_padded(3))
        call test_index_reversing(pass_Y, i, j, k, DIR_Y, &
                                  mesh%get_sz(), dims_padded(1), &
                                  dims_padded(2), dims_padded(3))
        call test_index_reversing(pass_Z, i, j, k, DIR_Z, &
                                  mesh%get_sz(), dims_padded(1), &
                                  dims_padded(2), dims_padded(3))
      end do
    end do
  end do
  if (.not. pass_X) print *, "Error in X direction for index reversing"
  if (.not. pass_Y) print *, "Error in Y direction for index reversing"
  if (.not. pass_Z) print *, "Error in Z direction for index reversing"

  allpass = (pass_X .and. pass_Y .and. pass_Z)

  ! Test reordering
  u_x => allocator%get_block(DIR_X)
  u_y => allocator%get_block(DIR_Y)
  u_z => allocator%get_block(DIR_Z)
  u_x_original => allocator%get_block(DIR_X)

  dims(:) = mesh%get_padded_dims(DIR_X)
  allocate (u_array(dims(1), dims(2), dims(3)))

  call random_number(u_array)

#ifdef CUDA
  allocate (temp_1(dims(1), dims(2), dims(3)))
  allocate (temp_2(dims(1), dims(2), dims(3)))

  select type (u_x_original)
  type is (cuda_field_t)
    u_x_original%data_d = u_array
  end select
#else
  select type (u_x_original)
  type is (field_t)
    u_x_original%data = u_array
  end select
#endif

  call backend%reorder(u_y, u_x_original, RDR_X2Y)
  call backend%reorder(u_x, u_y, RDR_Y2X)
  call check_reorder(allpass, u_x, u_x_original, "testing X2Y and Y2X failed")

  call backend%reorder(u_z, u_x, RDR_X2Z)
  call backend%reorder(u_x, u_z, RDR_Z2X)
  call check_reorder(allpass, u_x, u_x_original, "testing X2Z and Z2X failed")

  call backend%reorder(u_z, u_y, RDR_Y2Z)
  call backend%reorder(u_x, u_z, RDR_Z2X)
  call check_reorder(allpass, u_x, u_x_original, "testing Y2Z and Z2X failed")

  call backend%reorder(u_y, u_z, RDR_Z2Y)
  call backend%reorder(u_x, u_y, RDR_Y2X)
  call check_reorder(allpass, u_x, u_x_original, "testing Z2Y and Y2X failed")

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine test_index_reversing(pass, i, j, k, dir, SZ, nx, ny, nz)
    logical, intent(inout) :: pass
    integer, intent(in) :: i, j, k    ! original indices in the cartesian space
    integer, intent(in) :: dir
    integer, intent(in) :: SZ, nx, ny, nz
    integer :: dir_i, dir_j, dir_k    ! indices in the applicatin storage direction
    integer :: cart_i, cart_j, cart_k ! newly computed indices in the cartesian space

    call get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, SZ, nx, ny, nz)
    call get_index_ijk(cart_i, cart_j, cart_k, dir_i, dir_j, dir_k, dir, &
                       SZ, nx, ny, nz)

    if (i /= cart_i .or. j /= cart_j .or. k /= cart_k) then
      pass = .false.
    end if

  end subroutine

  subroutine check_reorder(allpass, a, b, message)
    logical, intent(inout) :: allpass
    class(field_t), intent(in) :: a, b
    character(len=*), intent(in) :: message
    real(dp) :: tol = 1d-8

#ifdef CUDA
    select type (a); type is (cuda_field_t); temp_1 = a%data_d; end select
    select type (b); type is (cuda_field_t); temp_2 = b%data_d; end select
    if (norm2(temp_1 - temp_2) > tol) then
      allpass = .false.
      write (stderr, '(a)') message
    end if
#else
    if (norm2(a%data(:, 1:mesh%get_n(DIR_X, VERT), :) - &
              b%data(:, 1:mesh%get_n(DIR_X, VERT), :)) > tol) then
      allpass = .false.
      write (stderr, '(a)') message
    end if
#endif

  end subroutine

end program test_reorder

