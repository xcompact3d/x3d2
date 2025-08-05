program test_fft
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_tdsops, only: dirps_t
  use m_solver, only: allocate_tdsops

  use m_common, only: dp, pi, MPI_X3D2_DP, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, CELL

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

  class(field_t), pointer :: input_field, output_field

  integer :: dims(3)

  integer :: nrank, nproc
  integer :: ierr, i, j, k

  real(dp) :: dx, dx_per

  class(base_backend_t), pointer :: backend
  class(mesh_t), allocatable :: mesh
  class(allocator_t), pointer :: allocator
  type(dirps_t), pointer :: xdirps, ydirps, zdirps
  integer, dimension(3) :: dims_padded, dims_global, nproc_dir
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  real(dp) :: x, y, z
  real(dp) :: error_norm
  real(dp), dimension(3) :: xloc
  real(dp), parameter :: tol = 1e-10
  logical :: use_2decomp
  real(dp), allocatable, dimension(:, :, :) :: input_data, output_data

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
  use_2decomp = .false.
#else
  use_2decomp = .true.
#endif

  ! Global number of cells in each direction
  dims_global = [64, 32, 128]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, nproc]

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, &
                BC_x, BC_y, BC_z, &
                use_2decomp=use_2decomp)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
  print *, 'CUDA allocator instantiated'

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  allocate (xdirps, ydirps, zdirps)
  xdirps%dir = DIR_X
  ydirps%dir = DIR_Y
  zdirps%dir = DIR_Z
  call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')
  call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')
  call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')

  input_field => allocator%get_block(DIR_C, CELL)
  output_field => allocator%get_block(DIR_C, CELL)

  call input_field%fill(0._dp)
  call output_field%fill(0._dp)

  dims = mesh%get_dims(CELL)
  allocate (input_data(dims(1), dims(2), dims(3)))

  ! Initialise field with some function
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        xloc = mesh%get_coordinates(i, j, k)
        x = xloc(1)
        y = xloc(2)
        z = xloc(3)
        input_data(i, j, k) = sin(x)*cos(y)*cos(z) + 2*x
      end do
    end do
  end do

  call backend%set_field_data(input_field, input_data, DIR_C)

  call backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)

  ! Compute FFT and back
  call backend%poisson_fft%fft_forward(input_field)
  call backend%poisson_fft%fft_backward(output_field)

  allocate (output_data(dims(1), dims(2), dims(3)))
  call backend%get_field_data(output_data, output_field, DIR_C)
  ! The output scaled with number of cells in domain, hence the first '/product(dims_global)'.
  ! RMS value is used for the norm, hence the second '/product(dims_global)'
  error_norm = norm2( &
               input_data - output_data/product(dims_global) &
               )**2/product(dims_global)
  call MPI_Allreduce(MPI_IN_PLACE, error_norm, 1, MPI_X3D2_DP, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)
  error_norm = sqrt(error_norm)

  if (error_norm > tol) then
    if (mesh%par%is_root()) then
      print *, "error in FFT result, error norm=", error_norm
    end if
    error stop 'TEST FAILED.'
  else
    if (mesh%par%is_root()) then
      print *, "TEST PASS"
    end if
  end if

  call MPI_Finalize(ierr)

end program test_fft
