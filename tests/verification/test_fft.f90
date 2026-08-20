program test_fft
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda
  use m_tdsops, only: dirps_t
  use m_solver, only: allocate_tdsops

  use m_common, only: dp, pi, MPI_X3D2_DP, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, CELL

  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  class(field_t), pointer :: input_field, output_field

  integer :: dims(3)

  integer :: nrank, nproc
  integer :: ierr, i, j, k

  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend
  type(mesh_t), target :: mesh
  class(allocator_t), pointer :: allocator
  type(dirps_t), pointer :: xdirps, ydirps, zdirps
  integer, dimension(3) :: dims_padded, dims_global, nproc_dir
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  real(dp) :: x, y, z
  real(dp) :: error_norm
  real(dp), dimension(3) :: xloc
#ifdef SINGLE_PREC
  real(dp), parameter :: tol = merge(1e-6_dp, 1e-5_dp, backend_is_cuda)
#else
  real(dp), parameter :: tol = 1e-10_dp
#endif
  logical :: use_2decomp
  logical :: allpass
  real(dp), allocatable, dimension(:, :, :) :: input_data, output_data

  ! Initialise variables and arrays
  call initialise_mpi(nrank, nproc)

  use_2decomp = .not. backend_is_cuda

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

  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend

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

  allpass = (error_norm <= tol)
  if (.not. allpass) then
    if (mesh%par%is_root()) then
      print *, "error in FFT result, error norm=", error_norm
    end if
  else
    if (mesh%par%is_root()) then
      print *, "TEST PASS"
    end if
  end if

  call allocator%release_block(input_field)
  call allocator%release_block(output_field)
  call finalise_test(allpass, nrank)

end program test_fft
