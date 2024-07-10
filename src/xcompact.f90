program xcompact
  use mpi

  use m_allocator
  use m_base_backend
  use m_common, only: pi
  use m_solver, only: solver_t
  use m_tdsops, only: tdsops_t
  use m_mesh

#ifdef CUDA
  use m_cuda_allocator
  use m_cuda_backend
  use m_cuda_common, only: SZ
  use m_cuda_tdsops, only: cuda_tdsops_t
#else
  use m_omp_backend
  use m_omp_common, only: SZ
#endif

  implicit none

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(mesh_t) :: mesh
  type(allocator_t), pointer :: host_allocator
  type(solver_t) :: solver

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
#endif

  type(allocator_t), target :: omp_allocator

  real(dp) :: t_start, t_end

  character(len=200) :: input_file
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  integer, dimension(3) :: dims_global
  integer, dimension(3) :: nproc_dir = 0
  real(dp), dimension(3) :: L_global
  integer :: nrank, nproc, ierr

  namelist /domain_params/ L_global, dims_global, nproc_dir, BC_x, BC_y, BC_z

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
#endif

  if (command_argument_count() >= 1) then
    call get_command_argument(1, input_file)
    open (100, file=input_file)
    read (100, nml=domain_params)
    close (100)
  else
    error stop 'Input file is not provided.'
  end if

  if (product(nproc_dir) /= nproc) then
    if (nrank == 0) print *, 'nproc_dir specified in the input file does &
                              &not match the total number of ranks, falling &
                              &back to a 1D decomposition along Z-dir instead.'
    nproc_dir = [1, 1, nproc]
  end if

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh, SZ)
  allocator => cuda_allocator
  if (nrank == 0) print *, 'CUDA allocator instantiated'

  omp_allocator = allocator_t(mesh, SZ)
  host_allocator => omp_allocator

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  if (nrank == 0) print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator
  if (nrank == 0) print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  if (nrank == 0) print *, 'OpenMP backend instantiated'
#endif

  solver = solver_t(backend, mesh, host_allocator)
  if (nrank == 0) print *, 'solver instantiated'

  call cpu_time(t_start)

  call solver%run()

  call cpu_time(t_end)

  if (nrank == 0) print *, 'Time: ', t_end - t_start

  call MPI_Finalize(ierr)

end program xcompact
