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
  type(allocator_t), pointer :: host_allocator
  type(solver_t) :: solver
  type(mesh_t), target :: mesh

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
  character(3) :: poisson_solver_type
  character(32) :: backend_name
  integer :: nrank, nproc, ierr
  logical :: use_2decomp

  namelist /domain_params/ L_global, dims_global, nproc_dir, BC_x, BC_y, BC_z
  namelist /solver_params/ poisson_solver_type

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
  backend_name = "CUDA"
#else 
  backend_name = "OMP"
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

  ! Decide whether 2decomp is used or not
  use_2decomp = (poisson_solver_type == 'FFT' .and.  trim(backend_name) == 'OMP')

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, use_2decomp=use_2decomp)

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
