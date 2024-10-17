program xcompact
  use mpi

  use m_allocator
  use m_base_backend
  use m_common, only: pi, globs_t, POISSON_SOLVER_FFT, POISSON_SOLVER_CG, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C
  use m_solver, only: solver_t
  use m_time_integrator, only: time_intg_t
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

  type(globs_t) :: globs
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), pointer :: host_allocator
  type(solver_t) :: solver
  type(time_intg_t) :: time_integrator
  type(dirps_t) :: xdirps, ydirps, zdirps
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
  integer, dimension(3) :: dims_global
  integer, dimension(3) :: nproc_dir
  real(dp), dimension(3) :: L_global
  integer :: nrank, nproc, ierr

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
#endif

  ! Global number of cells in each direction
  dims_global = [96, 96, 96]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, nproc]

  mesh = mesh_t(dims_global, nproc_dir, L_global)

  globs%dt = 0.001_dp
  globs%nu = 1._dp/1600._dp
  globs%n_iters = 20000
  globs%n_output = 10

  globs%poisson_solver_type = POISSON_SOLVER_FFT

  xdirps%dir = DIR_X; ydirps%dir = DIR_Y; zdirps%dir = DIR_Z

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

  time_integrator = time_intg_t(allocator=allocator, backend=backend)
  if (nrank == 0) print *, 'time integrator instantiated'
  solver = solver_t(backend, mesh, time_integrator, host_allocator, &
                    xdirps, ydirps, zdirps, globs)
  if (nrank == 0) print *, 'solver instantiated'

  call cpu_time(t_start)

  call solver%run()

  call cpu_time(t_end)

  if (nrank == 0) print *, 'Time: ', t_end - t_start

  call MPI_Finalize(ierr)

end program xcompact
