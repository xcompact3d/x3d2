program xcompact
  use mpi

  use m_allocator
  use m_base_backend
  use m_common, only: pi, globs_t, POISSON_SOLVER_FFT, POISSON_SOLVER_CG, &
                      DIR_X, DIR_Y, DIR_Z
  use m_domain, only: domain_decomposition
  use m_solver, only: solver_t
  use m_time_integrator, only: time_intg_t
  use m_tdsops, only: tdsops_t

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
  type(solver_t) :: solver
  type(time_intg_t) :: time_integrator
  type(dirps_t) :: xdirps, ydirps, zdirps

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
#endif

  real(dp), allocatable, dimension(:, :, :) :: u, v, w

  real(dp) :: t_start, t_end
  integer :: dims(3)
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

  ! read L_x/y/z from the input file
  globs%Lx = 2*pi; globs%Ly = 2*pi; globs%Lz = 2*pi
  xdirps%L = globs%Lx; ydirps%L = globs%Ly; zdirps%L = globs%Lz

  ! read ns from the input file
  globs%nx = 256; globs%ny = 256; globs%nz = 256

  globs%dt = 0.001_dp
  globs%nu = 1._dp/1600._dp
  globs%n_iters = 20000
  globs%n_output = 100

  ! set nprocs based on run time arguments
  globs%nproc_x = 1; globs%nproc_y = 1; globs%nproc_z = 1

  globs%poisson_solver_type = POISSON_SOLVER_FFT

  xdirps%nproc = globs%nproc_x
  ydirps%nproc = globs%nproc_y
  zdirps%nproc = globs%nproc_z

  ! lets assume simple cases for now
  globs%nx_loc = globs%nx/globs%nproc_x
  globs%ny_loc = globs%ny/globs%nproc_y
  globs%nz_loc = globs%nz/globs%nproc_z

  globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
  globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
  globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

  globs%dx = globs%Lx/globs%nx
  globs%dy = globs%Ly/globs%ny
  globs%dz = globs%Lz/globs%nz

  xdirps%d = globs%dx; ydirps%d = globs%dy; zdirps%d = globs%dz

  xdirps%n = globs%nx_loc
  ydirps%n = globs%ny_loc
  zdirps%n = globs%nz_loc

  xdirps%n_blocks = globs%n_groups_x
  ydirps%n_blocks = globs%n_groups_y
  zdirps%n_blocks = globs%n_groups_z

  xdirps%dir = DIR_X; ydirps%dir = DIR_Y; zdirps%dir = DIR_Z

  call domain_decomposition(xdirps, ydirps, zdirps, nrank, nproc)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(globs%nx_loc, globs%ny_loc, &
                                    globs%nz_loc, SZ)
  allocator => cuda_allocator
  if (nrank == 0) print *, 'CUDA allocator instantiated'

  cuda_backend = cuda_backend_t(globs, allocator)
  backend => cuda_backend
  if (nrank == 0) print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(globs%nx_loc, globs%ny_loc, globs%nz_loc, SZ)
  allocator => omp_allocator
  if (nrank == 0) print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(globs, allocator)
  backend => omp_backend
  if (nrank == 0) print *, 'OpenMP backend instantiated'
#endif

  dims(:) = allocator%cdims_padded
  allocate (u(dims(1), dims(2), dims(3)))
  allocate (v(dims(1), dims(2), dims(3)))
  allocate (w(dims(1), dims(2), dims(3)))

  time_integrator = time_intg_t(allocator=allocator, backend=backend)
  if (nrank == 0) print *, 'time integrator instantiated'
  solver = solver_t(backend, time_integrator, xdirps, ydirps, zdirps, globs)
  if (nrank == 0) print *, 'solver instantiated'

  call cpu_time(t_start)

  call solver%run(u, v, w)

  call cpu_time(t_end)

  if (nrank == 0) print *, 'Time: ', t_end - t_start

  if (nrank == 0) print *, 'norms', norm2(u), norm2(v), norm2(w)

  call MPI_Finalize(ierr)

end program xcompact
