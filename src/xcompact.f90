program xcompact
  !! Main program for X3D2 CFD solver.
  !!
  !! X3D2 is a high-order finite-difference incompressible Navier-Stokes
  !! solver based on Xcompact3D/Incompact3D. It solves the incompressible
  !! Navier-Stokes equations using:
  !!
  !! - **Compact finite differences** for spatial derivatives (4th-6th order)
  !! - **Fractional-step method** for pressure-velocity coupling
  !! - **FFT-based or iterative Poisson solvers** for pressure
  !! - **Explicit time integration** (Runge-Kutta or Adams-Bashforth)
  !!
  !! **Program Flow:**
  !!
  !! 1. Initialise MPI and determine rank/size
  !! 2. Select computational backend (CUDA GPU or OpenMP CPU)
  !! 3. Read configuration from input file (domain and solver parameters)
  !! 4. Create mesh with domain decomposition (pencil decomposition)
  !! 5. Instantiate allocator and backend for the selected platform
  !! 6. Select and instantiate flow case (channel, TGV, generic, etc.)
  !! 7. Run simulation via flow_case%run()
  !! 8. Report timing and finalise MPI
  !!
  !! **Backend Options:**
  !!
  !! - **CUDA**: GPU acceleration via NVIDIA CUDA (compile with -DCUDA)
  !! - **OMP**: CPU parallelism via OpenMP threading
  !!
  !! **Input:** Namelist file specified as command-line argument (e.g., input.x3d)
  !!
  !! **Domain Decomposition:**
  !!
  !! X3D2 supports two decomposition strategies:
  !!
  !! - **2DECOMP&FFT**: External library used when FFT Poisson solver + OMP backend.
  !!   Provides optimised pencil decomposition and FFT transforms. Cannot decompose
  !!   in X-direction (`nproc_dir(1)` must be 1).
  !! - **Generic**: Built-in X3D2 decomposition used for CUDA backend or when
  !!   2DECOMP&FFT is unavailable. Can decompose in any direction (X, Y, Z).
  !!
  !! The decomposition is selected automatically based on backend and solver type.
  use mpi

  use m_allocator
  use m_base_backend
  use m_base_case, only: base_case_t
  use m_common, only: pi, get_argument, VERT
  use m_config, only: domain_config_t, solver_config_t
  use m_mesh
  use m_case_channel, only: case_channel_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t

#ifdef CUDA
  use m_cuda_allocator
  use m_cuda_backend
  use m_cuda_common, only: SZ
#else
  use m_omp_backend
  use m_omp_common, only: SZ
#endif

  implicit none

  class(base_backend_t), pointer :: backend       !! Active computational backend (CUDA or OMP)
  class(allocator_t), pointer :: allocator        !! Memory allocator for device/host
  type(allocator_t), pointer :: host_allocator    !! Host memory allocator (for I/O, etc.)
  type(mesh_t), target :: mesh                    !! Computational mesh with decomposition
  class(base_case_t), allocatable :: flow_case    !! Flow case instance (polymorphic)

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend    !! CUDA backend implementation
  type(cuda_allocator_t), target :: cuda_allocator !! CUDA device memory allocator
  integer :: ndevs, devnum                         !! Number of GPUs, assigned device number
#else
  type(omp_backend_t), target :: omp_backend       !! OpenMP backend implementation
#endif

  type(allocator_t), target :: omp_allocator       !! Host/CPU memory allocator

  real(dp) :: t_start, t_end                       !! CPU timing for performance measurement

  type(domain_config_t) :: domain_cfg              !! Domain configuration from input file
  type(solver_config_t) :: solver_cfg              !! Solver configuration from input file
  character(32) :: backend_name                    !! Backend name string ("CUDA" or "OMP")
  integer :: dims(3), nrank, nproc, ierr           !! Dimensions, MPI rank/size, error code
  logical :: use_2decomp                           !! Whether to use 2DECOMP&FFT library

  ! Initialise MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) then
    print *, 'Parallel run with', nproc, 'ranks'
    print *, 'Data precision is', dp
  end if

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
  backend_name = "CUDA"
#else
  backend_name = "OMP"
#endif

  call domain_cfg%read(nml_file=get_argument(1))
  call solver_cfg%read(nml_file=get_argument(1))

  if (product(domain_cfg%nproc_dir) /= nproc) then
    if (nrank == 0) print *, 'nproc_dir specified in the input file does &
                              &not match the total number of ranks, falling &
                              &back to a 1D decomposition along Z-dir instead.'
    domain_cfg%nproc_dir = [1, 1, nproc]
  end if

  ! Select decomposition strategy:
  ! - 2DECOMP&FFT: Used for FFT Poisson solver with OMP backend (optimised)
  ! - Generic: Used for CUDA backend or non-FFT solvers (more flexible)
  use_2decomp = solver_cfg%poisson_solver_type == 'FFT' &
                .and. trim(backend_name) == 'OMP'

  mesh = mesh_t(domain_cfg%dims_global, domain_cfg%nproc_dir, &
                domain_cfg%L_global, domain_cfg%BC_x, domain_cfg%BC_y, &
                domain_cfg%BC_z, domain_cfg%stretching, domain_cfg%beta, &
                use_2decomp=use_2decomp)

  ! get local vertex dimensions
  dims = mesh%get_dims(VERT)
#ifdef CUDA
  cuda_allocator = cuda_allocator_t(dims, SZ)
  allocator => cuda_allocator
  if (nrank == 0) print *, 'CUDA allocator instantiated'

  omp_allocator = allocator_t(dims, SZ)
  host_allocator => omp_allocator

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  if (nrank == 0) print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(dims, SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator
  if (nrank == 0) print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  if (nrank == 0) print *, 'OpenMP backend instantiated'
#endif

  if (nrank == 0) print *, 'Flow case: ', domain_cfg%flow_case_name

  select case (trim(domain_cfg%flow_case_name))
  case ('channel')
    allocate (case_channel_t :: flow_case)
    flow_case = case_channel_t(backend, mesh, host_allocator)
  case ('generic')
    allocate (case_generic_t :: flow_case)
    flow_case = case_generic_t(backend, mesh, host_allocator)
  case ('tgv')
    allocate (case_tgv_t :: flow_case)
    flow_case = case_tgv_t(backend, mesh, host_allocator)
  case default
    error stop 'Undefined flow_case.'
  end select
  if (nrank == 0) print *, 'solver instantiated'

  call cpu_time(t_start)

  call flow_case%run()

  call cpu_time(t_end)

  if (nrank == 0) print *, 'Time: ', t_end - t_start

  call MPI_Finalize(ierr)

end program xcompact
