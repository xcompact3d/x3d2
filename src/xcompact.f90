program xcompact
  use mpi

  use m_allocator
  use m_base_backend
  use m_base_case, only: base_case_t
  use m_common, only: pi
  use m_solver, only: read_solver_input
  use m_mesh
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

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  class(base_case_t), allocatable :: flow_case

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
  character(len=20) :: flow_case_name
  integer, dimension(3) :: dims_global
  integer, dimension(3) :: nproc_dir = 0
  real(dp), dimension(3) :: L_global
  character(3) :: poisson_solver_type
  character(32) :: backend_name
  integer :: nrank, nproc, ierr
  logical :: use_2decomp

  namelist /domain_settings/ flow_case_name, L_global, dims_global, &
    nproc_dir, BC_x, BC_y, BC_z

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
    read (100, nml=domain_settings)
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

  call read_solver_input(i_poisson_solver_type=poisson_solver_type)

  ! Decide whether 2decomp is used or not
  use_2decomp = poisson_solver_type == 'FFT' .and. trim(backend_name) == 'OMP'

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                use_2decomp=use_2decomp)

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

  if (nrank == 0) print *, 'Flow case: ', flow_case_name

  select case (trim(flow_case_name))
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
