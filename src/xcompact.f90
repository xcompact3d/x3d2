program xcompact
  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size, &
                   MPI_Finalize, MPI_Init

  use m_allocator
  use m_base_backend
  use m_base_case, only: base_case_t
  use m_common, only: pi, get_argument, VERT
  use m_config, only: domain_config_t, solver_config_t
  use m_logger, only: log_info, log_warn, abort_run, str
  use m_mesh
  use m_case_abl, only: case_abl_t
  use m_case_channel, only: case_channel_t
  use m_case_cylinder, only: case_cylinder_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t
  use m_case_wind_turbine, only: case_wind_turbine_t

#ifdef CUDA
  use m_cuda_allocator
  use m_cuda_backend
  use m_cuda_common, only: SZ
  use m_cuda_errors, only: check_cuda
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

  type(domain_config_t) :: domain_cfg
  type(solver_config_t) :: solver_cfg
  character(32) :: backend_name
  integer :: dims(3), nrank, nproc, ierr
  logical :: use_2decomp

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  call log_info('Parallel run with '//str(nproc)//' ranks')
  call log_info('Data precision is '//str(dp))

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  call check_cuda(ierr, 'cudaGetDeviceCount')
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  call check_cuda(ierr, 'cudaSetDevice')
  ierr = cudaGetDevice(devnum)
  call check_cuda(ierr, 'cudaGetDevice')
  backend_name = "CUDA"
#else
  backend_name = "OMP"
#endif

  call domain_cfg%read(nml_file=get_argument(1))
  call solver_cfg%read(nml_file=get_argument(1))

  if (product(domain_cfg%nproc_dir) /= nproc) then
    call log_warn('nproc_dir specified in the input file does not match the &
                  &total number of ranks, falling back to a 1D decomposition &
                  &along Z-dir instead.')
    domain_cfg%nproc_dir = [1, 1, nproc]
  end if

  ! Decide whether 2decomp is used or not
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
  call log_info('CUDA allocator instantiated')

  omp_allocator = allocator_t(dims, SZ)
  host_allocator => omp_allocator

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  call log_info('CUDA backend instantiated')
#else
  omp_allocator = allocator_t(dims, SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator
  call log_info('OpenMP allocator instantiated')

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  call log_info('OpenMP backend instantiated')
#endif

  call log_info('Flow case: '//trim(domain_cfg%flow_case_name))

  select case (trim(domain_cfg%flow_case_name))
  case ('abl')
    allocate (case_abl_t :: flow_case)
    flow_case = case_abl_t(backend, mesh, host_allocator)
  case ('channel')
    allocate (case_channel_t :: flow_case)
    flow_case = case_channel_t(backend, mesh, host_allocator)
  case ('cylinder')
    allocate (case_cylinder_t :: flow_case)
    flow_case = case_cylinder_t(backend, mesh, host_allocator)
  case ('generic')
    allocate (case_generic_t :: flow_case)
    flow_case = case_generic_t(backend, mesh, host_allocator)
  case ('tgv')
    allocate (case_tgv_t :: flow_case)
    flow_case = case_tgv_t(backend, mesh, host_allocator)
  case ('wind_turbine')
    allocate (case_wind_turbine_t :: flow_case)
    flow_case = case_wind_turbine_t(backend, mesh, host_allocator)
  case default
    call abort_run('Undefined flow_case: ' &
                   //trim(domain_cfg%flow_case_name), collective=.true.)
  end select
  call log_info('solver instantiated')

  call flow_case%run()

  call MPI_Finalize(ierr)

end program xcompact
