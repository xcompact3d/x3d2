module m_backend_setup
  !! Creates the backend and allocators for the backend this build targets,
  !! so that the solver driver (xcompact.f90) needs no backend-specific code.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_rank

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: VERT
  use m_logger, only: log_info, abort_run
  use m_mesh, only: mesh_t

#ifdef CUDA
  use cudafor
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
  use m_cuda_errors, only: check_cuda
#elif defined(OMP_TGT)
  use omp_lib, only: omp_get_num_devices, omp_set_default_device
  use m_omptgt_allocator, only: omptgt_allocator_t
  use m_omptgt_backend, only: omptgt_backend_t
  use m_omptgt_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  private
  public :: backend_setup_t, backend_name

#ifdef CUDA
  character(len=*), parameter :: backend_name = 'CUDA'
#elif defined(OMP_TGT)
  character(len=*), parameter :: backend_name = 'OMP_TGT'
#else
  character(len=*), parameter :: backend_name = 'OMP'
#endif

  type :: backend_setup_t
    !! Owns the backend and allocators.  The pointer components alias
    !! concrete members of this same object, so instances must be declared
    !! with the TARGET attribute and must never be copied or assigned.
    class(base_backend_t), pointer :: backend => null()
    class(allocator_t), pointer :: allocator => null()
    type(allocator_t), pointer :: host_allocator => null()
#ifdef CUDA
    type(cuda_backend_t) :: cuda_backend
    type(cuda_allocator_t) :: cuda_allocator
#elif defined(OMP_TGT)
    type(omptgt_backend_t) :: omptgt_backend
    type(omptgt_allocator_t) :: omptgt_allocator
#else
    type(omp_backend_t) :: omp_backend
#endif
    type(allocator_t) :: omp_allocator
  contains
    procedure :: init
  end type backend_setup_t

contains

  subroutine init(self, mesh)
    class(backend_setup_t), target, intent(inout) :: self
    type(mesh_t), target, intent(inout) :: mesh

    integer :: dims(3), nrank, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call select_device(nrank)

    ! get local vertex dimensions
    dims = mesh%get_dims(VERT)
#ifdef CUDA
    self%cuda_allocator = cuda_allocator_t(dims, SZ)
    self%allocator => self%cuda_allocator
    call log_info('CUDA allocator instantiated')

    self%omp_allocator = allocator_t(dims, SZ)
    self%host_allocator => self%omp_allocator

    self%cuda_backend = cuda_backend_t(mesh, self%allocator)
    self%backend => self%cuda_backend
    call log_info('CUDA backend instantiated')
#elif defined(OMP_TGT)
    self%omptgt_allocator = omptgt_allocator_t(dims, SZ)
    self%allocator => self%omptgt_allocator
    call log_info('OpenMP target allocator instantiated')

    self%omp_allocator = allocator_t(dims, SZ)
    self%host_allocator => self%omp_allocator

    self%omptgt_backend = omptgt_backend_t(mesh, self%allocator)
    self%backend => self%omptgt_backend
    call log_info('OpenMP target backend instantiated')
#else
    self%omp_allocator = allocator_t(dims, SZ)
    self%allocator => self%omp_allocator
    self%host_allocator => self%omp_allocator
    call log_info('OpenMP allocator instantiated')

    self%omp_backend = omp_backend_t(mesh, self%allocator)
    self%backend => self%omp_backend
    call log_info('OpenMP backend instantiated')
#endif
  end subroutine init

  subroutine select_device(nrank)
    !! Select a device round-robin by MPI rank; a no-op on the CPU backend.
    integer, intent(in) :: nrank

#ifdef CUDA
    integer :: ierr, ndevs

    ierr = cudaGetDeviceCount(ndevs)
    call check_cuda(ierr, 'cudaGetDeviceCount')
    if (ndevs < 1) call abort_run('select_device: no CUDA devices available')
    ierr = cudaSetDevice(mod(nrank, ndevs))
    call check_cuda(ierr, 'cudaSetDevice')
#elif defined(OMP_TGT)
    integer :: ndevs

    ndevs = omp_get_num_devices()
    if (ndevs < 1) call abort_run('select_device: no OpenMP target devices &
                                  &available')
    call omp_set_default_device(mod(nrank, ndevs))
#endif
  end subroutine select_device

end module m_backend_setup
