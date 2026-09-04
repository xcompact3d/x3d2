module m_backend_runtime
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: VERT
  use m_mesh, only: mesh_t

#ifdef CUDA
  use cudafor
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#elif defined(OMP_TGT)
  use m_omp_common, only: SZ
  use m_omptgt_allocator, only: omptgt_allocator_t
  use m_omptgt_backend, only: omptgt_backend_t
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  integer, parameter :: backend_sz = SZ

#ifdef CUDA
  logical, parameter :: backend_is_cuda = .true.
#else
  logical, parameter :: backend_is_cuda = .false.
#endif

  type :: backend_runtime_t
    !! Test helper that creates the correct backend and allocator
    !! for the current build.  The pointer components (backend, allocator,
    !! host_allocator) alias concrete members of this same object, so
    !! instances must be declared with the TARGET attribute and must
    !! never be copied or assigned - the pointers would dangle.
    class(base_backend_t), pointer :: backend => null()
    class(allocator_t), pointer :: allocator => null()
    type(allocator_t), pointer :: host_allocator => null()
    character(len=8) :: backend_name = ''
#ifdef CUDA
    type(cuda_backend_t) :: cuda_backend
    type(cuda_allocator_t) :: cuda_allocator
#elif defined(OMP_TGT)
    type(omptgt_backend_t) :: omptgt_backend
    type(omptgt_allocator_t) :: omptgt_allocator
#else
    type(omp_backend_t) :: omp_backend
    type(allocator_t) :: omp_allocator
#endif
    type(allocator_t) :: host_allocator_storage
  contains
    procedure :: init
  end type backend_runtime_t

contains

#ifdef CUDA
  subroutine select_cuda_device(nrank, devnum)
    !! Select a CUDA device round-robin by MPI rank.  The optional result
    !! returns the device that CUDA reports as current after selection.
    integer, intent(in) :: nrank
    integer, optional, intent(out) :: devnum

    integer :: ierr, ndevs, selected_device

    ierr = cudaGetDeviceCount(ndevs)
    if (ndevs < 1) then
      error stop 'select_cuda_device: no CUDA devices available'
    end if

    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(selected_device)

    if (present(devnum)) devnum = selected_device
  end subroutine select_cuda_device
#endif

  subroutine init(self, mesh, separate_host_allocator)
    class(backend_runtime_t), target, intent(inout) :: self
    type(mesh_t), target, intent(inout) :: mesh
    logical, optional, intent(in) :: separate_host_allocator

    integer :: dims(3)
#if !defined(CUDA) && !defined(OMP_TGT)
    logical :: need_separate_host_allocator
#endif

#ifdef CUDA
    integer :: ierr, nrank
#endif

    dims = mesh%get_dims(VERT)
#if !defined(CUDA) && !defined(OMP_TGT)
    need_separate_host_allocator = .false.
    if (present(separate_host_allocator)) then
      need_separate_host_allocator = separate_host_allocator
    end if
#endif

#ifdef CUDA
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call select_cuda_device(nrank)

    self%backend_name = 'CUDA'
    self%cuda_allocator = cuda_allocator_t(dims, SZ)
    self%allocator => self%cuda_allocator
    ! CUDA always needs a separate host allocator because the main
    ! allocator lives in device memory.  The solver and case objects
    ! (xcompact.f90) unconditionally pass host_allocator, so it must
    ! be valid regardless of whether the caller requested one.
    self%host_allocator_storage = allocator_t(dims, SZ)
    self%host_allocator => self%host_allocator_storage
    self%cuda_backend = cuda_backend_t(mesh, self%allocator)
    self%backend => self%cuda_backend
#elif defined(OMP_TGT)
    self%backend_name = 'OMP_TGT'
    self%omptgt_allocator = omptgt_allocator_t(dims, SZ)
    self%allocator => self%omptgt_allocator
    ! omp_tgt manages offload memory, so always use a separate host
    ! allocator (the omptgt allocator cannot be aliased by a plain
    ! type(allocator_t) pointer).
    self%host_allocator_storage = allocator_t(dims, SZ)
    self%host_allocator => self%host_allocator_storage
    self%omptgt_backend = omptgt_backend_t(mesh, self%allocator)
    self%backend => self%omptgt_backend
#else
    self%backend_name = 'OMP'
    self%omp_allocator = allocator_t(dims, SZ)
    self%allocator => self%omp_allocator

    if (need_separate_host_allocator) then
      self%host_allocator_storage = allocator_t(dims, SZ)
      self%host_allocator => self%host_allocator_storage
    else
      self%host_allocator => self%omp_allocator
    end if

    self%omp_backend = omp_backend_t(mesh, self%allocator)
    self%backend => self%omp_backend
#endif

  end subroutine init

end module m_backend_runtime
