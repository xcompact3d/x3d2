module m_backend_runtime
  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Wtime

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, VERT
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t

#ifdef CUDA
  use cudafor
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: cuda_penta_compact => exec_dist_penta_compact, &
                              cuda_penta_periodic => exec_dist_penta_periodic
  use m_cuda_tdsops, only: cuda_tdsops_t
#elif defined(OMP_TGT)
  use omp_lib, only: omp_get_num_devices, omp_set_default_device, &
                     omp_get_default_device
  use m_omptgt_common, only: SZ
  use m_omptgt_allocator, only: omptgt_allocator_t
  use m_omptgt_backend, only: omptgt_backend_t
  use m_omp_exec_dist, only: exec_dist_penta_compact, exec_dist_penta_periodic
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_penta_compact, exec_dist_penta_periodic
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
  subroutine select_device(nrank, devnum)
    !! Select a CUDA device round-robin by MPI rank.  The optional result
    !! returns the device that CUDA reports as current after selection.
    integer, intent(in) :: nrank
    integer, optional, intent(out) :: devnum

    integer :: ierr, ndevs, selected_device

    ierr = cudaGetDeviceCount(ndevs)
    if (ndevs < 1) then
      error stop 'select_device: no CUDA devices available'
    end if

    ierr = cudaSetDevice(mod(nrank, ndevs))
    ierr = cudaGetDevice(selected_device)

    if (present(devnum)) devnum = selected_device
  end subroutine select_device
#elif defined(OMP_TGT)
  subroutine select_device(nrank, devnum)
    !! Select an OpenMP target device round-robin by MPI rank.  The optional
    !! result returns the device OpenMP reports as default after selection.
    integer, intent(in) :: nrank
    integer, optional, intent(out) :: devnum

    integer :: ndevs, selected_device

    ndevs = omp_get_num_devices()
    if (ndevs < 1) then
      error stop 'select_device: no TGT devices available'
    end if

    call omp_set_default_device(mod(nrank, ndevs))
    selected_device = omp_get_default_device()

    if (present(devnum)) devnum = selected_device
  end subroutine select_device
#endif

  subroutine init(self, mesh, separate_host_allocator)
    class(backend_runtime_t), target, intent(inout) :: self
    type(mesh_t), target, intent(inout) :: mesh
    logical, optional, intent(in) :: separate_host_allocator

    integer :: dims(3)
#if !defined(CUDA) && !defined(OMP_TGT)
    logical :: need_separate_host_allocator
#endif

#if defined(CUDA) || defined(OMP_TGT)
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
    call select_device(nrank)

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
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call select_device(nrank)

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

  subroutine penta_solve(du, u, u_recv_s, u_recv_e, tdsops, periodic)
    !! Run the single-subdomain pentadiagonal compact solve on host arrays of
    !! shape (backend_sz, n, n_groups) with the halos supplied by the caller.
    !! The pentadiagonal kernels are not reachable through base_backend_t yet,
    !! so this hides the backend-specific call (and, on CUDA, the device
    !! copies) from the tests.  tdsops must come from backend%alloc_tdsops.
    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e
    class(tdsops_t), intent(in) :: tdsops
    logical, intent(in) :: periodic

    real(dp) :: elapsed

    call run_penta(du, u, u_recv_s, u_recv_e, tdsops, periodic, 0, 1, &
                   elapsed)
  end subroutine penta_solve

  real(dp) function penta_benchmark(u, u_recv_s, u_recv_e, tdsops, periodic, &
                                    n_warmup, n_iters) result(elapsed)
    !! As penta_solve, but repeat the solve and return the wall time of the
    !! n_iters timed repetitions, with the data resident on the device.
    real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e
    class(tdsops_t), intent(in) :: tdsops
    logical, intent(in) :: periodic
    integer, intent(in) :: n_warmup, n_iters

    real(dp), allocatable :: du(:, :, :)

    allocate (du, mold=u)
    call run_penta(du, u, u_recv_s, u_recv_e, tdsops, periodic, n_warmup, &
                   n_iters, elapsed)
  end function penta_benchmark

  subroutine run_penta(du, u, u_recv_s, u_recv_e, tdsops, periodic, &
                       n_warmup, n_iters, elapsed)
    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e
    class(tdsops_t), intent(in) :: tdsops
    logical, intent(in) :: periodic
    integer, intent(in) :: n_warmup, n_iters
    real(dp), intent(out) :: elapsed

    real(dp) :: tstart
    integer :: iter
#ifdef CUDA
    real(dp), device, allocatable, dimension(:, :, :) :: du_d, u_d, &
                                                         u_recv_s_d, u_recv_e_d
    type(dim3) :: blocks, threads
    integer :: ierr

    allocate (du_d, mold=du)
    u_d = u
    u_recv_s_d = u_recv_s
    u_recv_e_d = u_recv_e
    blocks = dim3(size(u, 3), 1, 1)
    threads = dim3(SZ, 1, 1)

    select type (tdsops)
    type is (cuda_tdsops_t)
      do iter = 1, n_warmup
        call launch(tdsops)
      end do
      ierr = cudaDeviceSynchronize()
      tstart = MPI_Wtime()
      do iter = 1, n_iters
        call launch(tdsops)
      end do
      ierr = cudaDeviceSynchronize()
      elapsed = MPI_Wtime() - tstart
    class default
      error stop 'penta_solve: tdsops is not a cuda_tdsops_t'
    end select

    du = du_d
#else
    select type (tdsops)
    type is (tdsops_t)
      do iter = 1, n_warmup
        call launch(tdsops)
      end do
      tstart = MPI_Wtime()
      do iter = 1, n_iters
        call launch(tdsops)
      end do
      elapsed = MPI_Wtime() - tstart
    class default
      error stop 'penta_solve: tdsops is not a tdsops_t'
    end select
#endif

  contains

#ifdef CUDA
    subroutine launch(tdsops_dev)
      type(cuda_tdsops_t), intent(in) :: tdsops_dev

      if (periodic) then
        call cuda_penta_periodic(du_d, u_d, u_recv_s_d, u_recv_e_d, &
                                 tdsops_dev, blocks, threads)
      else
        call cuda_penta_compact(du_d, u_d, u_recv_s_d, u_recv_e_d, &
                                tdsops_dev, blocks, threads)
      end if
    end subroutine launch
#else
    subroutine launch(tdsops_host)
      type(tdsops_t), intent(in) :: tdsops_host

      if (periodic) then
        call exec_dist_penta_periodic(du, u, u_recv_s, u_recv_e, &
                                      tdsops_host, size(u, 3))
      else
        call exec_dist_penta_compact(du, u, u_recv_s, u_recv_e, &
                                     tdsops_host, size(u, 3))
      end if
    end subroutine launch
#endif

  end subroutine run_penta

end module m_backend_runtime
