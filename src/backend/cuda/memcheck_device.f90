module m_cuda_memcheck_device
  !! CUDA implementation of the x3d2-memcheck device seam: the only
  !! memcheck file that uses cudafor or the m_cuda_* modules.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice, cudaMemGetInfo, &
                     cuda_count_kind
  use m_common, only: dp, i8
  use m_cuda_common, only: SZ
  use m_cuda_memory_estimate, only: check_status, context_floor_bytes, &
                                    fft_workspace_bytes_query
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_memcheck_device, only: memcheck_device_t
  use m_memcheck_estimate, only: to_gib
  implicit none

  private
  public :: cuda_memcheck_device_t

  type, extends(memcheck_device_t) :: cuda_memcheck_device_t
    private
    !> Queried once, at this process's own rank count (always 1 - see
    !> ensure_fft_query). These are real, live cudaMemGetInfo measurements
    !> for a SINGLE-RANK communicator; there is no way to measure the true
    !> ng>1 cost without actually running under mpirun -n ng, which this
    !> single-process estimator does not do. overhead_query_bytes instead
    !> extrapolates ng>1 by the 1/ng law confirmed by a real mpirun -n 2
    !> xcompact run on examples/TGV/input.x3d (2026-09-10: 3339 MiB/GPU
    !> measured vs 3.21 GiB predicted this way) - ng1_worksize_bytes (plain
    !> cuFFT path only) and ng1_heap_bytes (cuFFTMp's NVSHMEM heap) both
    !> scale this way; ng1_xtdesc_bytes (cuFFTMp's distributed data buffer)
    !> does not - it is added only at ng=1, since it was measured to cost
    !> 0 extra at ng=2 (drawn from heap space the plan already reserved).
    integer(i8) :: ng1_worksize_bytes, ng1_heap_bytes, ng1_xtdesc_bytes
    logical :: ng1_used_cufftmp
    !> Set by ensure_fft_query, the cudaMemGetInfo delta across the whole
    !> throwaway plan query (create+destroy) - if cufftDestroy fully released
    !> the NVSHMEM heap this reserved, this should read ~0.
    real(dp) :: ng1_query_residual_gib
    logical :: ng1_query_ran = .false.
    type(cuda_allocator_t) :: cuda_allocator
    type(cuda_backend_t) :: cuda_backend
  contains
    procedure :: init
    procedure :: mem_info
    procedure :: sz => pencil_size
    procedure :: overhead_floor_bytes
    procedure :: overhead_query_bytes
    procedure :: make_backend
  end type cuda_memcheck_device_t

contains

  subroutine init(self, ok, msg)
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    integer :: ierr, ndevs

    ok = .false.
    ierr = cudaGetDeviceCount(ndevs)
    if (ierr /= 0 .or. ndevs < 1) then
      write (msg, '(a,i0,a,i0,a)') 'x3d2-memcheck: no usable CUDA device &
        &(cudaGetDeviceCount code ', ierr, ', devices ', ndevs, ')'
      return
    end if
    ierr = cudaSetDevice(0)
    if (ierr /= 0) then
      write (msg, '(a,i0,a)') 'x3d2-memcheck: no usable CUDA device &
        &(cudaSetDevice failed, code ', ierr, ')'
      return
    end if
    ok = .true.
  end subroutine init

  subroutine mem_info(self, total_bytes, free_bytes)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8), intent(out) :: total_bytes, free_bytes
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo')
    total_bytes = int(total_b, i8)
    free_bytes = int(free_b, i8)
  end subroutine mem_info

  function pencil_size(self) result(n)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer :: n

    n = SZ
  end function pencil_size

  function overhead_floor_bytes(self, ng, bc_is_110) result(nbytes8)
    !! Assumes cuFFTMp (the realistic case for every BC the solver attempts
    !! it for) except for 110, which always uses plain cuFFT.
    class(cuda_memcheck_device_t), intent(in) :: self
    integer, intent(in) :: ng
    logical, intent(in) :: bc_is_110
    integer(i8) :: nbytes8
    logical :: uses_cufftmp

    uses_cufftmp = .not. bc_is_110
    nbytes8 = context_floor_bytes(ng, uses_cufftmp)
  end function overhead_floor_bytes

  subroutine ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    !! Runs the real (single-rank) Tier 2 plan query at most once per
    !! process and caches the result in ng1_worksize_bytes/ng1_heap_bytes/
    !! ng1_xtdesc_bytes/ng1_used_cufftmp, plus the cudaMemGetInfo delta
    !! across the whole call in ng1_query_residual_gib (see its
    !! declaration).
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    integer, intent(in) :: cdims(3)
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b

    if (self%ng1_query_ran) return
    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before FFT probe)')
    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   is_root, self%ng1_worksize_bytes, &
                                   self%ng1_heap_bytes, &
                                   self%ng1_xtdesc_bytes, &
                                   self%ng1_used_cufftmp)
    ierr = cudaMemGetInfo(free_after, total_b)
    call check_status(ierr, 'cudaMemGetInfo (after FFT probe)')
    self%ng1_query_residual_gib = to_gib(int(free_before, i8) - &
                                         int(free_after, i8))
    self%ng1_query_ran = .true.
  end subroutine ensure_fft_query

  function overhead_query_bytes(self, ng, bc_is_100, bc_is_110, cdims, &
                                is_root, residual_gib) result(nbytes8)
    class(cuda_memcheck_device_t), intent(inout) :: self
    integer, intent(in) :: ng, cdims(3)
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    real(dp), intent(out) :: residual_gib
    integer(i8) :: nbytes8

    integer(i8) :: worksize_bytes, xtdesc_bytes, heap_bytes, context_bytes

    call ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    residual_gib = self%ng1_query_residual_gib
    ! context_bytes: the CUDA-context-alone baseline (no ng argument
    ! matters here - context_floor_bytes only adds its hardcoded NVSHMEM
    ! heap when uses_cufftmp=.true., so passing .false. always yields just
    ! the context). The cuFFTMp heap itself now comes from the live
    ! ng1_heap_bytes measurement below instead of that hardcoded constant.
    context_bytes = context_floor_bytes(ng, .false.)
    if (self%ng1_used_cufftmp) then
      ! Live-measured (fft_workspace_bytes_query): worksize is carved out
      ! of the NVSHMEM heap for cuFFTMp, not a separate allocation - do
      ! not add it on top of heap_bytes (see that function's docstring).
      worksize_bytes = 0_i8
      heap_bytes = int(real(self%ng1_heap_bytes, dp)/real(ng, dp), i8)
      ! xtdesc cost was measured to be 0 extra at ng=2 (drawn from heap
      ! space the plan already reserved) - add the ng=1 value only there.
      xtdesc_bytes = merge(self%ng1_xtdesc_bytes, 0_i8, ng == 1)
    else
      ! See ng1_worksize_bytes' declaration: this is an extrapolation for
      ! ng>1, not an independent per-ng measurement (though in practice
      ! this branch only ever runs at ng=1 - plain cuFFT is only reached
      ! by the 110 BC or a cuFFTMp fallback, neither of which supports
      ! nproc>1 in this solver).
      worksize_bytes = int(real(self%ng1_worksize_bytes, dp)/real(ng, dp), &
                           i8)
      heap_bytes = 0_i8
      xtdesc_bytes = 0_i8
    end if
    nbytes8 = worksize_bytes + xtdesc_bytes + heap_bytes + context_bytes
  end function overhead_query_bytes

  subroutine make_backend(self, mesh, dims, allocator, backend)
    class(cuda_memcheck_device_t), target, intent(inout) :: self
    type(mesh_t), target, intent(inout) :: mesh
    integer, intent(in) :: dims(3)
    class(allocator_t), pointer, intent(out) :: allocator
    class(base_backend_t), pointer, intent(out) :: backend

    self%cuda_allocator = cuda_allocator_t(dims, SZ)
    allocator => self%cuda_allocator
    self%cuda_backend = cuda_backend_t(mesh, allocator)
    backend => self%cuda_backend
  end subroutine make_backend

end module m_cuda_memcheck_device
