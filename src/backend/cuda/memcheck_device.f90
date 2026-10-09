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
  use m_base_case, only: base_case_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t
  use m_memcheck_device, only: memcheck_device_t, busy_card_gib
  use m_memcheck_estimate, only: to_gib
  implicit none

  private
  public :: cuda_memcheck_device_t

  type, extends(memcheck_device_t) :: cuda_memcheck_device_t
    private
    !> Card used right after cudaSetDevice (total - free): this process's
    !> context, plus any other process on the card. Set by init.
    integer(i8) :: start_used_bytes = 0_i8
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
    !> Whether a --build's real build actually used cuFFTMp (only known
    !> once the real solver is constructed) - the only real "cuFFTMp
    !> available" signal in this codebase, more authoritative than Tier 2's
    !> throwaway plan query (ng1_used_cufftmp above) for gating the
    !> measured table's ng>1 rows.
    logical :: use_cufftmp = .false., use_cufftmp_known = .false.
    type(cuda_allocator_t) :: cuda_allocator
    type(cuda_backend_t) :: cuda_backend
  contains
    procedure :: init
    procedure :: mem_info
    procedure :: context_bytes
    procedure :: context_term_bytes
    procedure :: sz => pencil_size
    procedure :: overhead_floor_bytes
    procedure :: overhead_query_bytes
    procedure :: query_residual_gib
    procedure :: make_backend
    procedure :: capture_solver_info
    procedure :: fft_path_known
    procedure :: uses_distributed_fft
  end type cuda_memcheck_device_t

contains

  subroutine init(self, ok, msg)
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    integer :: ierr, ndevs
    integer(kind=cuda_count_kind) :: free_b, total_b

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
    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo (context size)')
    self%start_used_bytes = int(total_b, i8) - int(free_b, i8)
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

  function context_bytes(self) result(nbytes8)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8) :: nbytes8

    nbytes8 = self%start_used_bytes
  end function context_bytes

  function context_term_bytes(self) result(nbytes8)
    !! The live context size (context_bytes) when the card was otherwise
    !! idle at start, else the A100 constant (context_floor_bytes without
    !! the heap).
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8) :: nbytes8

    if (to_gib(self%start_used_bytes) <= busy_card_gib) then
      nbytes8 = self%start_used_bytes
    else
      nbytes8 = context_floor_bytes(1, .false.)
    end if
  end function context_term_bytes

  function pencil_size(self) result(n)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer :: n

    n = SZ
  end function pencil_size

  function overhead_floor_bytes(self, ng, bc_is_110) result(nbytes8)
    !! Assumes cuFFTMp (the realistic case for every BC the solver attempts
    !! it for) except for 110, which always uses plain cuFFT.
    !! Formula: context_floor_bytes(ng, uses_cufftmp) with its context
    !! constant swapped for context_term_bytes, i.e.
    !! context_floor_bytes(ng, uses_cufftmp) - context_floor_bytes(ng,
    !! .false.) + context_term_bytes (the heap part is unchanged).
    class(cuda_memcheck_device_t), intent(in) :: self
    integer, intent(in) :: ng
    logical, intent(in) :: bc_is_110
    integer(i8) :: nbytes8
    logical :: uses_cufftmp

    uses_cufftmp = .not. bc_is_110
    nbytes8 = context_floor_bytes(ng, uses_cufftmp) - &
              context_floor_bytes(ng, .false.) + self%context_term_bytes()
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
                                is_root) result(nbytes8)
    class(cuda_memcheck_device_t), intent(inout) :: self
    integer, intent(in) :: ng, cdims(3)
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    integer(i8) :: nbytes8

    integer(i8) :: worksize_bytes, xtdesc_bytes, heap_bytes, context_bytes

    call ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    ! context_bytes: the CUDA-context-alone baseline, the live start-of-run
    ! measurement when the card was idle (see context_term_bytes). The
    ! cuFFTMp heap itself now comes from the live ng1_heap_bytes
    ! measurement below instead of that hardcoded constant.
    context_bytes = self%context_term_bytes()
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

  subroutine capture_solver_info(self, flow_case)
    class(cuda_memcheck_device_t), intent(inout) :: self
    class(base_case_t), intent(in) :: flow_case

    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      self%use_cufftmp = pf%use_cufftmp
      self%use_cufftmp_known = .true.
    end select
  end subroutine capture_solver_info

  function fft_path_known(self) result(known)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: known

    known = self%use_cufftmp_known
  end function fft_path_known

  function query_residual_gib(self, residual_gib) result(ran)
    class(cuda_memcheck_device_t), intent(in) :: self
    real(dp), intent(out) :: residual_gib
    logical :: ran

    ran = self%ng1_query_ran
    residual_gib = 0._dp
    if (ran) residual_gib = self%ng1_query_residual_gib
  end function query_residual_gib

  function uses_distributed_fft(self) result(distributed)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: distributed

    distributed = self%use_cufftmp
  end function uses_distributed_fft

end module m_cuda_memcheck_device
