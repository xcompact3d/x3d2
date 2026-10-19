module m_cuda_memcheck_device
  !! The CUDA device of x3d2-memcheck: the FFT plan probe, the memory
  !! readings and the backend of the real build.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice, cudaMemGetInfo, &
                     cuda_count_kind, int_ptr_kind
  use cufft, only: cufftDestroy, CUFFT_D2Z, CUFFT_Z2D, CUFFT_R2C, CUFFT_C2R
  use cufftXt, only: cufftXtMalloc, cufftXtFree, CUFFT_XT_FORMAT_INPLACE, &
                     cudaLibXtDesc
  use m_common, only: dp, i8, is_sp
  use m_cuda_common, only: SZ
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t, create_fft_plan
  implicit none

  private
  public :: cuda_memcheck_device_t, busy_card_gib

  !> Card usage at start (GiB) above which another process is assumed to
  !> share the card, so the measured context size is not trusted.
  real(dp), parameter :: busy_card_gib = 0.5_dp

  type :: cuda_memcheck_device_t
    private
    !> Card memory used right after cudaSetDevice: this process's context,
    !> plus any other process on the card.
    integer(i8) :: start_used_bytes = 0_i8
    !> The single-rank FFT plan query, run once by ensure_fft_query. The
    !> heap scales with 1/num_gpus for num_gpus>1. The worksize (plain
    !> cuFFT, num_gpus=1 only) and the xtdesc buffer (added at num_gpus=1
    !> only) do not.
    integer(i8) :: ng1_worksize_bytes = 0_i8, ng1_heap_bytes = 0_i8, &
                   ng1_xtdesc_bytes = 0_i8
    logical :: ng1_used_cufftmp = .false.
    !> The free memory change across the plan query, ~0 when cufftDestroy
    !> released the heap.
    real(dp) :: ng1_query_residual_gib = 0._dp
    !> Whether the plan query ran, and whether it returned a failure status.
    !> The ng1_* values are then not valid and overhead_query_bytes uses the
    !> formula.
    logical :: ng1_query_ran = .false.
    logical :: ng1_query_failed = .false.
    !> Whether the real build used cuFFTMp, once it is known.
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

    if (real(self%start_used_bytes, dp)/1024._dp**3 <= busy_card_gib) then
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

  function overhead_floor_bytes(self, num_gpus, bc_is_110) result(nbytes8)
    !! The formula overhead at num_gpus GPUs: context_floor_bytes with the
    !! context term swapped for context_term_bytes. It assumes cuFFTMp,
    !! except for the 110 BC (plain cuFFT).
    class(cuda_memcheck_device_t), intent(in) :: self
    integer, intent(in) :: num_gpus
    logical, intent(in) :: bc_is_110
    integer(i8) :: nbytes8
    logical :: uses_cufftmp

    uses_cufftmp = .not. bc_is_110
    nbytes8 = context_floor_bytes(num_gpus, uses_cufftmp) - &
              context_floor_bytes(num_gpus, .false.) + &
              self%context_term_bytes()
  end function overhead_floor_bytes

  subroutine ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    !! Runs the single-rank FFT plan query once and caches its result and the
    !! free memory change across it.
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    integer, intent(in) :: cdims(3)
    integer :: ierr, stat
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b

    if (self%ng1_query_ran) return
    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before FFT probe)')
    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   is_root, self%ng1_worksize_bytes, &
                                   self%ng1_heap_bytes, &
                                   self%ng1_xtdesc_bytes, &
                                   self%ng1_used_cufftmp, stat)
    self%ng1_query_failed = stat /= 0
    ierr = cudaMemGetInfo(free_after, total_b)
    call check_status(ierr, 'cudaMemGetInfo (after FFT probe)')
    self%ng1_query_residual_gib = real(int(free_before, i8) - &
                                       int(free_after, i8), dp)/1024._dp**3
    self%ng1_query_ran = .true.
  end subroutine ensure_fft_query

  function overhead_query_bytes(self, num_gpus, bc_is_100, bc_is_110, cdims, &
                                is_root) result(nbytes8)
    class(cuda_memcheck_device_t), intent(inout) :: self
    integer, intent(in) :: num_gpus, cdims(3)
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    integer(i8) :: nbytes8

    integer(i8) :: worksize_bytes, xtdesc_bytes, heap_bytes, context_bytes

    call ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    if (self%ng1_query_failed) then
      nbytes8 = self%overhead_floor_bytes(num_gpus, bc_is_110)
      return
    end if
    ! The context is the live reading on an idle card (see
    ! context_term_bytes). The heap comes from the plan query.
    context_bytes = self%context_term_bytes()
    if (self%ng1_used_cufftmp) then
      ! The worksize is part of the heap, so it is not added.
      worksize_bytes = 0_i8
      heap_bytes = int(real(self%ng1_heap_bytes, dp)/real(num_gpus, dp), i8)
      ! The xtdesc buffer costs nothing extra at num_gpus=2: add it at
      ! num_gpus=1 only.
      xtdesc_bytes = merge(self%ng1_xtdesc_bytes, 0_i8, num_gpus == 1)
    else
      ! Plain cuFFT: the 110 BC or a cuFFTMp fallback. Neither supports
      ! nproc>1. The first plan query can run inside an estimate for
      ! num_gpus>1, so the callers ask num_gpus_unsupported_reason again and
      ! drop that value.
      worksize_bytes = self%ng1_worksize_bytes
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

    known = self%use_cufftmp_known .or. &
            (self%ng1_query_ran .and. .not. self%ng1_query_failed)
  end function fft_path_known

  function query_residual_gib(self, residual_gib) result(ran)
    class(cuda_memcheck_device_t), intent(in) :: self
    real(dp), intent(out) :: residual_gib
    logical :: ran

    ran = self%ng1_query_ran .and. .not. self%ng1_query_failed
    residual_gib = 0._dp
    if (ran) residual_gib = self%ng1_query_residual_gib
  end function query_residual_gib

  function uses_distributed_fft(self) result(distributed)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: distributed

    if (self%use_cufftmp_known) then
      distributed = self%use_cufftmp
    else
      distributed = self%ng1_used_cufftmp
    end if
  end function uses_distributed_fft

  subroutine check_status(ierr, what, stat)
    !! Stops with a labelled message on a non-zero CUDA/cuFFT return code.
    !! With stat present, returns the code in stat instead.
    integer, intent(in) :: ierr
    character(*), intent(in) :: what
    integer, optional, intent(out) :: stat
    character(len=128) :: msg

    if (present(stat)) stat = ierr
    if (ierr == 0 .or. present(stat)) return
    write (msg, '(a,a,a,i0,a)') 'x3d2-memcheck: ', trim(what), &
      ' failed (code ', ierr, ')'
    error stop trim(msg)
  end subroutine check_status

  pure function context_floor_bytes(num_gpus, uses_cufftmp) result(nbytes8)
    !! Grid-independent memory reserved before any field/FFT-data bytes:
    !! the CUDA context itself, plus (when cuFFTMp is used) an estimate of
    !! the NVSHMEM symmetric heap. Two different callers use this function
    !! for two different purposes:
    !!   - memcheck.f90's formula floor (before any GPU plan is queried) and
    !!     its final num_gpus>1 estimate for the plain-cuFFT path (which never
    !!     runs cuFFTMp, so uses_cufftmp=.false. there and only the
    !!     grid-independent cuda_context_gib term applies).
    !!   - memcheck.f90's --build path (report_measured_table), which
    !!     already has a REAL measured num_gpus=1 overhead (from an actual case
    !!     build) and only needs this function's num_gpus-scaling RATIO
    !!     (context_floor_bytes(num_gpus,.true.) -
    !!     context_floor_bytes(1,.true.))
    !!     to correct that real measurement down for num_gpus>1 - confirmed
    !!     correct by a real mpirun -n 2 xcompact run on
    !!     examples/TGV/input.x3d, 2026-09-10 (3339 MiB/GPU measured vs
    !!     3.24 GiB predicted using this ratio, vs the 4.69 GiB this table
    !!     used to report when it wrongly held the num_gpus=1 overhead
    !!     constant).
    !! memcheck.f90's own cuFFTMp probe estimate does NOT use this
    !! function's heap constant any more - it live-measures the real
    !! num_gpus=1 heap delta itself via fft_workspace_bytes_query (bracketing
    !! the real cufftMakePlan3D calls with cudaMemGetInfo, same as it
    !! already did for worksize), then applies the SAME 1/num_gpus ratio this
    !! function encodes to extrapolate to num_gpus>1. That keeps the num_gpus=1
    !! absolute number self-calibrating against NVHPC/CUDA/NVSHMEM version
    !! drift, while both callers still agree on the same measured
    !! num_gpus-scaling shape.
    !! Constants below measured directly with cudaMemGetInfo on kolmogorov
    !! (A100-PCIE-40GB, NVHPC 25.3/CUDA 12.8, NVSHMEM_SYMMETRIC_SIZE unset
    !! i.e. NVSHMEM's built-in default), 2026-09-09:
    !!   (memcheck now uses the live context size on an idle card; this
    !!   constant is its shared-card fallback and the --static floor.)
    !!   CUDA context alone: ~0.42 GiB (39.493 GiB card total - 39.078 GiB
    !!     free immediately after cudaSetDevice, before any cuFFT/cuFFTMp
    !!     call - see memcheck.f90's own "Card memory: 39.49 GiB" line,
    !!     query_card_gib, for the card-total cross-check)
    !!   + NVSHMEM heap at num_gpus=1: 2972 MiB from cufftMakePlan3D itself,
    !!     confirmed IDENTICAL at two very different grid sizes (64^3 and
    !!     256x256x512) AND identical regardless of the plan's own reported
    !!     worksize (2 MiB vs 258 MiB for those two grids) - grid- and
    !!     worksize-independent. This means the worksize cufftMakePlan3D
    !!     reports for a cuFFTMp plan does NOT correspond to a separate,
    !!     additional allocation - it is carved out of this same heap, so
    !!     the cuFFTMp probe's worksize must NOT be added on top of this
    !!     for the cuFFTMp path (see fft_workspace_bytes_query - worksize is
    !!     only added for .not. used_cufftmp). An earlier version of this
    !!     estimator added both and looked accurate on TGV by coincidence
    !!     (a second, separately-missing term happened to approximately
    !!     cancel the double-count); adding a real second FFT plan's
    !!     worksize on top made that cancellation worse (6.10 GiB -> 6.23
    !!     GiB estimate, further over the 6.007 GiB measured value), which
    !!     is what exposed the double-count.
    !!   num_gpus=1/num_gpus=2 heap: 2972/1484 MiB, i.e. ~2972/num_gpus -
    !!     confirmed at num_gpus=1,2 only (kolmogorov has 2 A100s);
    !!     num_gpus=4/8 extrapolate the same 1/num_gpus law, unconfirmed
    !!     beyond num_gpus=2.
    !! This whole function is hardware/library-specific and must be
    !! re-measured if the target GPU model or NVHPC/CUDA/NVSHMEM version
    !! changes - or, for memcheck.f90's own num_gpus=1 cuFFTMp number, simply
    !! re-run (it live-measures its own num_gpus=1 baseline).
    integer, intent(in) :: num_gpus
    logical, intent(in) :: uses_cufftmp
    integer(i8) :: nbytes8
    real(dp), parameter :: cuda_context_gib = 0.42_dp
    real(dp), parameter :: nvshmem_heap_ng1_gib = 2972._dp/1024._dp

    nbytes8 = int(cuda_context_gib*1024._dp**3, i8)
    if (uses_cufftmp) &
      nbytes8 = nbytes8 + &
                int(nvshmem_heap_ng1_gib/real(num_gpus, dp)*1024._dp**3, i8)
  end function context_floor_bytes

  subroutine fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, &
                                       try_cufftmp, is_root, &
                                       worksize_bytes, heap_bytes, &
                                       xtdesc_bytes, used_cufftmp, stat)
    !! Creates and destroys the forward and backward FFT plans the way the
    !! solver does, for the cell dimensions cdims, and for cuFFTMp also the
    !! xtdesc buffer, each bracketed with cudaMemGetInfo. worksize_bytes is
    !! the plan workspace (used only on the plain cuFFT path, in cuFFTMp it
    !! is part of the heap), heap_bytes the memory the plans took,
    !! xtdesc_bytes the xtdesc buffer. The 110 BC never uses cuFFTMp. stat is
    !! the first non-zero cuFFT/CUDA code, or 0, and the other outputs are
    !! then not valid. A plan that aborts inside the library never returns.
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    logical, intent(in) :: try_cufftmp
    logical, intent(in) :: is_root
    integer(i8), intent(out) :: worksize_bytes, heap_bytes, xtdesc_bytes
    logical, intent(out) :: used_cufftmp
    integer, intent(out) :: stat

    integer :: plan_fw, plan_bw, fft_n1, fft_n2, fft_n3
    integer :: fw_plan_type, bw_plan_type, ierr
    integer(int_ptr_kind()) :: worksize_fw, worksize_bw
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b
    type(cudaLibXtDesc), pointer :: xtdesc
    logical :: fw_was_cufftmp

    ! The plan dimensions follow the permutation of the solver's init().
    if (bc_is_100) then
      fft_n1 = cdims(2); fft_n2 = cdims(1); fft_n3 = cdims(3)
    else if (bc_is_110) then
      fft_n1 = cdims(3); fft_n2 = cdims(1); fft_n3 = cdims(2)
    else
      fft_n1 = cdims(1); fft_n2 = cdims(2); fft_n3 = cdims(3)
    end if

    fw_plan_type = merge(CUFFT_R2C, CUFFT_D2Z, is_sp)
    bw_plan_type = merge(CUFFT_C2R, CUFFT_Z2D, is_sp)

    used_cufftmp = try_cufftmp .and. (.not. bc_is_110)
    stat = 0

    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before plans)')
    call create_fft_plan(plan_fw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         fw_plan_type, is_root, 'memcheck probe fwd', &
                         worksize_fw, stat)
    if (stat /= 0) return
    fw_was_cufftmp = used_cufftmp
    ! A fallback to plain cuFFT sets used_cufftmp to false. A forward plan
    ! built on cuFFTMp is then rebuilt on plain cuFFT, as init() does.
    call create_fft_plan(plan_bw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         bw_plan_type, is_root, 'memcheck probe bwd', &
                         worksize_bw, stat)
    if (stat /= 0) then
      ierr = cufftDestroy(plan_fw)
      return
    end if
    if (fw_was_cufftmp .and. (.not. used_cufftmp)) then
      ierr = cufftDestroy(plan_fw)
      call check_status(ierr, 'cufftDestroy (fwd, cuFFTMp fallback rebuild)')
      call create_fft_plan(plan_fw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                           fw_plan_type, is_root, &
                           'memcheck probe fwd (cuFFT rebuild)', worksize_fw, &
                           stat)
      if (stat /= 0) then
        ierr = cufftDestroy(plan_bw)
        return
      end if
    end if
    ierr = cudaMemGetInfo(free_after, total_b)
    call check_status(ierr, 'cudaMemGetInfo (after plans)')
    worksize_bytes = int(worksize_fw, i8) + int(worksize_bw, i8)
    heap_bytes = int(free_before, i8) - int(free_after, i8)

    if (used_cufftmp) then
      ierr = cudaMemGetInfo(free_before, total_b)
      call check_status(ierr, 'cudaMemGetInfo (before xtdesc)')
      ierr = cufftXtMalloc(plan_fw, xtdesc, CUFFT_XT_FORMAT_INPLACE)
      call check_status(ierr, 'cufftXtMalloc', stat)
      if (stat /= 0) then
        ierr = cufftDestroy(plan_fw)
        ierr = cufftDestroy(plan_bw)
        return
      end if
      ierr = cudaMemGetInfo(free_after, total_b)
      call check_status(ierr, 'cudaMemGetInfo (after xtdesc)')
      xtdesc_bytes = int(free_before, i8) - int(free_after, i8)
      ierr = cufftXtFree(xtdesc)
      call check_status(ierr, 'cufftXtFree')
    else
      xtdesc_bytes = 0_i8
    end if

    ierr = cufftDestroy(plan_fw)
    call check_status(ierr, 'cufftDestroy (fwd)')
    ierr = cufftDestroy(plan_bw)
    call check_status(ierr, 'cufftDestroy (bwd)')
  end subroutine fft_workspace_bytes_query

end module m_cuda_memcheck_device
