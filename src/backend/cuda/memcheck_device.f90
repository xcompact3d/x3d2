module m_cuda_memcheck_device
  !! The CUDA device of x3d2-memcheck: the only memcheck file that uses
  !! cudafor, cufft or the m_cuda_* modules. It holds the FFT plan probe,
  !! the memory readings and the backend of the real build.
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
    !> measured vs 3.21 GiB predicted this way) - ng1_heap_bytes (cuFFTMp's
    !> NVSHMEM heap) scales this way, ng1_worksize_bytes (plain cuFFT path
    !> only, ng=1 only) does not; ng1_xtdesc_bytes (cuFFTMp's distributed
    !> data buffer) does not either - it is added only at ng=1, since it was
    !> measured to cost 0 extra at ng=2 (drawn from heap space the plan
    !> already reserved).
    integer(i8) :: ng1_worksize_bytes = 0_i8, ng1_heap_bytes = 0_i8, &
                   ng1_xtdesc_bytes = 0_i8
    logical :: ng1_used_cufftmp = .false.
    !> Set by ensure_fft_query, the cudaMemGetInfo delta across the whole
    !> throwaway plan query (create+destroy) - if cufftDestroy fully released
    !> the NVSHMEM heap this reserved, this should read ~0.
    real(dp) :: ng1_query_residual_gib = 0._dp
    logical :: ng1_query_ran = .false.
    !> Set by ensure_fft_query when the throwaway plan query returned a
    !> failure status; the ng1_* values are then not valid and
    !> overhead_query_bytes uses the formula (overhead_floor_bytes).
    logical :: ng1_query_failed = .false.
    !> Whether a --build's real build actually used cuFFTMp (only known
    !> once the real solver is constructed) - the only real "cuFFTMp
    !> available" signal in this codebase, more authoritative than the cuFFTMp
    !> probe's throwaway plan query (ng1_used_cufftmp above) for gating the
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
    !! Runs the real (single-rank) cuFFTMp probe plan query at most once per
    !! process and caches the result in ng1_worksize_bytes/ng1_heap_bytes/
    !! ng1_xtdesc_bytes/ng1_used_cufftmp, plus the cudaMemGetInfo delta
    !! across the whole call in ng1_query_residual_gib (see its
    !! declaration).
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

  function overhead_query_bytes(self, ng, bc_is_100, bc_is_110, cdims, &
                                is_root) result(nbytes8)
    class(cuda_memcheck_device_t), intent(inout) :: self
    integer, intent(in) :: ng, cdims(3)
    logical, intent(in) :: bc_is_100, bc_is_110, is_root
    integer(i8) :: nbytes8

    integer(i8) :: worksize_bytes, xtdesc_bytes, heap_bytes, context_bytes

    call ensure_fft_query(self, bc_is_100, bc_is_110, cdims, is_root)
    if (self%ng1_query_failed) then
      nbytes8 = self%overhead_floor_bytes(ng, bc_is_110)
      return
    end if
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
      ! Plain cuFFT is only reached by the 110 BC or a cuFFTMp fallback.
      ! Neither supports nproc>1 (ng_unsupported_reason). The first probe
      ! can run inside an estimate for ng>1, so this branch can see ng>1
      ! once. The callers ask ng_unsupported_reason again after the
      ! estimate and drop that value.
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
    !! Abort with a labelled message if a CUDA/cuFFT return code is
    !! non-zero (0 = success) - used after every raw cudaMemGetInfo,
    !! cufftXtMalloc, cufftXtFree and cufftDestroy call in
    !! fft_workspace_bytes_query, none of which check their own ierr.
    !! When stat is present the code is returned in it instead and the
    !! program does not stop.
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

  pure function context_floor_bytes(ng, uses_cufftmp) result(nbytes8)
    !! Grid-independent memory reserved before any field/FFT-data bytes:
    !! the CUDA context itself, plus (when cuFFTMp is used) an estimate of
    !! the NVSHMEM symmetric heap. Two different callers use this function
    !! for two different purposes:
    !!   - memcheck.f90's formula floor (before any GPU plan is queried) and
    !!     its final ng>1 estimate for the plain-cuFFT path (which never
    !!     runs cuFFTMp, so uses_cufftmp=.false. there and only the
    !!     grid-independent cuda_context_gib term applies).
    !!   - memcheck.f90's --build path (report_measured_table), which
    !!     already has a REAL measured ng=1 overhead (from an actual case
    !!     build) and only needs this function's ng-scaling RATIO
    !!     (context_floor_bytes(ng,.true.) - context_floor_bytes(1,.true.))
    !!     to correct that real measurement down for ng>1 - confirmed
    !!     correct by a real mpirun -n 2 xcompact run on
    !!     examples/TGV/input.x3d, 2026-09-10 (3339 MiB/GPU measured vs
    !!     3.24 GiB predicted using this ratio, vs the 4.69 GiB this table
    !!     used to report when it wrongly held the ng=1 overhead constant).
    !! memcheck.f90's own cuFFTMp probe estimate does NOT use this
    !! function's heap constant any more - it live-measures the real
    !! ng=1 heap delta itself via fft_workspace_bytes_query (bracketing
    !! the real cufftMakePlan3D calls with cudaMemGetInfo, same as it
    !! already did for worksize), then applies the SAME 1/ng ratio this
    !! function encodes to extrapolate to ng>1. That keeps the ng=1
    !! absolute number self-calibrating against NVHPC/CUDA/NVSHMEM version
    !! drift, while both callers still agree on the same measured
    !! ng-scaling shape.
    !! Constants below measured directly with cudaMemGetInfo on kolmogorov
    !! (A100-PCIE-40GB, NVHPC 25.3/CUDA 12.8, NVSHMEM_SYMMETRIC_SIZE unset
    !! i.e. NVSHMEM's built-in default), 2026-09-09:
    !!   (memcheck now uses the live context size on an idle card; this
    !!   constant is its shared-card fallback and the --static floor.)
    !!   CUDA context alone: ~0.42 GiB (39.493 GiB card total - 39.078 GiB
    !!     free immediately after cudaSetDevice, before any cuFFT/cuFFTMp
    !!     call - see memcheck.f90's own "Card memory: 39.49 GiB" line,
    !!     query_card_gib, for the card-total cross-check)
    !!   + NVSHMEM heap at ng=1: 2972 MiB from cufftMakePlan3D itself,
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
    !!   ng=1/ng=2 heap: 2972/1484 MiB, i.e. ~2972/ng - confirmed at ng=1,2
    !!     only (kolmogorov has 2 A100s); ng=4/8 extrapolate the same 1/ng
    !!     law, unconfirmed beyond ng=2.
    !! This whole function is hardware/library-specific and must be
    !! re-measured if the target GPU model or NVHPC/CUDA/NVSHMEM version
    !! changes - or, for memcheck.f90's own ng=1 cuFFTMp number, simply
    !! re-run (it live-measures its own ng=1 baseline).
    integer, intent(in) :: ng
    logical, intent(in) :: uses_cufftmp
    integer(i8) :: nbytes8
    real(dp), parameter :: cuda_context_gib = 0.42_dp
    real(dp), parameter :: nvshmem_heap_ng1_gib = 2972._dp/1024._dp

    nbytes8 = int(cuda_context_gib*1024._dp**3, i8)
    if (uses_cufftmp) &
      nbytes8 = nbytes8 + &
                int(nvshmem_heap_ng1_gib/real(ng, dp)*1024._dp**3, i8)
  end function context_floor_bytes

  subroutine fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, &
                                       try_cufftmp, is_root, &
                                       worksize_bytes, heap_bytes, &
                                       xtdesc_bytes, used_cufftmp, stat)
    !! Standalone plan-only probe for the FFT term of the static estimator.
    !! Creates+destroys the real forward+backward plans (mirroring init()),
    !! and for cuFFTMp also creates+frees the real xtdesc data buffer
    !! (mirroring init()'s single cufftXtMalloc call on the forward plan -
    !! src/backend/cuda/poisson_fft.f90:446-447). Both the plan creation
    !! and the xtdesc allocation are bracketed with cudaMemGetInfo to
    !! measure their true cost directly on this hardware/library version,
    !! rather than relying on hardcoded constants or guessed multipliers
    !! that can drift with the NVHPC/CUDA/NVSHMEM version or the grid
    !! shape.
    !!
    !! FFT plan dimension permutation mirrors
    !! src/backend/cuda/poisson_fft.f90:299-312 (the same rule init() uses
    !! to turn global CELL dims into the plan's (n1,n2,n3) - cdims here is
    !! global CELL dims, i.e. the same quantity as poisson_fft%nx_glob/
    !! ny_glob/nz_glob, see src/poisson_fft.f90:133-134):
    !!   100: (ny, nx, nz)
    !!   110: (nz, nx, ny)
    !!   000/010: (nx, ny, nz)
    !!
    !! Both plans are created and held simultaneously, then both
    !! destroyed, mirroring init() exactly (src/backend/cuda/poisson_fft.
    !! f90:406-421 for the 110 case, :413-429 for the general case): the
    !! real solver keeps plan3D_fw and plan3D_bw alive for the whole run,
    !! each with its own auto-allocated workspace. If creating the backward
    !! plan forces a cuFFTMp->cuFFT fallback after the forward plan already
    !! built on cuFFTMp, the forward plan is destroyed and re-created on
    !! plain cuFFT, mirroring init()'s own fw_was_cufftmp rebuild
    !! (src/backend/cuda/poisson_fft.f90:429,437-442).
    !!
    !! Outputs:
    !!   worksize_bytes - forward+backward plan workspace, cuFFT's own
    !!     reported size. Used by the caller ONLY for the plain-cuFFT path
    !!     (110, or a cuFFTMp fallback): confirmed to match a real
    !!     cudaMemGetInfo delta almost exactly there. For cuFFTMp this
    !!     number is carved out of heap_bytes below, not a separate
    !!     allocation - the caller must not add both (see
    !!     context_floor_bytes' docstring for how an earlier version of
    !!     this estimator double-counted this).
    !!   heap_bytes - the real cudaMemGetInfo delta across both
    !!     create_fft_plan calls. For cuFFTMp this is dominated by the
    !!     NVSHMEM symmetric heap; for plain cuFFT it is expected to equal
    !!     worksize_bytes (both were confirmed to match a real measurement
    !!     independently during this feature's design).
    !!   xtdesc_bytes - the real cudaMemGetInfo delta from cufftXtMalloc
    !!     alone (0 for plain cuFFT, which never calls it). Grid-dependent -
    !!     this replaces an earlier hardcoded "2x spectral_slab_bytes"
    !!     guess that overshot on a triply-periodic grid (TGV: predicted
    !!     6.23 GiB vs 6.007 GiB measured) because it was calibrated off a
    !!     single, differently-shaped 100-case grid; a live measurement
    !!     removes the need to guess the multiplier per BC/grid shape.
    !! All three are single-rank (ng=1) measurements only - see their
    !! callers in memcheck.f90 for how ng>1 is handled (this solver's own
    !! nproc>1 cuFFTMp path cannot be probed from a single-rank process).
    !!
    !! The 110 case never uses cuFFTMp in this solver ("no cuFFTMp for
    !! non-periodic BCs" - src/backend/cuda/poisson_fft.f90:389-391): this
    !! is enforced here too, regardless of try_cufftmp, so the query can't
    !! silently diverge from what init() actually does.
    !!
    !! stat is the first non-zero cuFFT/CUDA code, or 0. A failure then
    !! returns instead of stopping, with the plans and buffers made so far
    !! destroyed, and the other outputs are not valid. The caller falls
    !! back to the formula. A plan that aborts inside the library (cuFFTMp
    !! on NVSHMEM) never returns, so the caller must not start a probe
    !! that cannot fit (see estimate_for_ng).
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
    ! If the forward plan fell back from cuFFTMp to plain cuFFT,
    ! used_cufftmp is now .false. - match it for the backward plan too, the
    ! same way init() rebuilds the forward plan on a backward-plan fallback
    ! (src/backend/cuda/poisson_fft.f90:424-428). The reverse case - the
    ! backward plan's own cuFFTMp attach fails after the forward plan
    ! already built on cuFFTMp - leaves the forward plan mismatched with
    ! the now-false used_cufftmp flag, so it is destroyed and re-created on
    ! plain cuFFT below, mirroring init()'s own fw_was_cufftmp rebuild.
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
