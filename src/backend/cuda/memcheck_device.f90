module m_cuda_memcheck_device
  !! The CUDA device of x3d2-memcheck: the FFT plan probe, the memory
  !! readings and the backend of the real build.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice, cudaMemGetInfo, &
                     cuda_count_kind, int_ptr_kind
  use cufft, only: cufftDestroy, CUFFT_D2Z, CUFFT_Z2D, CUFFT_R2C, CUFFT_C2R
  use cufftXt, only: cufftXtMalloc, cufftXtFree, CUFFT_XT_FORMAT_INPLACE, &
                     cudaLibXtDesc
  use m_common, only: dp, i8, is_sp
  use m_memory_estimate, only: predicted_fft_heap_bytes, round_up_bytes, &
                               fft_xtdesc_step_bytes
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
  public :: cuda_memcheck_device_t, fft_workspace_bytes_query

  !> What one FFT probe returned. ok is false when the probe returned a
  !> failure status.
  type :: probe_result_t
    logical :: ok = .false.
    integer(i8) :: worksize = 0_i8, heap = 0_i8, xtdesc = 0_i8, context = 0_i8
    logical :: used_cufftmp = .false.
  end type probe_result_t

  type :: cuda_memcheck_device_t
    private
    !> The probe whose context and heap the estimate uses.
    type(probe_result_t), public :: probe
    !> The full-grid probe returned a failure status, and that status. The
    !> small-grid result then stays in place.
    logical, public :: full_failed = .false.
    integer, public :: full_failed_stat = 0
    !> The memory the FFT plans of the full grid need, from the query-only
    !> mode (see query_plan_need). For cuFFTMp it is the predicted heap plus
    !> the xtdesc buffer rounded up to 0.5 GiB.
    integer(i8), public :: plan_need_bytes = 0_i8
    !> The query-only mode used cuFFTMp, so the two sizes below are valid.
    logical :: plan_need_cufftmp = .false.
    !> The predicted heap without the xtdesc buffer, the same definition as
    !> probe%heap, and the xtdesc buffer apart, as probe%xtdesc.
    integer(i8) :: plan_heap_bytes = 0_i8, plan_xtdesc_bytes = 0_i8
    !> The estimate takes the heap and the xtdesc buffer from the two sizes
    !> above instead of the probe that ran (see use_predicted_heap).
    logical, public :: heap_predicted = .false.
    !> Card memory used and free right after cudaSetDevice.
    integer(i8) :: start_used_bytes = 0_i8, start_free_bytes = 0_i8
    !> The full-grid probe was tried, and it returned a result.
    logical :: full_tried = .false., full_ok = .false.
    !> Whether the real build used cuFFTMp, once it is known.
    logical :: use_cufftmp = .false., use_cufftmp_known = .false.
    type(cuda_allocator_t) :: cuda_allocator
    type(cuda_backend_t) :: cuda_backend
  contains
    procedure :: init
    procedure :: mem_info
    procedure :: context_bytes
    procedure :: sz => pencil_size
    procedure :: overhead_floor_bytes
    procedure :: overhead_query_bytes
    procedure :: probe_ran
    procedure :: make_backend
    procedure :: ensure_probe
    procedure :: query_plan_need
    procedure :: use_predicted_heap
    procedure :: heap_bytes
    procedure :: capture_solver_info
    procedure :: footprint_bytes
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
    self%start_free_bytes = int(free_b, i8)
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

  function pencil_size(self) result(n)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer :: n

    n = SZ
  end function pencil_size

  function overhead_floor_bytes(self, num_gpus, bc_is_110) result(nbytes8)
    !! The formula overhead at num_gpus GPUs: the context plus, unless the BC
    !! is 110 (plain cuFFT), the heap over num_gpus, both from the probe that
    !! ran.
    class(cuda_memcheck_device_t), intent(in) :: self
    integer, intent(in) :: num_gpus
    logical, intent(in) :: bc_is_110
    integer(i8) :: nbytes8

    nbytes8 = self%probe%context
    ! The worksize of a cuFFTMp plan is carved out of the heap.
    if (.not. bc_is_110) &
      nbytes8 = nbytes8 + &
                int(real(self%heap_bytes(), dp)/real(num_gpus, dp), i8)
    ! The predicted xtdesc buffer counts at num_gpus=1 only, as in
    ! overhead_query_bytes.
    if (self%heap_predicted .and. num_gpus == 1) &
      nbytes8 = nbytes8 + self%plan_xtdesc_bytes
  end function overhead_floor_bytes

  subroutine ensure_probe(self, bc_is_100, bc_is_110, cdims, try_full)
    !! Runs the FFT probe in this process. The probe on a 32^3 grid runs
    !! first, in every mode, and gives the context and the heap. When
    !! try_full is set, the probe on the cell dimensions cdims runs once
    !! after it. A full-grid probe that returns a failure status leaves the
    !! small-grid result in place and sets full_failed. A plan that does not
    !! fit aborts inside the library, so the caller sets try_full only when
    !! the plan fits. Stops the program when the small-grid probe fails.
    !! Each probe measures its own context just before its query.
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    logical, intent(in) :: try_full
    type(probe_result_t) :: res
    integer :: stat
    character(len=128) :: msg

    if (.not. self%probe%ok) then
      call run_probe(self, bc_is_100, bc_is_110, [32, 32, 32], res, stat)
      if (.not. res%ok) then
        write (msg, '(a,i0,a)') 'x3d2-memcheck: the calibration FFT probe &
          &returned status ', stat, '.'
        error stop trim(msg)
      end if
      self%probe = res
    end if
    if (try_full .and. .not. self%full_tried) then
      self%full_tried = .true.
      call run_probe(self, bc_is_100, bc_is_110, cdims, res, stat)
      self%full_ok = res%ok
      if (res%ok) then
        self%probe = res
      else
        self%full_failed = .true.
        self%full_failed_stat = stat
      end if
    end if
  end subroutine ensure_probe

  subroutine query_plan_need(self, bc_is_100, bc_is_110, cdims, &
                             xtdesc_bytes, stat)
    !! Runs the query-only mode of the FFT workspace query on the cell
    !! dimensions cdims. It allocates no workspace and no xtdesc buffer, so
    !! the plan cannot abort for lack of memory. xtdesc_bytes is the size of
    !! the xtdesc buffer, which the caller knows from the spectral slab.
    !! Sets plan_need_bytes when stat is zero.
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    integer(i8), intent(in) :: xtdesc_bytes
    integer, intent(out) :: stat
    integer(i8) :: worksize, heap, xtdesc
    logical :: used_cufftmp

    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   .true., worksize, heap, xtdesc, &
                                   used_cufftmp, stat, query_only=.true.)
    if (stat /= 0) return
    ! A plain cuFFT plan holds its workspace in the plan memory.
    self%plan_need_bytes = heap
    self%plan_need_cufftmp = used_cufftmp
    if (used_cufftmp) then
      self%plan_heap_bytes = predicted_fft_heap_bytes(heap, worksize)
      self%plan_xtdesc_bytes = round_up_bytes(xtdesc_bytes, &
                                              fft_xtdesc_step_bytes)
      self%plan_need_bytes = self%plan_heap_bytes + self%plan_xtdesc_bytes
    end if
  end subroutine query_plan_need

  subroutine use_predicted_heap(self)
    !! Makes the estimate take the heap predicted by query_plan_need, for the
    !! case that the full-grid probe does not run. A plain cuFFT plan has no
    !! heap, so it stays as it is.
    class(cuda_memcheck_device_t), intent(inout) :: self

    self%heap_predicted = self%plan_need_cufftmp
  end subroutine use_predicted_heap

  function heap_bytes(self) result(nbytes8)
    !! The NVSHMEM heap the estimate uses: the predicted heap when
    !! heap_predicted is set, else the heap of the probe that ran. Both are
    !! without the xtdesc buffer.
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8) :: nbytes8

    nbytes8 = self%probe%heap
    if (self%heap_predicted) nbytes8 = self%plan_heap_bytes
  end function heap_bytes

  subroutine run_probe(self, bc_is_100, bc_is_110, cdims, res, stat)
    !! The single-rank FFT probe on the cell dimensions cdims. The context
    !! is the card memory used just before the query. For the full-grid probe
    !! that includes what the small-grid probe's libraries keep allocated
    !! after cufftDestroy.
    class(cuda_memcheck_device_t), intent(in) :: self
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    type(probe_result_t), intent(out) :: res
    integer, intent(out) :: stat
    integer(i8) :: total_b, free_b

    call self%mem_info(total_b, free_b)
    call fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, .true., &
                                   .true., res%worksize, res%heap, &
                                   res%xtdesc, res%used_cufftmp, stat)
    res%ok = stat == 0
    if (res%ok) res%context = total_b - free_b
  end subroutine run_probe

  function overhead_query_bytes(self, num_gpus) result(nbytes8)
    !! The overhead at num_gpus GPUs from the full-grid probe, once probe_ran.
    class(cuda_memcheck_device_t), intent(in) :: self
    integer, intent(in) :: num_gpus
    integer(i8) :: nbytes8

    nbytes8 = self%probe%context
    if (self%probe%used_cufftmp) then
      ! The worksize is part of the heap. The xtdesc buffer costs nothing
      ! extra at num_gpus>1.
      nbytes8 = nbytes8 + &
                int(real(self%probe%heap, dp)/real(num_gpus, dp), i8) + &
                merge(self%probe%xtdesc, 0_i8, num_gpus == 1)
    else
      ! Plain cuFFT: the 110 BC or a cuFFTMp fallback, num_gpus=1 only.
      nbytes8 = nbytes8 + self%probe%worksize
    end if
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

  function footprint_bytes(self, free_end_bytes) result(nbytes8)
    !! The card memory held at the end of the real build: the context at the
    !! start plus the fall of the free memory since then.
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8), intent(in) :: free_end_bytes
    integer(i8) :: nbytes8

    nbytes8 = self%start_used_bytes + (self%start_free_bytes - free_end_bytes)
  end function footprint_bytes

  function fft_path_known(self) result(known)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: known

    known = self%use_cufftmp_known .or. self%full_ok
  end function fft_path_known

  function probe_ran(self) result(ran)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: ran

    ran = self%full_ok
  end function probe_ran

  function uses_distributed_fft(self) result(distributed)
    class(cuda_memcheck_device_t), intent(in) :: self
    logical :: distributed

    if (self%use_cufftmp_known) then
      distributed = self%use_cufftmp
    else
      distributed = self%probe%used_cufftmp
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

  subroutine fft_workspace_bytes_query(bc_is_100, bc_is_110, cdims, &
                                       try_cufftmp, is_root, &
                                       worksize_bytes, heap_bytes, &
                                       xtdesc_bytes, used_cufftmp, stat, &
                                       query_only)
    !! Creates and destroys the forward and backward FFT plans the way the
    !! solver does, for the cell dimensions cdims, and for cuFFTMp also the
    !! xtdesc buffer, each bracketed with cudaMemGetInfo. worksize_bytes is
    !! the plan workspace (used only on the plain cuFFT path, in cuFFTMp it
    !! is part of the heap), heap_bytes the memory the plans took,
    !! xtdesc_bytes the xtdesc buffer. The 110 BC never uses cuFFTMp. stat is
    !! the first non-zero cuFFT/CUDA code, or 0, and the other outputs are
    !! then not valid. A plan that aborts inside the library never returns.
    !! With query_only set, a cuFFTMp plan is made with auto-allocation off.
    !! heap_bytes is then the plan memory without the workspace, worksize_bytes
    !! the workspace of both plans (not allocated), xtdesc_bytes is zero and no
    !! xtdesc buffer is allocated. A plain cuFFT plan is not changed by it.
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    logical, intent(in) :: try_cufftmp
    logical, intent(in) :: is_root
    integer(i8), intent(out) :: worksize_bytes, heap_bytes, xtdesc_bytes
    logical, intent(out) :: used_cufftmp
    integer, intent(out) :: stat
    logical, optional, intent(in) :: query_only

    integer :: plan_fw, plan_bw, fft_n1, fft_n2, fft_n3
    integer :: fw_plan_type, bw_plan_type, ierr
    integer(int_ptr_kind()) :: worksize_fw, worksize_bw
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b
    type(cudaLibXtDesc), pointer :: xtdesc
    logical :: fw_was_cufftmp, no_auto_alloc

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
    no_auto_alloc = .false.
    if (present(query_only)) no_auto_alloc = query_only

    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before plans)')
    call create_fft_plan(plan_fw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         fw_plan_type, is_root, 'memcheck probe fwd', &
                         worksize_fw, stat, no_auto_alloc)
    if (stat /= 0) return
    fw_was_cufftmp = used_cufftmp
    ! A fallback to plain cuFFT sets used_cufftmp to false. A forward plan
    ! built on cuFFTMp is then rebuilt on plain cuFFT, as init() does.
    call create_fft_plan(plan_bw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         bw_plan_type, is_root, 'memcheck probe bwd', &
                         worksize_bw, stat, no_auto_alloc)
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

    if (used_cufftmp .and. (.not. no_auto_alloc)) then
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
