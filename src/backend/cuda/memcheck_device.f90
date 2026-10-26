module m_cuda_memcheck_device
  !! The CUDA device of x3d2-memcheck: the FFT plan probe and the memory
  !! readings.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice, cudaMemGetInfo, &
                     cuda_count_kind, int_ptr_kind
  use cufft, only: cufftDestroy
  use cufftXt, only: cufftXtMalloc, cufftXtFree, CUFFT_XT_FORMAT_INPLACE, &
                     cudaLibXtDesc
  use m_common, only: dp, i8
  use m_logger, only: abort_run
  use iso_c_binding, only: c_int, c_long_long, c_ptr, c_char
  use m_memory_estimate, only: predicted_fft_heap_bytes, round_up_bytes, &
                               fft_xtdesc_step_bytes
  use m_cuda_common, only: SZ
  use m_mesh, only: FFT_BC_000, FFT_BC_100, FFT_BC_110
  use m_base_case, only: base_case_t
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t, create_fft_plan_pair, &
                                cuda_fft_plan_setup
  implicit none

  private
  public :: cuda_memcheck_device_t, fft_workspace_bytes_query, &
            foreign_start_ok, foreign_start_nproc, foreign_start_bytes, &
            foreign_worst, foreign_nproc_last, foreign_valid

  !> The v3 layout of nvmlProcessInfo_t.
  type, bind(C) :: nvml_proc_t
    integer(c_int) :: pid
    integer(c_long_long) :: used_bytes
    integer(c_int) :: gpu_instance_id, compute_instance_id
  end type nvml_proc_t
  integer, parameter :: nvml_max_procs = 256
  !> NVML_VALUE_NOT_AVAILABLE: all bits set.
  integer(c_long_long), parameter :: nvml_not_available = -1_c_long_long

  !> The other processes on the card at the start (false: NVML gave no usable
  !> list), the largest change from the start over all samples so far, and the
  !> process count of the last sample.
  logical, save, protected :: foreign_start_ok = .false.
  integer, save, protected :: foreign_start_nproc = 0, foreign_nproc_last = 0
  integer(i8), save, protected :: foreign_start_bytes = 0_i8, &
                                  foreign_worst = 0_i8
  !> False after a failed sample, which stops the sampling.
  logical, save, protected :: foreign_valid = .false.

  interface
    integer(c_int) function nvmlInit_v2() bind(C, name='nvmlInit_v2')
      import :: c_int
    end function nvmlInit_v2

    integer(c_int) function nvmlDeviceGetHandleByPciBusId_v2(busid, device) &
      bind(C, name='nvmlDeviceGetHandleByPciBusId_v2')
      import :: c_int, c_char, c_ptr
      character(kind=c_char), dimension(*), intent(in) :: busid
      type(c_ptr), intent(out) :: device
    end function nvmlDeviceGetHandleByPciBusId_v2

    !> NVML_ERROR_INSUFFICIENT_SIZE (7) is a failure like any other.
    integer(c_int) function nvmlDeviceGetComputeRunningProcesses_v3( &
      device, count, infos) &
      bind(C, name='nvmlDeviceGetComputeRunningProcesses_v3')
      import :: c_int, c_ptr, nvml_proc_t, nvml_max_procs
      type(c_ptr), value :: device
      integer(c_int), intent(inout) :: count
      type(nvml_proc_t), intent(out) :: infos(nvml_max_procs)
    end function nvmlDeviceGetComputeRunningProcesses_v3

    !> The CUDA runtime call: cudafor has no Fortran binding for it.
    integer(c_int) function cudaDeviceGetPCIBusId(busid, len, device) &
      bind(C, name='cudaDeviceGetPCIBusId')
      import :: c_int, c_char
      character(kind=c_char), dimension(*), intent(out) :: busid
      integer(c_int), value :: len, device
    end function cudaDeviceGetPCIBusId
  end interface

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
  contains
    procedure :: init
    procedure :: mem_info
    procedure :: context_bytes
    procedure :: sz => pencil_size
    procedure :: overhead_floor_bytes
    procedure :: overhead_query_bytes
    procedure :: probe_ran
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
    call foreign_gpu_bytes(foreign_start_ok, foreign_start_nproc, &
                           foreign_start_bytes)
    foreign_valid = foreign_start_ok
    foreign_worst = 0_i8
    foreign_nproc_last = foreign_start_nproc
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
    call foreign_note()
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
        call abort_run(trim(msg), collective=.true.)
      end if
      call require_cufftmp(bc_is_110, res%used_cufftmp, 'calibration')
      self%probe = res
    end if
    if (try_full .and. .not. self%full_tried) then
      self%full_tried = .true.
      call run_probe(self, bc_is_100, bc_is_110, cdims, res, stat)
      self%full_ok = res%ok
      if (res%ok) then
        call require_cufftmp(bc_is_110, res%used_cufftmp, 'full-grid')
        self%probe = res
      else
        self%full_failed = .true.
        self%full_failed_stat = stat
      end if
    end if
  end subroutine ensure_probe

  subroutine require_cufftmp(bc_is_110, used_cufftmp, what)
    !! Aborts x3d2-memcheck when an FFT plan that asked for cuFFTMp got plain
    !! cuFFT. fft_workspace_bytes_query asks for cuFFTMp for every BC except
    !! 110, so only a 110 plan is allowed to be plain cuFFT. This is the
    !! 110 exception to the solver rule: the solver uses cuFFTMp for 110 on
    !! nproc > 1, but x3d2-memcheck keeps 110 single GPU until the multi-rank
    !! buffers of #417 are in the estimate.
    logical, intent(in) :: bc_is_110, used_cufftmp
    character(*), intent(in) :: what

    if (bc_is_110 .or. used_cufftmp) return
    call abort_run('x3d2-memcheck needs cuFFTMp, but the '//what// &
                   ' FFT plan fell back to plain cuFFT.', collective=.true.)
  end subroutine require_cufftmp

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
    call require_cufftmp(bc_is_110, used_cufftmp, 'plan query')
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
      ! Plain cuFFT: the 110 BC, num_gpus=1 only.
      nbytes8 = nbytes8 + self%probe%worksize
    end if
  end function overhead_query_bytes

  subroutine capture_solver_info(self, flow_case)
    class(cuda_memcheck_device_t), intent(inout) :: self
    class(base_case_t), intent(in) :: flow_case

    select type (pf => flow_case%solver%backend%poisson_fft)
    type is (cuda_poisson_fft_t)
      call require_cufftmp(pf%is_110_case, pf%use_cufftmp, 'solver')
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

  subroutine foreign_gpu_bytes(ok, nproc, nbytes)
    !! The processes other than this one on the device, and the memory they
    !! hold, from NVML. ok is false when a call fails, a memory value is not
    !! available or our own PID is not in the list (a container shows host
    !! PIDs). NVML picks the card by PCI bus id, not by CUDA_VISIBLE_DEVICES.
    logical, intent(out) :: ok
    integer, intent(out) :: nproc
    integer(i8), intent(out) :: nbytes
    character(kind=c_char, len=32) :: busid
    type(c_ptr) :: handle
    type(nvml_proc_t) :: procs(nvml_max_procs)
    integer(c_int) :: count
    integer :: i, getpid
    logical :: found_self

    ok = .false.
    nproc = 0
    nbytes = 0_i8
    count = nvml_max_procs
    ! The driver cleans NVML up at process exit, so there is no shutdown.
    if (nvmlInit_v2() /= 0) return
    ! Device 0 is the one init set.
    if (cudaDeviceGetPCIBusId(busid, int(len(busid), c_int), 0_c_int) /= 0) &
      return
    if (nvmlDeviceGetHandleByPciBusId_v2(busid, handle) /= 0) return
    if (nvmlDeviceGetComputeRunningProcesses_v3(handle, count, procs) /= 0) &
      return
    found_self = .false.
    do i = 1, count
      if (procs(i)%used_bytes == nvml_not_available) return
      if (procs(i)%pid == getpid()) then
        found_self = .true.
      else
        nproc = nproc + 1
        nbytes = nbytes + int(procs(i)%used_bytes, i8)
      end if
    end do
    ok = found_self
  end subroutine foreign_gpu_bytes

  subroutine foreign_note()
    !! Samples the other processes next to a cudaMemGetInfo reading and keeps
    !! the largest change from the start.
    logical :: ok
    integer :: nproc
    integer(i8) :: nbytes

    if (.not. foreign_valid) return
    call foreign_gpu_bytes(ok, nproc, nbytes)
    foreign_valid = ok
    if (.not. ok) return
    foreign_nproc_last = nproc
    if (abs(nbytes - foreign_start_bytes) > abs(foreign_worst)) &
      foreign_worst = nbytes - foreign_start_bytes
  end subroutine foreign_note

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
    call abort_run(trim(msg), collective=.true.)
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

    integer :: plan_fw, plan_bw, fft_n(3), fft_n_loc(2), bc_code
    integer :: fw_plan_type, bw_plan_type, ierr
    integer(int_ptr_kind()) :: worksize_fw, worksize_bw
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b
    type(cudaLibXtDesc), pointer :: xtdesc
    logical :: no_auto_alloc

    ! The plan dimensions and types are the solver's, one rank.
    if (bc_is_100) then
      bc_code = FFT_BC_100
    else if (bc_is_110) then
      bc_code = FFT_BC_110
    else
      bc_code = FFT_BC_000
    end if
    call cuda_fft_plan_setup(bc_code, cdims, cdims, fft_n, fft_n_loc, &
                             fw_plan_type, bw_plan_type)

    used_cufftmp = try_cufftmp .and. (.not. bc_is_110)
    stat = 0
    no_auto_alloc = .false.
    if (present(query_only)) no_auto_alloc = query_only

    ierr = cudaMemGetInfo(free_before, total_b)
    call check_status(ierr, 'cudaMemGetInfo (before plans)')
    call foreign_note()
    call create_fft_plan_pair(plan_fw, plan_bw, used_cufftmp, fft_n, &
                              fw_plan_type, bw_plan_type, ' memcheck probe', &
                              worksize_fw, worksize_bw, stat, no_auto_alloc)
    if (stat /= 0) return
    ierr = cudaMemGetInfo(free_after, total_b)
    call check_status(ierr, 'cudaMemGetInfo (after plans)')
    call foreign_note()
    worksize_bytes = int(worksize_fw, i8) + int(worksize_bw, i8)
    heap_bytes = int(free_before, i8) - int(free_after, i8)

    if (used_cufftmp .and. (.not. no_auto_alloc)) then
      ierr = cudaMemGetInfo(free_before, total_b)
      call check_status(ierr, 'cudaMemGetInfo (before xtdesc)')
      call foreign_note()
      ierr = cufftXtMalloc(plan_fw, xtdesc, CUFFT_XT_FORMAT_INPLACE)
      call check_status(ierr, 'cufftXtMalloc', stat)
      if (stat /= 0) then
        ierr = cufftDestroy(plan_fw)
        ierr = cufftDestroy(plan_bw)
        return
      end if
      ierr = cudaMemGetInfo(free_after, total_b)
      call check_status(ierr, 'cudaMemGetInfo (after xtdesc)')
      call foreign_note()
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
