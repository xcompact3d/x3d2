module m_cuda_memory_estimate
  !! CUDA-only Tier 2 of the static memory estimator: create a throwaway
  !! FFT plan for the case's actual BC/dims/rank-count via the real
  !! solver's own create_fft_plan (src/backend/cuda/poisson_fft.f90),
  !! capture cuFFT/cuFFTMp's real reported workspace size, then destroy it
  !! - no mesh, no fields, no case object, no substep. This is what lets
  !! x3d2-memcheck report an exact FFT-workspace number for the one term
  !! Tier 1 (m_memory_estimate, all-backend, no-GPU) cannot get from array
  !! shapes alone.
  !!
  !! This lives under backend/cuda (CUDASRC, not SRC) rather than inside
  !! m_memory_estimate because it depends on m_cuda_poisson_fft, cudafor,
  !! and cufft, which are only compiled for the NVHPC/PGI backend;
  !! m_memory_estimate itself must stay buildable on every backend (it is
  !! in SRC, not CUDASRC).
  use cudafor, only: int_ptr_kind, cudaMemGetInfo, cuda_count_kind
  use cufft, only: cufftDestroy, CUFFT_D2Z, CUFFT_Z2D, CUFFT_R2C, CUFFT_C2R
  use cufftXt, only: cufftXtMalloc, cufftXtFree, CUFFT_XT_FORMAT_INPLACE, &
                     cudaLibXtDesc
  use m_common, only: dp, i8, is_sp
  use m_cuda_poisson_fft, only: create_fft_plan

  implicit none

  private
  public :: fft_workspace_bytes_query, context_floor_bytes

contains

  pure function context_floor_bytes(ng, uses_cufftmp) result(nbytes8)
    !! Grid-independent memory reserved before any field/FFT-data bytes:
    !! the CUDA context itself, plus (when cuFFTMp is used) an estimate of
    !! the NVSHMEM symmetric heap. Two different callers use this function
    !! for two different purposes:
    !!   - memcheck.f90's Tier 1 floor (before any GPU plan is queried) and
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
    !! memcheck.f90's own Tier 2 cuFFTMp estimate does NOT use this
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
    !!     the Tier 2 query's worksize must NOT be added on top of this
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
                                       xtdesc_bytes, used_cufftmp)
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
    !! each with its own auto-allocated workspace.
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
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    logical, intent(in) :: try_cufftmp
    logical, intent(in) :: is_root
    integer(i8), intent(out) :: worksize_bytes, heap_bytes, xtdesc_bytes
    logical, intent(out) :: used_cufftmp

    integer :: plan_fw, plan_bw, fft_n1, fft_n2, fft_n3
    integer :: fw_plan_type, bw_plan_type, ierr
    integer(int_ptr_kind()) :: worksize_fw, worksize_bw
    integer(kind=cuda_count_kind) :: free_before, free_after, total_b
    type(cudaLibXtDesc), pointer :: xtdesc

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

    ierr = cudaMemGetInfo(free_before, total_b)
    call create_fft_plan(plan_fw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         fw_plan_type, is_root, 'memcheck probe fwd', &
                         worksize_fw)
    ! If the forward plan fell back from cuFFTMp to plain cuFFT,
    ! used_cufftmp is now .false. - match it for the backward plan too, the
    ! same way init() rebuilds the forward plan on a backward-plan fallback
    ! (src/backend/cuda/poisson_fft.f90:424-428) rather than running one
    ! plan on cuFFTMp and the other on plain cuFFT.
    call create_fft_plan(plan_bw, used_cufftmp, fft_n1, fft_n2, fft_n3, &
                         bw_plan_type, is_root, 'memcheck probe bwd', &
                         worksize_bw)
    ierr = cudaMemGetInfo(free_after, total_b)
    worksize_bytes = int(worksize_fw, i8) + int(worksize_bw, i8)
    heap_bytes = int(free_before, i8) - int(free_after, i8)

    if (used_cufftmp) then
      ierr = cudaMemGetInfo(free_before, total_b)
      ierr = cufftXtMalloc(plan_fw, xtdesc, CUFFT_XT_FORMAT_INPLACE)
      ierr = cudaMemGetInfo(free_after, total_b)
      xtdesc_bytes = int(free_before, i8) - int(free_after, i8)
      ierr = cufftXtFree(xtdesc)
    else
      xtdesc_bytes = 0_i8
    end if

    ierr = cufftDestroy(plan_fw)
    ierr = cufftDestroy(plan_bw)
  end subroutine fft_workspace_bytes_query

end module m_cuda_memory_estimate
