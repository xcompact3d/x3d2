module m_cuda_common
  !! Common constants for CUDA backend.
  !!
  !! CUDA GPUs execute threads in groups of 32 called warps. Setting the
  !! pencil size to 32 ensures coalesced memory access patterns, where all
  !! threads in a warp access consecutive memory locations simultaneously.
  !! This is critical for GPU memory bandwidth efficiency.
  !!
  !! **Performance impact:** Matching the hardware warp size eliminates
  !! divergence and maximises memory throughput, typically improving
  !! performance by 2-3x compared to non-coalesced access.
  implicit none

  integer, parameter :: SZ = 32  !! Pencil size matching GPU warp width

end module m_cuda_common
