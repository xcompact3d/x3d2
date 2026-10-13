module m_omptgt_common
  implicit none

  !> Number of points a pencil group stacks together in the leading
  !! dimension, set to the width the target hardware issues a wavefront at so
  !! that one group fills one wavefront. The last branch covers a build that
  !! selects no vendor.
#if defined(OMP_TGT_AMD)
  integer, parameter :: SZ = 64
#elif defined(OMP_TGT_NVIDIA)
  integer, parameter :: SZ = 32
#else
  integer, parameter :: SZ = 16
#endif

end module m_omptgt_common
