module m_omp_common
  implicit none

#if defined(OMP_TGT_AMD)
  integer, parameter :: SZ = 64
#elif defined(OMP_TGT_NVIDIA)
  integer, parameter :: SZ = 32
#else
  integer, parameter :: SZ = 16
#endif

end module m_omp_common
