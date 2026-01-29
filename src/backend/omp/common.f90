module m_omp_common
  !! Common constants for OpenMP backend implementation.
  !!
  !! Defines compile-time constants used throughout the OMP backend
  !! for performance tuning and buffer sizing.
  !!
  !! **SZ (pencil size):** Maximum pencil dimension for data reordering
  !! operations. Set to 16 for optimal cache utilisation and vectorisation
  !! on typical CPU architectures. Larger values may improve performance
  !! for very large problems but increase memory overhead.
  implicit none

  integer, parameter :: SZ = 16 !! Maximum pencil size for reordering buffers

end module m_omp_common
