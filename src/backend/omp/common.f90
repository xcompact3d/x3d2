module m_omp_common
  implicit none

  !> Number of points a pencil group stacks together in the leading
  !! dimension, sized so that a group's working set stays in cache.
  integer, parameter :: SZ = 16

end module m_omp_common
