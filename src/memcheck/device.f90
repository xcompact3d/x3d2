module m_memcheck_device
  !! The backend seam of x3d2-memcheck: everything the generic estimate,
  !! report and build modules need from a compute backend, so that only one
  !! implementation (src/backend/cuda/memcheck_device.f90 for CUDA) touches
  !! backend-specific libraries.
  implicit none

  type, abstract :: memcheck_device_t
  end type memcheck_device_t

end module m_memcheck_device
