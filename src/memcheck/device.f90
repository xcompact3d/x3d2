module m_memcheck_device
  !! The backend seam of x3d2-memcheck: everything the generic estimate,
  !! report and build modules need from a compute backend, so that only one
  !! implementation (src/backend/cuda/memcheck_device.f90 for CUDA) touches
  !! backend-specific libraries.
  implicit none

  type, abstract :: memcheck_device_t
  contains
    !> Select the device. ok=.false. and the line to print in msg if there
    !> is no usable one.
    procedure(init_iface), deferred :: init
  end type memcheck_device_t

  abstract interface
    subroutine init_iface(self, ok, msg)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(inout) :: self
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
    end subroutine init_iface
  end interface

end module m_memcheck_device
