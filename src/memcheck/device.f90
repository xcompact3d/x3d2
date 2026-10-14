module m_memcheck_device
  !! The backend seam of x3d2-memcheck: everything the generic estimate,
  !! report and build modules need from a compute backend, so that only one
  !! implementation (src/backend/cuda/memcheck_device.f90 for CUDA) touches
  !! backend-specific libraries.
  use m_common, only: i8
  implicit none

  type, abstract :: memcheck_device_t
  contains
    !> Select the device. ok=.false. and the line to print in msg if there
    !> is no usable one.
    procedure(init_iface), deferred :: init
    !> Total and currently free device memory.
    procedure(mem_info_iface), deferred :: mem_info
  end type memcheck_device_t

  abstract interface
    subroutine init_iface(self, ok, msg)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(inout) :: self
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
    end subroutine init_iface

    subroutine mem_info_iface(self, total_bytes, free_bytes)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(in) :: self
      integer(i8), intent(out) :: total_bytes, free_bytes
    end subroutine mem_info_iface
  end interface

end module m_memcheck_device
