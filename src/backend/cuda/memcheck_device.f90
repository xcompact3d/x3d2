module m_cuda_memcheck_device
  !! CUDA implementation of the x3d2-memcheck device seam: the only
  !! memcheck file that uses cudafor or the m_cuda_* modules.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice
  use m_memcheck_device, only: memcheck_device_t
  implicit none

  private
  public :: cuda_memcheck_device_t

  type, extends(memcheck_device_t) :: cuda_memcheck_device_t
  contains
    procedure :: init
  end type cuda_memcheck_device_t

contains

  subroutine init(self, ok, msg)
    class(cuda_memcheck_device_t), intent(inout) :: self
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    integer :: ierr, ndevs

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
    ok = .true.
  end subroutine init

end module m_cuda_memcheck_device
