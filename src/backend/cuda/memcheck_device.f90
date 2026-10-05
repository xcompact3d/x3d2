module m_cuda_memcheck_device
  !! CUDA implementation of the x3d2-memcheck device seam: the only
  !! memcheck file that uses cudafor or the m_cuda_* modules.
  use cudafor, only: cudaGetDeviceCount, cudaSetDevice, cudaMemGetInfo, &
                     cuda_count_kind
  use m_common, only: i8
  use m_cuda_memory_estimate, only: check_status
  use m_memcheck_device, only: memcheck_device_t
  implicit none

  private
  public :: cuda_memcheck_device_t

  type, extends(memcheck_device_t) :: cuda_memcheck_device_t
  contains
    procedure :: init
    procedure :: mem_info
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

  subroutine mem_info(self, total_bytes, free_bytes)
    class(cuda_memcheck_device_t), intent(in) :: self
    integer(i8), intent(out) :: total_bytes, free_bytes
    integer :: ierr
    integer(kind=cuda_count_kind) :: free_b, total_b

    ierr = cudaMemGetInfo(free_b, total_b)
    call check_status(ierr, 'cudaMemGetInfo')
    total_bytes = int(total_b, i8)
    free_bytes = int(free_b, i8)
  end subroutine mem_info

end module m_cuda_memcheck_device
