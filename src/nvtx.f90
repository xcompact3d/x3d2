module m_nvtx
!! NVTX ranges for profiling I/O with Nsight Systems.
!!
!! Compiled in for CUDA builds with ADIOS2 (X3D2_NVTX), where they label both
!! the GPU-aware and the host-staged write paths; a no-op otherwise. Set
!! X3D2_NVTX=0 at runtime to turn them off.
  use iso_c_binding, only: c_char, c_null_char, c_ptr, c_loc, c_int

  implicit none

  private
  public :: nvtx_push, nvtx_pop, nvtx_enabled

  logical, save :: enabled = .true.
  logical, save :: initialised = .false.

#ifdef X3D2_NVTX
  interface
    function nvtx_range_push_a(name) bind(C, name='nvtxRangePushA') &
      result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: name
      integer(c_int) :: status
    end function nvtx_range_push_a

    function nvtx_range_pop() bind(C, name='nvtxRangePop') result(status)
      import :: c_int
      integer(c_int) :: status
    end function nvtx_range_pop
  end interface
#endif

contains

  logical function nvtx_enabled()
    !! Whether ranges are emitted: X3D2_NVTX (default on), read on first use.
    character(len=16) :: raw_value
    integer :: status, i

    if (.not. initialised) then
      initialised = .true.
      call get_environment_variable("X3D2_NVTX", raw_value, status=status)
      if (status == 0) then
        do i = 1, len(raw_value)
          if (raw_value(i:i) >= "A" .and. raw_value(i:i) <= "Z") then
            raw_value(i:i) = achar(iachar(raw_value(i:i)) + 32)
          end if
        end do
        select case (trim(adjustl(raw_value)))
        case ("0", "false", "no", "off")
          enabled = .false.
        end select
      end if
    end if
#ifdef X3D2_NVTX
    nvtx_enabled = enabled
#else
    nvtx_enabled = .false.
#endif
  end function nvtx_enabled

  subroutine nvtx_push(range_name)
    character(len=*), intent(in) :: range_name
#ifdef X3D2_NVTX
    integer :: nvtx_status
    character(kind=c_char), allocatable, target :: c_name(:)
    integer :: i, n
    if (nvtx_enabled()) then
      n = len_trim(range_name)
      allocate (c_name(n + 1))
      do i = 1, n
        c_name(i) = achar(iachar(range_name(i:i)), kind=c_char)
      end do
      c_name(n + 1) = c_null_char
      nvtx_status = nvtx_range_push_a(c_loc(c_name(1)))
    end if
#endif
  end subroutine nvtx_push

  subroutine nvtx_pop()
#ifdef X3D2_NVTX
    integer :: nvtx_status
    if (nvtx_enabled()) nvtx_status = nvtx_range_pop()
#endif
  end subroutine nvtx_pop

end module m_nvtx
