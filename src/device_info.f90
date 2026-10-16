module m_omptgt_device_info
  !! Peak memory clock and bus width of an OpenMP target device, read from
  !! the vendor runtime.
  !!
  !! OpenMP offers no query for either, so this asks the CUDA driver on
  !! NVIDIA hardware and HIP on AMD hardware, the same attributes the CUDA
  !! backend reads. The vendor library is opened at run time with dlopen
  !! instead of being linked: no compiler's offload runtime links it for
  !! us, the backend builds with NVHPC, Flang and Cray alike, and a machine
  !! without the library reports the figures as unavailable rather than
  !! failing to link or start.
  !!
  !! The OpenMP device number is passed as the vendor device ordinal. Both
  !! enumerate the devices visible to the process in the same order.
  use iso_c_binding, only: c_ptr, c_funptr, c_int, c_char, c_null_char, &
                           c_associated, c_f_procpointer

  implicit none

  private
  public :: query_device_memory

  integer(c_int), parameter :: RTLD_NOW = 2

  ! CU_DEVICE_ATTRIBUTE_* in cuda.h (CUDA 13 still provides both)
  integer(c_int), parameter :: CU_MEMORY_CLOCK_RATE = 36
  integer(c_int), parameter :: CU_GLOBAL_MEMORY_BUS_WIDTH = 37
  ! hipDeviceAttribute_t in hip_runtime_api.h (unchanged from ROCm 5.7 on)
  integer(c_int), parameter :: HIP_MEMORY_BUS_WIDTH = 59
  integer(c_int), parameter :: HIP_MEMORY_CLOCK_RATE = 60

  interface
    type(c_ptr) function dlopen(filename, flag) bind(c, name='dlopen')
      import :: c_ptr, c_char, c_int
      character(kind=c_char), dimension(*), intent(in) :: filename
      integer(c_int), value :: flag
    end function dlopen

    type(c_funptr) function dlsym(handle, symbol) bind(c, name='dlsym')
      import :: c_ptr, c_funptr, c_char
      type(c_ptr), value :: handle
      character(kind=c_char), dimension(*), intent(in) :: symbol
    end function dlsym

    integer(c_int) function dlclose(handle) bind(c, name='dlclose')
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
    end function dlclose
  end interface

  abstract interface
    integer(c_int) function cu_init_t(flags) bind(c)
      import :: c_int
      integer(c_int), value :: flags
    end function cu_init_t

    integer(c_int) function cu_device_get_t(device, ordinal) bind(c)
      import :: c_int
      integer(c_int), intent(out) :: device
      integer(c_int), value :: ordinal
    end function cu_device_get_t

    integer(c_int) function get_attribute_t(value, attribute, device) &
      bind(c)
      !! cuDeviceGetAttribute and hipDeviceGetAttribute share this shape
      import :: c_int
      integer(c_int), intent(out) :: value
      integer(c_int), value :: attribute, device
    end function get_attribute_t
  end interface

contains

  subroutine query_device_memory(device, mem_clock_khz, bus_width_bits, &
                                 available)
    !! Memory clock (kHz) and bus width (bits) of OpenMP device `device`.
    !! The vendor this build targets is tried first, then the other one, as
    !! the build only guesses the vendor for some compilers.
    integer, intent(in) :: device
    integer, intent(out) :: mem_clock_khz, bus_width_bits
    logical, intent(out) :: available

#if defined(OMP_TGT_AMD)
    call query_hip(device, mem_clock_khz, bus_width_bits, available)
    if (.not. available) then
      call query_cuda(device, mem_clock_khz, bus_width_bits, available)
    end if
#else
    call query_cuda(device, mem_clock_khz, bus_width_bits, available)
    if (.not. available) then
      call query_hip(device, mem_clock_khz, bus_width_bits, available)
    end if
#endif

  end subroutine query_device_memory

  subroutine query_cuda(device, mem_clock_khz, bus_width_bits, available)
    integer, intent(in) :: device
    integer, intent(out) :: mem_clock_khz, bus_width_bits
    logical, intent(out) :: available

    type(c_ptr) :: lib
    type(c_funptr) :: fp_init, fp_get, fp_attr
    procedure(cu_init_t), pointer :: cu_init
    procedure(cu_device_get_t), pointer :: cu_device_get
    procedure(get_attribute_t), pointer :: get_attribute
    integer(c_int) :: dev, clock, width, ierr

    mem_clock_khz = 0
    bus_width_bits = 0
    available = .false.

    lib = dlopen('libcuda.so.1'//c_null_char, RTLD_NOW)
    if (.not. c_associated(lib)) return

    fp_init = dlsym(lib, 'cuInit'//c_null_char)
    fp_get = dlsym(lib, 'cuDeviceGet'//c_null_char)
    fp_attr = dlsym(lib, 'cuDeviceGetAttribute'//c_null_char)
    if (c_associated(fp_init) .and. c_associated(fp_get) .and. &
        c_associated(fp_attr)) then
      call c_f_procpointer(fp_init, cu_init)
      call c_f_procpointer(fp_get, cu_device_get)
      call c_f_procpointer(fp_attr, get_attribute)
      ! The offload runtime has initialised the driver already; cuInit is
      ! then a no-op, but required before any other driver call.
      if (cu_init(0_c_int) == 0) then
        if (cu_device_get(dev, int(device, c_int)) == 0) then
          if (get_attribute(clock, CU_MEMORY_CLOCK_RATE, dev) == 0 .and. &
              get_attribute(width, CU_GLOBAL_MEMORY_BUS_WIDTH, dev) == 0) &
            then
            mem_clock_khz = clock
            bus_width_bits = width
            available = clock > 0 .and. width > 0
          end if
        end if
      end if
    end if
    ierr = dlclose(lib)

  end subroutine query_cuda

  subroutine query_hip(device, mem_clock_khz, bus_width_bits, available)
    integer, intent(in) :: device
    integer, intent(out) :: mem_clock_khz, bus_width_bits
    logical, intent(out) :: available

    character(len=*), parameter :: names(3) = [character(len=17) :: &
                                               'libamdhip64.so', &
                                               'libamdhip64.so.6', &
                                               'libamdhip64.so.5']
    type(c_ptr) :: lib
    type(c_funptr) :: fp
    procedure(get_attribute_t), pointer :: get_attribute
    integer(c_int) :: clock, width, ierr
    integer :: i

    mem_clock_khz = 0
    bus_width_bits = 0
    available = .false.

    do i = 1, size(names)
      lib = dlopen(trim(names(i))//c_null_char, RTLD_NOW)
      if (c_associated(lib)) exit
    end do
    if (.not. c_associated(lib)) return

    fp = dlsym(lib, 'hipDeviceGetAttribute'//c_null_char)
    if (c_associated(fp)) then
      call c_f_procpointer(fp, get_attribute)
      if (get_attribute(clock, HIP_MEMORY_CLOCK_RATE, &
                        int(device, c_int)) == 0 .and. &
          get_attribute(width, HIP_MEMORY_BUS_WIDTH, &
                        int(device, c_int)) == 0) then
        mem_clock_khz = clock
        bus_width_bits = width
        available = clock > 0 .and. width > 0
      end if
    end if
    ierr = dlclose(lib)

  end subroutine query_hip

end module m_omptgt_device_info
