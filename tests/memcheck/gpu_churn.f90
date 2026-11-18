program gpu_churn
  !! Allocates and frees a device array in a loop, to change the memory of
  !! this process on the card. It checks that x3d2-memcheck --build skips its
  !! CHECK line on a shared card. Usage: gpu_churn [seconds] [gib].
  use cudafor
  use iso_fortran_env, only: int64
  use iso_c_binding, only: c_int
  implicit none

  interface
    integer(c_int) function usleep(usec) bind(C, name='usleep')
      import :: c_int
      integer(c_int), value :: usec
    end function usleep
  end interface

  real, parameter :: GIB = 1024.**3
  real :: seconds, gib_alloc
  integer :: cycle_count, stat, ierr
  integer(int64) :: nbytes, t_start, t_now, t_rate
  type(c_devptr) :: buf
  character(len=32) :: arg

  seconds = 60.
  gib_alloc = 1.5
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read (arg, *) seconds
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read (arg, *) gib_alloc
  end if
  nbytes = int(gib_alloc*GIB, int64)

  ierr = cudaSetDevice(0)
  call system_clock(t_start, t_rate)
  cycle_count = 0
  do
    call system_clock(t_now)
    if (real(t_now - t_start)/real(t_rate) >= seconds) exit
    cycle_count = cycle_count + 1
    ! cudaMalloc and cudaFree return the memory to the card. An allocate of
    ! a device array does not: nvfortran keeps it in a pool.
    stat = cudaMalloc(buf, nbytes)
    if (stat /= 0) then
      print '(a,i0,a,i0)', 'gpu_churn cycle ', cycle_count, &
        ': allocate failed, stat ', stat
    else
      ierr = cudaMemset(buf, 1, nbytes)
      ierr = cudaDeviceSynchronize()
      print '(a,i0,a,f0.2,a)', 'gpu_churn cycle ', cycle_count, &
        ': holding ', gib_alloc, ' GiB'
    end if
    ierr = usleep(500000_c_int)
    if (stat == 0) ierr = cudaFree(buf)
    print '(a,i0,a)', 'gpu_churn cycle ', cycle_count, ': freed'
    ierr = usleep(500000_c_int)
  end do

end program gpu_churn
