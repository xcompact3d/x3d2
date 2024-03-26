  program test_cuda_reorder
    use iso_fortran_env, only: stderr => error_unit
    use cudafor

    use m_common, only: dp
    use m_cuda_common, only: SZ
    use m_cuda_kernels_reorder, only: reorder_x2y, reorder_x2z, reorder_y2x, &
      reorder_y2z, reorder_z2x, reorder_z2y, &
      reorder_c2x, reorder_x2c

    implicit none

    logical :: allpass = .true.
    real(dp), allocatable, dimension(:, :, :) :: u_i, u_o, u_temp, u_c
    real(dp), device, allocatable, dimension(:, :, :) :: u_i_d, u_o_d, &
      u_temp_d, u_c_d

    integer :: n_block, i, n_iters
    integer :: nx, ny, nz, ndof

    type(dim3) :: blocks, threads
    real(dp) :: norm_u, tol = 1d-8, tstart, tend


    nx = 512; ny = 512; nz = 512
    n_block = ny*nz/SZ
    ndof = nx*ny*nz
    n_iters = 100

    allocate (u_i(SZ, nx, n_block), u_o(SZ, nx, n_block))
    allocate (u_temp(SZ, nx, n_block))
    allocate (u_i_d(SZ, nx, n_block), u_o_d(SZ, nx, n_block))
    allocate (u_temp_d(SZ, nx, n_block))

    ! Cartesian order storage
    allocate (u_c_d(nx, ny, nz))
    allocate (u_c(nx, ny, nz))

    ! set a random field
    call random_number(u_i)

    ! move data to device
    u_i_d = u_i

    ! do a x to y reordering first and then a y to x
    blocks = dim3(nx/SZ, nz, ny/SZ)
    threads = dim3(SZ, SZ, 1)
    call reorder_x2y<<<blocks, threads>>>(u_temp_d, u_i_d, nz)

    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call reorder_y2x<<<blocks, threads>>>(u_o_d, u_temp_d, nz)

    ! move the result back to host
    u_o = u_o_d

    ! and check whether it matches the initial random field
    norm_u = norm2(u_o - u_i)
    if ( norm_u > tol ) then
      allpass = .false.
      write(stderr, '(a)') 'Check reorder x2y and y2x... failed'
    else
      write(stderr, '(a)') 'Check reorder x2y and y2x... passed'
    end if

    ! we reuse u_o_d so zeroize in any case
    u_o_d = 0

    ! u_temp_d is in y orientation, use y2z to reorder it into z direction
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call reorder_y2z<<<blocks, threads>>>(u_o_d, u_temp_d, nx, nz)

    ! store initial z oriented field
    u_temp = u_temp_d

    ! z oriented field into y
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call reorder_z2y<<<blocks, threads>>>(u_temp_d, u_o_d, nx, nz)

    u_o = u_temp_d

    ! compare two y oriented fields
    norm_u = norm2(u_o - u_temp)
    if ( norm_u > tol ) then
      allpass = .false.
      write(stderr, '(a)') 'Check reorder y2z and y2z... failed'
    else
      write(stderr, '(a)') 'Check reorder y2z and y2z... passed'
    end if

    ! reorder initial random field into z orientation
    blocks = dim3(nx, ny/SZ, 1)
    threads = dim3(SZ, 1, 1)
    call reorder_x2z<<<blocks, threads>>>(u_o_d, u_i_d, nz)

    ! z oriented field into x
    blocks = dim3(nx, ny/SZ, 1)
    threads = dim3(SZ, 1, 1)
    call reorder_z2x<<<blocks, threads>>>(u_temp_d, u_o_d, nz)
    u_o = u_temp_d

    ! compare two z oriented fields
    norm_u = norm2(u_o - u_i)
    if ( norm_u > tol ) then
      allpass = .false.
      write(stderr, '(a)') 'Check reorder x2z and z2x... failed'
    else
      write(stderr, '(a)') 'Check reorder x2z and z2x... passed'
    end if

    ! x ordering into Cartesian ordering
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call reorder_x2c<<<blocks, threads>>>(u_c_d, u_i_d, nz)

    ! sanitise u_o_d
    u_o_d = 0

    ! Cartesian ordering back to x ordering
    call reorder_c2x<<<blocks, threads>>>(u_o_d, u_c_d, nz)
    u_o = u_o_d

    ! now both u_o and u_i in x ordering, compare them
    norm_u = norm2(u_o - u_i)
    if ( norm_u > tol ) then
      allpass = .false.
      write(stderr, '(a)') 'Check reorder x2c and c2x... failed'
    else
      write(stderr, '(a)') 'Check reorder x2c and c2x... passed'
    end if

    if (allpass) then
      write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if

    ! Now the performance checks

    print*, 'Performance test: reorder_x2y'
    blocks = dim3(nx/SZ, nz, ny/SZ)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_x2y<<<blocks, threads>>>(u_o_d, u_i_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_x2z'
    blocks = dim3(nx, ny/SZ, 1)
    threads = dim3(SZ, 1, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_x2z<<<blocks, threads>>>(u_o_d, u_i_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_y2x'
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_y2x<<<blocks, threads>>>(u_o_d, u_i_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_y2z'
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_y2z<<<blocks, threads>>>(u_o_d, u_i_d, nx, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_z2x'
    blocks = dim3(nx, ny/SZ, 1)
    threads = dim3(SZ, 1, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_z2x<<<blocks, threads>>>(u_o_d, u_i_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_z2y'
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_z2y<<<blocks, threads>>>(u_o_d, u_i_d, nx, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_x2c'
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_x2c<<<blocks, threads>>>(u_c_d, u_i_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

    print*, 'Performance test: reorder_c2x'
    blocks = dim3(nx/SZ, ny/SZ, nz)
    threads = dim3(SZ, SZ, 1)
    call cpu_time(tstart)
    do i = 1, n_iters
      call reorder_c2x<<<blocks, threads>>>(u_o_d, u_c_d, nz)
    end do
    call cpu_time(tend)

    call checkperf(tend - tstart, n_iters, ndof, 2._dp)

  contains

    subroutine checkperf(t_tot, n_iters, ndof, consumed_bw)
      implicit none

      real(dp), intent(in) :: t_tot, consumed_bw
      integer, intent(in) :: n_iters, ndof

      real(dp) :: achievedBW, devBW
      integer :: ierr, memClockRt, memBusWidth

      ! BW utilisation and performance checks
      achievedBW = consumed_bw*n_iters*ndof*dp/t_tot

      print'(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'

      ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
      ierr = cudaDeviceGetAttribute(memBusWidth, &
        cudaDevAttrGlobalMemoryBusWidth, 0)
      devBW = 2*memBusWidth/8._dp*memClockRt*1000

      print'(a, f8.3, a)', 'Device BW:   ', devBW/2**30, ' GiB/s'
      print'(a, f5.2)', 'Effective BW util: %', achievedBW/devBW*100

    end subroutine checkperf

  end program test_cuda_reorder

