program perf_cuda_reorder
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ
  use m_cuda_kernels_reorder, only: reorder_x2y, reorder_x2z, reorder_y2x, &
                                    reorder_y2z, reorder_z2x, reorder_z2y, &
                                    reorder_c2x, reorder_x2c
  use m_test_utils, only: write_perf_metric, write_perf_summary, &
                          write_device_bw_metric

  implicit none

  integer, parameter :: CASE_X2Y = 1
  integer, parameter :: CASE_X2Z = 2
  integer, parameter :: CASE_Y2X = 3
  integer, parameter :: CASE_Y2Z = 4
  integer, parameter :: CASE_Z2X = 5
  integer, parameter :: CASE_Z2Y = 6
  integer, parameter :: CASE_X2C = 7
  integer, parameter :: CASE_C2X = 8

  integer, parameter :: nx = 512
  integer, parameter :: ny = 512
  integer, parameter :: nz = 512
  integer, parameter :: n_iters = 100
  integer, parameter :: n_warmup = 10
  real(dp), parameter :: consumed_bw = 2.0_dp

  integer :: ierr
  integer :: n_block, ndof
  integer :: memClockRt, memBusWidth
  real(dp), allocatable :: u_i(:, :, :)
  real(dp), device, allocatable :: u_i_d(:, :, :), u_o_d(:, :, :), &
                                   u_temp_d(:, :, :)
  real(dp), device, allocatable :: u_c_d(:, :, :)
  type(dim3) :: blocks, threads

  n_block = ny*nz/SZ
  ndof = nx*ny*nz

  allocate (u_i(SZ, nx, n_block))
  allocate (u_i_d(SZ, nx, n_block), u_o_d(SZ, nx, n_block))
  allocate (u_temp_d(SZ, nx, n_block))
  allocate (u_c_d(nx, ny, nz))

  call random_number(u_i)
  u_i_d = u_i

  ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
  ierr = cudaDeviceGetAttribute(memBusWidth, &
                                cudaDevAttrGlobalMemoryBusWidth, 0)

  call run_case('cuda_reorder_x2y', CASE_X2Y)
  call run_case('cuda_reorder_x2z', CASE_X2Z)
  call run_case('cuda_reorder_y2x', CASE_Y2X)
  call run_case('cuda_reorder_y2z', CASE_Y2Z)
  call run_case('cuda_reorder_z2x', CASE_Z2X)
  call run_case('cuda_reorder_z2y', CASE_Z2Y)
  call run_case('cuda_reorder_x2c', CASE_X2C)
  call run_case('cuda_reorder_c2x', CASE_C2X)

  call write_device_bw_metric(memClockRt, memBusWidth)

contains

  subroutine run_case(label, case_id)
    character(len=*), intent(in) :: label
    integer, intent(in) :: case_id

    integer :: iter
    real(dp) :: tstart, tend

    print *, 'Performance test:', trim(label)

    call prepare_input(case_id)

    do iter = 1, n_warmup
      call launch_kernel(case_id)
    end do
    call sync_device()

    call cpu_time(tstart)
    do iter = 1, n_iters
      call launch_kernel(case_id)
    end do
    call sync_device()
    call cpu_time(tend)

    call write_perf_metric(label, tend - tstart, n_iters, ndof, consumed_bw)
    call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw, &
                            memClockRt, memBusWidth)
  end subroutine run_case

  subroutine prepare_input(case_id)
    integer, intent(in) :: case_id

    select case (case_id)
    case (CASE_Y2X, CASE_Y2Z)
      blocks = dim3(nx/SZ, nz, ny/SZ)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_x2y <  <  < blocks, threads >  >  > (u_temp_d, u_i_d, nz)
      call sync_device()
    case (CASE_Z2X, CASE_Z2Y)
      blocks = dim3(nx, ny/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_x2z <  <  < blocks, threads >  >  > (u_temp_d, u_i_d, nz)
      call sync_device()
    case (CASE_C2X)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_x2c <  <  < blocks, threads >  >  > (u_c_d, u_i_d, nz)
      call sync_device()
    case default
      continue
    end select
  end subroutine prepare_input

  subroutine launch_kernel(case_id)
    integer, intent(in) :: case_id

    select case (case_id)
    case (CASE_X2Y)
      blocks = dim3(nx/SZ, nz, ny/SZ)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_x2y <  <  < blocks, threads >  >  > (u_o_d, u_i_d, nz)
    case (CASE_X2Z)
      blocks = dim3(nx, ny/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_x2z <  <  < blocks, threads >  >  > (u_o_d, u_i_d, nz)
    case (CASE_Y2X)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_y2x <  <  < blocks, threads >  >  > (u_o_d, u_temp_d, nz)
    case (CASE_Y2Z)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
     call reorder_y2z <  <  < blocks, threads >  >  > (u_o_d, u_temp_d, nx, nz)
    case (CASE_Z2X)
      blocks = dim3(nx, ny/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_z2x <  <  < blocks, threads >  >  > (u_o_d, u_temp_d, nz)
    case (CASE_Z2Y)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
     call reorder_z2y <  <  < blocks, threads >  >  > (u_o_d, u_temp_d, nx, nz)
    case (CASE_X2C)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_x2c <  <  < blocks, threads >  >  > (u_c_d, u_i_d, nz)
    case (CASE_C2X)
      blocks = dim3(nx/SZ, ny/SZ, nz)
      threads = dim3(min(SZ, 32), min(SZ, 32), 1)
      call reorder_c2x <  <  < blocks, threads >  >  > (u_o_d, u_c_d, nz)
    end select
  end subroutine launch_kernel

  subroutine sync_device()
    ierr = cudaDeviceSynchronize()
  end subroutine sync_device

end program perf_cuda_reorder
