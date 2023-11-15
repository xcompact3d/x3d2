program test_cuda_tridiag
   use iso_fortran_env, only: stderr => error_unit
   use cudafor
   use mpi

   use m_common, only: dp, pi
   use m_cuda_common, only: SZ
   use m_cuda_exec_dist, only: exec_dist_transeq_3fused
   use m_cuda_sendrecv, only: sendrecv_fields, sendrecv_3fields
   use m_cuda_tdsops, only: cuda_tdsops_t

   implicit none

   logical :: allpass = .true.
   real(dp), allocatable, dimension(:, :, :) :: u, v, r_u
   real(dp), device, allocatable, dimension(:, :, :) :: &
      u_dev, v_dev, r_u_dev, & ! main fields u, v and result r_u
      du_dev, dud_dev, d2u_dev ! intermediate solution arrays
   real(dp), device, allocatable, dimension(:, :, :) :: &
      du_recv_s_dev, du_recv_e_dev, du_send_s_dev, du_send_e_dev, &
      dud_recv_s_dev, dud_recv_e_dev, dud_send_s_dev, dud_send_e_dev, &
      d2u_recv_s_dev, d2u_recv_e_dev, d2u_send_s_dev, d2u_send_e_dev

   real(dp), device, allocatable, dimension(:, :, :) :: &
      u_send_s_dev, u_send_e_dev, u_recv_s_dev, u_recv_e_dev, &
      v_send_s_dev, v_send_e_dev, v_recv_s_dev, v_recv_e_dev

   type(cuda_tdsops_t) :: der1st, der2nd

   integer :: n, n_block, i, j, k, n_halo, n_iters
   integer :: n_glob
   integer :: nrank, nproc, pprev, pnext, tag1=1234, tag2=1234
   integer :: srerr(12), mpireq(12)
   integer :: ierr, ndevs, devnum, memClockRt, memBusWidth

   type(dim3) :: blocks, threads
   real(dp) :: dx, dx_per, nu, norm_du, tol = 1d-8, tstart, tend
   real(dp) :: achievedBW, deviceBW, achievedBWmax, achievedBWmin

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'

   ierr = cudaGetDeviceCount(ndevs)
   ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
   ierr = cudaGetDevice(devnum)

   !print*, 'I am rank', nrank, 'I am running on device', devnum
   pnext = modulo(nrank - nproc + 1, nproc)
   pprev = modulo(nrank - 1, nproc)

   n_glob = 512*4
   n = n_glob/nproc
   n_block = 512*512/SZ
   n_iters = 100

   nu = 1._dp

   allocate(u(SZ, n, n_block), v(SZ, n, n_block), r_u(SZ, n, n_block))

   ! main input fields
   allocate(u_dev(SZ, n, n_block), v_dev(SZ, n, n_block))
   ! field for storing the result
   allocate(r_u_dev(SZ, n, n_block))
   ! intermediate solution fields
   allocate(du_dev(SZ, n, n_block))
   allocate(dud_dev(SZ, n, n_block))
   allocate(d2u_dev(SZ, n, n_block))

   dx_per = 2*pi/n_glob
   dx = 2*pi/(n_glob - 1)

   do k = 1, n_block
      do j = 1, n
         do i = 1, SZ
            u(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
            v(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
         end do
      end do
   end do

   ! move data to device
   u_dev = u
   v_dev = v

   n_halo = 4

   ! arrays for exchanging data between ranks
   allocate(u_send_s_dev(SZ, n_halo, n_block))
   allocate(u_send_e_dev(SZ, n_halo, n_block))
   allocate(u_recv_s_dev(SZ, n_halo, n_block))
   allocate(u_recv_e_dev(SZ, n_halo, n_block))
   allocate(v_send_s_dev(SZ, n_halo, n_block))
   allocate(v_send_e_dev(SZ, n_halo, n_block))
   allocate(v_recv_s_dev(SZ, n_halo, n_block))
   allocate(v_recv_e_dev(SZ, n_halo, n_block))

   allocate(du_send_s_dev(SZ, 1, n_block), du_send_e_dev(SZ, 1, n_block))
   allocate(du_recv_s_dev(SZ, 1, n_block), du_recv_e_dev(SZ, 1, n_block))
   allocate(dud_send_s_dev(SZ, 1, n_block), dud_send_e_dev(SZ, 1, n_block))
   allocate(dud_recv_s_dev(SZ, 1, n_block), dud_recv_e_dev(SZ, 1, n_block))
   allocate(d2u_send_s_dev(SZ, 1, n_block), d2u_send_e_dev(SZ, 1, n_block))
   allocate(d2u_recv_s_dev(SZ, 1, n_block), d2u_recv_e_dev(SZ, 1, n_block))

   ! preprocess the operator and coefficient arrays
   der1st = cuda_tdsops_t(n, dx_per, operation='first-deriv', &
                          scheme='compact6')
   der2nd = cuda_tdsops_t(n, dx_per, operation='second-deriv', &
                          scheme='compact6')

   blocks = dim3(n_block, 1, 1)
   threads = dim3(SZ, 1, 1)

   call cpu_time(tstart)
   do i = 1, n_iters
      u_send_s_dev(:, :, :) = u_dev(:, 1:4, :)
      u_send_e_dev(:, :, :) = u_dev(:, n - n_halo + 1:n, :)
      v_send_s_dev(:, :, :) = v_dev(:, 1:4, :)
      v_send_e_dev(:, :, :) = v_dev(:, n - n_halo + 1:n, :)


      ! halo exchange
      call sendrecv_fields(u_recv_s_dev, u_recv_e_dev, &
                           u_send_s_dev, u_send_e_dev, &
                           SZ*4*n_block, nproc, pprev, pnext)

      call sendrecv_fields(v_recv_s_dev, v_recv_e_dev, &
                           v_send_s_dev, v_send_e_dev, &
                           SZ*4*n_block, nproc, pprev, pnext)

      call exec_dist_transeq_3fused( &
         r_u_dev, &
         u_dev, u_recv_s_dev, u_recv_e_dev, &
         v_dev, v_recv_s_dev, v_recv_e_dev, &
         du_dev, dud_dev, d2u_dev, &
         du_send_s_dev, du_send_e_dev, du_recv_s_dev, du_recv_e_dev, &
         dud_send_s_dev, dud_send_e_dev, dud_recv_s_dev, dud_recv_e_dev, &
         d2u_send_s_dev, d2u_send_e_dev, d2u_recv_s_dev, d2u_recv_e_dev, &
         der1st, der2nd, nu, nproc, pprev, pnext, blocks, threads &
      )
   end do

   call cpu_time(tend)
   if (nrank == 0) print*, 'Total time', tend - tstart

   ! BW utilisation and performance checks
   ! 11 in the first phase, 5 in the second phase, 16 in total
   achievedBW = 16._dp*n_iters*n*n_block*SZ*dp/(tend - tstart)
   call MPI_Allreduce(achievedBW, achievedBWmax, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(achievedBW, achievedBWmin, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MIN, MPI_COMM_WORLD, ierr)

   if (nrank == 0) then
      print'(a, f8.3, a)', 'Achieved BW min: ', achievedBWmin/2**30, ' GiB/s'
      print'(a, f8.3, a)', 'Achieved BW max: ', achievedBWmax/2**30, ' GiB/s'
   end if

   ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
   ierr = cudaDeviceGetAttribute(memBusWidth, &
                                 cudaDevAttrGlobalMemoryBusWidth, 0)
   deviceBW = 2*memBusWidth/8._dp*memClockRt*1000

   if (nrank == 0) then
      print'(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
      print'(a, f5.2)', 'Effective BW util min: %', achievedBWmin/deviceBW*100
      print'(a, f5.2)', 'Effective BW util max: %', achievedBWmax/deviceBW*100
   end if

   ! check error
   r_u = r_u_dev
   r_u = r_u - (-v*v + 0.5_dp*u*u - nu*u)
   norm_du = norm2(r_u)
   norm_du = norm_du*norm_du/n_glob/n_block/SZ
   call MPI_Allreduce(MPI_IN_PLACE, norm_du, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)
   norm_du = sqrt(norm_du)

   if (nrank == 0) print*, 'error norm', norm_du

   if (nrank == 0) then
      if ( norm_du > tol ) then
         allpass = .false.
         write(stderr, '(a)') 'Check second derivatives... failed'
      else
         write(stderr, '(a)') 'Check second derivatives... passed'
      end if
   end if

   if (allpass) then
      if (nrank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
   else
      error stop 'SOME TESTS FAILED.'
   end if

   call MPI_Finalize(ierr)

end program test_cuda_tridiag

