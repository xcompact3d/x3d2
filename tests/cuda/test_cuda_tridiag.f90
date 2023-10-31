program test_cuda_tridiag
   use iso_fortran_env, only: stderr => error_unit
   use cudafor
   use mpi

   use m_common, only: dp, pi
   use m_cuda_common, only: SZ
   use m_cuda_kernels_dist, only: der_univ_dist, der_univ_subs
   use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init

   implicit none

   logical :: allpass = .true.
   real(dp), allocatable, dimension(:, :, :) :: u, du
   real(dp), device, allocatable, dimension(:, :, :) :: u_dev, du_dev
   real(dp), device, allocatable, dimension(:, :, :) :: &
      u_recv_s_dev, u_recv_e_dev, u_send_s_dev, u_send_e_dev

   real(dp), device, allocatable, dimension(:, :, :) :: &
      send_s_dev, send_e_dev, recv_s_dev, recv_e_dev

   type(cuda_tdsops_t) :: tdsops

   integer :: n, n_block, i, j, k, n_halo, n_iters
   integer :: n_glob
   integer :: nrank, nproc, pprev, pnext, tag1=1234, tag2=1234
   integer :: srerr(4), mpireq(4)
   integer :: ierr, ndevs, devnum, memClockRt, memBusWidth

   type(dim3) :: blocks, threads
   real(dp) :: dx, dx_per, norm_du, tol = 1d-8, tstart, tend
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

   allocate(u(SZ, n, n_block), du(SZ, n, n_block))
   allocate(u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))

   dx_per = 2*pi/n_glob
   dx = 2*pi/(n_glob - 1)

   do k = 1, n_block
      do j = 1, n
         do i = 1, SZ
            u(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
         end do
      end do
   end do

   ! move data to device
   u_dev = u

   n_halo = 4

   ! arrays for exchanging data between ranks
   allocate(u_send_s_dev(SZ, n_halo, n_block))
   allocate(u_send_e_dev(SZ, n_halo, n_block))
   allocate(u_recv_s_dev(SZ, n_halo, n_block))
   allocate(u_recv_e_dev(SZ, n_halo, n_block))

   allocate(send_s_dev(SZ, 1, n_block), send_e_dev(SZ, 1, n_block))
   allocate(recv_s_dev(SZ, 1, n_block), recv_e_dev(SZ, 1, n_block))

   ! preprocess the operator and coefficient arrays
   tdsops = cuda_tdsops_init(n, dx_per, operation='second-deriv', &
                             scheme='compact6')

   blocks = dim3(n_block, 1, 1)
   threads = dim3(SZ, 1, 1)

   call cpu_time(tstart)
   do i = 1, n_iters
      u_send_s_dev(:, :, :) = u_dev(:, 1:4, :)
      u_send_e_dev(:, :, :) = u_dev(:, n - n_halo + 1:n, :)

      ! halo exchange
      if (nproc == 1) then
         u_recv_s_dev = u_send_e_dev
         u_recv_e_dev = u_send_s_dev
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(u_send_s_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(u_recv_e_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag1, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(u_send_e_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(u_recv_s_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag2, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if


      call der_univ_dist<<<blocks, threads>>>( &
         du_dev, send_s_dev, send_e_dev, u_dev, u_recv_s_dev, u_recv_e_dev, &
         tdsops%coeffs_s_dev, tdsops%coeffs_e_dev, tdsops%coeffs_dev, &
         n, tdsops%dist_fw_dev, tdsops%dist_bw_dev, tdsops%dist_af_dev &
      )

      ! halo exchange for 2x2 systems
      if (nproc == 1) then
         recv_s_dev = send_e_dev
         recv_e_dev = send_s_dev
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(send_s_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(recv_e_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(send_e_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(recv_s_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if

      call der_univ_subs<<<blocks, threads>>>( &
         du_dev, recv_s_dev, recv_e_dev, &
         n, tdsops%dist_sa_dev, tdsops%dist_sc_dev &
      )
   end do

   call cpu_time(tend)
   if (nrank == 0) print*, 'Total time', tend - tstart

   ! BW utilisation and performance checks
   ! 4 in the first phase, 2 in the second phase, 6 in total
   achievedBW = 6._dp*n_iters*n*n_block*SZ*dp/(tend - tstart)
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
   du = du_dev
   norm_du = norm2(u + du)
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

