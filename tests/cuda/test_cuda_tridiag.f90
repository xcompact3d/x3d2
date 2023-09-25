program test_cuda_tridiag
   use iso_fortran_env, only: stderr => error_unit
   use cudafor
   use mpi

   use m_common, only: dp, SZ, pi
   use m_cuda_kernels_dist, only: der_univ_dist, der_univ_subs
   use m_derparams, only: der_1_vv, der_2_vv

   implicit none

   logical :: allpass = .true.
   real(dp), allocatable, dimension(:,:,:) :: u, du, u_b, u_e
   real(dp), device, allocatable, dimension(:,:,:) :: u_dev, du_dev
   real(dp), device, allocatable, dimension(:,:,:) :: u_recv_b_dev, u_recv_e_dev, &
                                                      u_send_b_dev, u_send_e_dev

   real(dp), device, allocatable, dimension(:,:,:) :: send_b_dev, send_e_dev, &
                                                      recv_b_dev, recv_e_dev

   real(dp), allocatable, dimension(:,:) :: coeffs_b, coeffs_e
   real(dp), allocatable, dimension(:) :: coeffs, dist_fr, dist_bc, &
                                          dist_sa, dist_sc

   real(dp), device, allocatable, dimension(:,:) :: coeffs_b_dev, coeffs_e_dev
   real(dp), device, allocatable, dimension(:) :: coeffs_dev, &
                                                  dist_fr_dev, dist_bc_dev, &
                                                  dist_sa_dev, dist_sc_dev

   integer :: n, n_block, i, j, k, n_halo, n_stencil, n_iters
   integer :: n_glob
   integer :: nrank, nproc, pprev, pnext, tag1=1234, tag2=1234
   integer, allocatable :: srerr(:), mpireq(:)
   integer :: ierr, ndevs, devnum, memClockRt, memBusWidth

   type(dim3) :: blocks, threads
   real(dp) :: dx, dx2, alfa, norm_du, tol = 1d-8, tstart, tend
   real(dp) :: achievedBW, deviceBW

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'

   ierr = cudaGetDeviceCount(ndevs)
   ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
   ierr = cudaGetDevice(devnum)
   print*, 'I am rank', nrank, 'I am running on device', devnum
   pnext = modulo(nrank-nproc+1, nproc)
   pprev = modulo(nrank-1, nproc)
   print*, 'rank', nrank, 'pnext', pnext, 'pprev', pprev
   allocate(srerr(nproc), mpireq(nproc))

   n_glob = 512*4
   n = n_glob/nproc
   n_block = 512*512/SZ
   n_iters = 1000

   allocate(u(SZ, n, n_block), du(SZ, n, n_block))
   allocate(u_dev(SZ, n, n_block), du_dev(SZ, n, n_block))

   dx = 2*pi/n_glob
   dx2 = dx*dx

   do k = 1, n_block
      do j = 1, n
         do i = 1, SZ
            u(i, j, k) = sin((j-1+nrank*n)*dx)
         end do
      end do
   end do

   ! move data to device
   u_dev = u

   ! set up the tridiagonal solver coeffs
   call der_2_vv(coeffs, coeffs_b, coeffs_e, dist_fr, dist_bc, &
                 dist_sa, dist_sc, n_halo, alfa, dx2, n, 'periodic')

   n_stencil = n_halo*2 + 1

   allocate(coeffs_b_dev(n_halo, n_stencil), coeffs_e_dev(n_halo, n_stencil))
   allocate(coeffs_dev(n_stencil))
   coeffs_b_dev(:, :) = coeffs_b(:, :); coeffs_e_dev(:, :) = coeffs_e(:, :)
   coeffs_dev(:) = coeffs(:)

   allocate(dist_fr_dev(n), dist_bc_dev(n), dist_sa_dev(n), dist_sc_dev(n))
   dist_fr_dev(:) = dist_fr(:); dist_bc_dev(:) = dist_bc(:)
   dist_sa_dev(:) = dist_sa(:); dist_sc_dev(:) = dist_sc(:)

   ! arrays for exchanging data between ranks
   allocate(u_send_b_dev(SZ, n_halo, n_block))
   allocate(u_send_e_dev(SZ, n_halo, n_block))
   allocate(u_recv_b_dev(SZ, n_halo, n_block))
   allocate(u_recv_e_dev(SZ, n_halo, n_block))

   allocate(send_b_dev(SZ, 1, n_block), send_e_dev(SZ, 1, n_block))
   allocate(recv_b_dev(SZ, 1, n_block), recv_e_dev(SZ, 1, n_block))

   blocks = dim3(n_block, 1, 1)
   threads = dim3(SZ, 1, 1)

   call cpu_time(tstart)
   do i = 1, n_iters
      u_send_b_dev(:,:,:) = u_dev(:,1:4,:)
      u_send_e_dev(:,:,:) = u_dev(:,n-n_halo+1:n,:)

      ! halo exchange
      if (nproc == 1) then
         u_recv_b_dev = u_send_e_dev
         u_recv_e_dev = u_send_b_dev
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(u_send_b_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(u_recv_e_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag1, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(u_send_e_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(u_recv_b_dev, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag2, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if


      call der_univ_dist<<<blocks, threads>>>( &
         du_dev, send_b_dev, send_e_dev, u_dev, u_recv_b_dev, u_recv_e_dev, &
         coeffs_b_dev, coeffs_e_dev, coeffs_dev, &
         n, alfa, dist_fr_dev, dist_bc_dev &
      )

      ! halo exchange for 2x2 systems
      if (nproc == 1) then
         recv_b_dev = send_e_dev
         recv_e_dev = send_b_dev
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(send_b_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(recv_e_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(send_e_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(recv_b_dev, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if

      call der_univ_subs<<<blocks, threads>>>( &
         du_dev, recv_b_dev, recv_e_dev, &
         n, alfa, dist_sa_dev, dist_sc_dev &
      )
   end do
   call cpu_time(tend)
   print*, 'Total time', tend-tstart

   achievedBW = 6._dp*n_iters*n*n_block*SZ*dp/(tend-tstart)
   print'(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'

   ierr = cudaDeviceGetAttribute(memClockRt, cudaDevAttrMemoryClockRate, 0)
   ierr = cudaDeviceGetAttribute(memBusWidth, cudaDevAttrGlobalMemoryBusWidth, 0)
   deviceBW = 2*memBusWidth/8._dp*memClockRt*1000

   print'(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
   print'(a, f5.2)', 'Effective BW utilization: %', achievedBW/deviceBW*100

   ! check error
   du = du_dev
   norm_du = norm2(u+du)/n_block
   print*, 'error norm', norm_du

   if ( norm_du > tol ) then
      allpass = .false.
      write(stderr, '(a)') 'Check second derivatives... failed'
   else
      write(stderr, '(a)') 'Check second derivatives... passed'
   end if

   if (allpass) then
      write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
   else
      error stop 'SOME TESTS FAILED.'
   end if

end program test_cuda_tridiag

