program test_dist_tridiag
   use iso_fortran_env, only: stderr => error_unit
   use mpi
   use omp_lib

   use m_common, only: dp, pi
   use m_omp_common, only: SZ
   use m_omp_kernels_dist, only: der_univ_dist_omp, der_univ_subs_omp
   use m_derparams, only: der_1_vv, der_2_vv

   implicit none

   logical :: allpass = .true.
   real(dp), allocatable, dimension(:,:,:) :: u, du, u_b, u_e
   real(dp), allocatable, dimension(:,:,:) :: u_recv_b, u_recv_e, &
                                              u_send_b, u_send_e

   real(dp), allocatable, dimension(:,:,:) :: send_b, send_e, &
                                              recv_b, recv_e

   real(dp), allocatable, dimension(:,:) :: coeffs_s, coeffs_e
   real(dp), allocatable, dimension(:) :: coeffs, dist_fr, dist_bc, dist_af, &
                                          dist_sa, dist_sc

   integer :: n, n_block, i, j, k, n_halo, n_stencil, n_iters, iters
   integer :: n_glob
   integer :: nrank, nproc, pprev, pnext, tag1=1234, tag2=1234
   integer, allocatable :: srerr(:), mpireq(:)
   integer :: ierr, ndevs, devnum, memClockRt, memBusWidth

   real(dp) :: dx, dx2, norm_du, tol = 1d-8, tstart, tend
   real(dp) :: achievedBW, deviceBW

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'

   print*, 'I am rank', nrank, 'I am running on device'
   pnext = modulo(nrank-nproc+1, nproc)
   pprev = modulo(nrank-1, nproc)
   print*, 'rank', nrank, 'pnext', pnext, 'pprev', pprev
   allocate(srerr(nproc), mpireq(nproc))

   n_glob = 512
   n = n_glob/nproc
   n_block = 512*512/SZ
   n_iters = 1000

   allocate(u(SZ, n, n_block), du(SZ, n, n_block))

   dx = 2*pi/n_glob
   dx2 = dx*dx

   do k = 1, n_block
      do j = 1, n
         do i = 1, SZ
            u(i, j, k) = sin((j-1+nrank*n)*dx)
         end do
      end do
   end do

   ! set up the tridiagonal solver coeffs
   call der_2_vv(coeffs, coeffs_s, coeffs_e, dist_fr, dist_bc, dist_af, &
                 dist_sa, dist_sc, n_halo, dx2, n, 'periodic')

   n_stencil = n_halo*2 + 1

   ! arrays for exchanging data between ranks
   allocate(u_send_b(SZ, n_halo, n_block))
   allocate(u_send_e(SZ, n_halo, n_block))
   allocate(u_recv_b(SZ, n_halo, n_block))
   allocate(u_recv_e(SZ, n_halo, n_block))

   allocate(send_b(SZ, 1, n_block), send_e(SZ, 1, n_block))
   allocate(recv_b(SZ, 1, n_block), recv_e(SZ, 1, n_block))

   !call cpu_time(tstart)
   tstart = omp_get_wtime()
   do iters = 1, n_iters
      ! first copy halo data into buffers
      !$omp parallel do
      do k = 1, n_block
         do j = 1, 4
            !$omp simd
            do i = 1, SZ
               u_send_b(i,j,k) = u(i,j,k)
               u_send_e(i,j,k) = u(i,n-n_halo+j,k)
            end do
            !$omp end simd
         end do
      end do
      !$omp end parallel do

      ! halo exchange
      if (nproc == 1) then
         u_recv_b = u_send_e
         u_recv_e = u_send_b
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(u_send_b, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(u_recv_e, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag1, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(u_send_e, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(u_recv_b, SZ*n_halo*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag2, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if


      !$omp parallel do
      do k = 1, n_block
         call der_univ_dist_omp( &
            du(:, :, k), send_b(:, :, k), send_e(:, :, k), u(:, :, k), &
            u_recv_b(:, :, k), u_recv_e(:, :, k), &
            coeffs_s, coeffs_e, coeffs, n, dist_fr, dist_bc, dist_af &
         )
      end do
      !$omp end parallel do

      ! halo exchange for 2x2 systems
      if (nproc == 1) then
         recv_b = send_e
         recv_e = send_b
      else
         ! MPI send/recv for multi-rank simulations
         call MPI_Isend(send_b, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(1), srerr(1))
         call MPI_Irecv(recv_e, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(2), srerr(2))
         call MPI_Isend(send_e, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pnext, tag2, MPI_COMM_WORLD, &
                        mpireq(3), srerr(3))
         call MPI_Irecv(recv_b, SZ*n_block, &
                        MPI_DOUBLE_PRECISION, pprev, tag1, MPI_COMM_WORLD, &
                        mpireq(4), srerr(4))

         call MPI_Waitall(4, mpireq, MPI_STATUSES_IGNORE, ierr)
      end if

      !$omp parallel do
      do k = 1, n_block
         call der_univ_subs_omp(du(:, :, k), recv_b(:, :, k), recv_e(:, :, k), &
                                n, dist_sa, dist_sc)
      end do
      !$omp end parallel do
   end do

   !call cpu_time(tend)
   tend = omp_get_wtime()
   print*, 'Total time', tend-tstart

   ! 3 in the first phase, 2 in the second phase, so 5 in total
   achievedBW = 5._dp*n_iters*n*n_block*SZ*dp/(tend-tstart)
   print'(a, f8.3, a)', 'Achieved BW: ', achievedBW/2**30, ' GiB/s'

   memClockRt = 3200
   memBusWidth = 64
   deviceBW = 2*memBusWidth/8._dp*memClockRt*1000000

   print'(a, f8.3, a)', 'Device BW:   ', deviceBW/2**30, ' GiB/s'
   print'(a, f5.2)', 'Effective BW utilization: %', achievedBW/deviceBW*100

   ! check error
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

   call MPI_Finalize(ierr)

end program test_dist_tridiag

