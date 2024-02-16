program test_omp_transeq
   use iso_fortran_env, only: stderr => error_unit
   use mpi

   use m_allocator, only: allocator_t, field_t
   use m_common, only: dp, pi, globs_t, set_pprev_pnext
   use m_omp_common, only: SZ
   use m_omp_sendrecv, only: sendrecv_fields
   use m_omp_backend, only: omp_backend_t, transeq_x_omp, base_backend_t
   use m_tdsops, only: dirps_t, tdsops_t
   use m_solver, only: allocate_tdsops

   implicit none

   logical :: allpass = .true.
   class(field_t), pointer :: u, v, w
   class(field_t), pointer :: du, dv, dw
   real(dp), dimension(:, :, :), allocatable :: r_u

   integer :: n, n_block, i, j, k
   integer :: n_glob
   integer :: nrank, nproc
   integer :: ierr

   real(dp) :: dx, dx_per, nu, norm_du, tol = 1d-8, tstart, tend

   type(globs_t) :: globs
   class(base_backend_t), pointer :: backend
   class(allocator_t), pointer :: allocator

   type(omp_backend_t), target :: omp_backend
   type(allocator_t), target :: omp_allocator
   type(dirps_t) :: xdirps, ydirps, zdirps

   ! Initialise variables and arrays
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   globs%nx = 64
   globs%ny = 64
   globs%nz = 64

   globs%nx_loc = globs%nx/nproc
   globs%ny_loc = globs%ny/nproc
   globs%nz_loc = globs%nz/nproc

   globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
   globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
   globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

   xdirps%nproc = nproc
   ydirps%nproc = nproc
   zdirps%nproc = nproc

   call set_pprev_pnext( &
      xdirps%pprev, xdirps%pnext, &
      ydirps%pprev, ydirps%pnext, &
      zdirps%pprev, zdirps%pnext, &
      xdirps%nproc, ydirps%nproc, zdirps%nproc, nrank &
   ) 
   
   xdirps%n = globs%nx_loc
   ydirps%n = globs%ny_loc
   zdirps%n = globs%nz_loc

   xdirps%n_blocks = globs%n_groups_x
   ydirps%n_blocks = globs%n_groups_y
   zdirps%n_blocks = globs%n_groups_z

   omp_allocator = allocator_t([SZ, globs%nx_loc, globs%n_groups_x])
   allocator => omp_allocator
   print*, 'OpenMP allocator instantiated'

   omp_backend = omp_backend_t(globs, allocator)
   backend => omp_backend
   print*, 'OpenMP backend instantiated'


   if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'

   n_glob = globs%nx_loc
   n = n_glob/nproc
   n_block = xdirps%n_blocks 

   nu = 1._dp
   omp_backend%nu = nu

   
   u => allocator%get_block()
   v => allocator%get_block()
   w => allocator%get_block()

   du => allocator%get_block()
   dv => allocator%get_block()
   dw => allocator%get_block()

   dx_per = 2*pi/n_glob
   dx = 2*pi/(n_glob - 1)
   globs%dx = dx


   do k = 1, n_block
      do j = 1, n
         do i = 1, SZ
            u%data(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
            v%data(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
         end do
      end do
   end do
   w%data(:, :, :) = 0.d0

   call allocate_tdsops(xdirps, globs%nx_loc, dx_per, omp_backend)

   call cpu_time(tstart)
   call transeq_x_omp(omp_backend, du, dv, dw, u, v, w, xdirps)
   call cpu_time(tend)

   if (nrank == 0) print*, 'Total time', tend - tstart

   allocate(r_u(SZ, n, n_block))

   ! check error
   r_u = dv%data - (-v%data*v%data + 0.5_dp*u%data*u%data - nu*u%data)
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

end program test_omp_transeq

