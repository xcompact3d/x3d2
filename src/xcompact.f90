program xcompact
   use mpi

   use m_allocator
   use m_base_backend
   use m_common, only: pi, globs_t, set_pprev_pnext
   use m_solver, only: solver_t
   use m_time_integrator, only: time_intg_t
   use m_tdsops, only: tdsops_t

   use m_cuda_allocator
   use m_cuda_backend
   use m_cuda_common, only: SZ
   use m_cuda_tdsops, only: cuda_tdsops_t

   implicit none

   type(globs_t) :: globs
   class(base_backend_t), pointer :: backend
   class(allocator_t), pointer :: allocator
   type(solver_t) :: solver
   type(time_intg_t) :: time_integrator
   type(dirps_t) :: xdirps, ydirps, zdirps

   type(cuda_backend_t), target :: cuda_backend

   type(cuda_allocator_t), target :: cuda_allocator

   real(dp), allocatable, dimension(:, :, :) :: u, v, w

   real(dp) :: t_start, t_end
   integer :: nrank, nproc, ierr, ndevs, devnum

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'

   ierr = cudaGetDeviceCount(ndevs)
   ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
   ierr = cudaGetDevice(devnum)

   ! read L_x/y/z from the input file
   globs%Lx = 2*pi; globs%Ly = 2*pi; globs%Lz = 2*pi

   ! read ns from the input file
   globs%nx = 512; globs%ny = 512; globs%nz = 512

   ! set nprocs based on run time arguments
   globs%nproc_x = 1; globs%nproc_y = 1; globs%nproc_z = 1

   ! Lets allow a 1D decomposition for the moment
   !globs%nproc_x = nproc

   xdirps%nproc = globs%nproc_x
   ydirps%nproc = globs%nproc_y
   zdirps%nproc = globs%nproc_z

   ! Better if we move this somewhere else
   ! Set the pprev and pnext for each rank
   call set_pprev_pnext( &
      xdirps%pprev, xdirps%pnext, &
      ydirps%pprev, ydirps%pnext, &
      zdirps%pprev, zdirps%pnext, &
      xdirps%nproc, ydirps%nproc, zdirps%nproc, nrank &
   )

   print*, 'nrank:', nrank, 'xprev, xnext:', xdirps%pprev, xdirps%pnext
   print*, 'nrank:', nrank, 'yprev, ynext:', ydirps%pprev, ydirps%pnext
   print*, 'nrank:', nrank, 'zprev, znext:', zdirps%pprev, zdirps%pnext

   ! lets assume simple cases for now
   globs%nx_loc = globs%nx/globs%nproc_x
   globs%ny_loc = globs%ny/globs%nproc_y
   globs%nz_loc = globs%nz/globs%nproc_z

   globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
   globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
   globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

   globs%dx = globs%Lx/globs%nx
   globs%dy = globs%Ly/globs%ny
   globs%dz = globs%Lz/globs%nz

   xdirps%n = globs%nx_loc
   ydirps%n = globs%ny_loc
   zdirps%n = globs%nz_loc

   cuda_allocator = cuda_allocator_t([SZ, globs%nx_loc, globs%n_groups_x])
   allocator => cuda_allocator
   print*, 'allocator instantiated'

   cuda_backend = cuda_backend_t(globs, allocator, xdirps, ydirps, zdirps)
   backend => cuda_backend
   print*, 'backend instantiated'
   backend%nu = 1._dp

   allocate(u(SZ, globs%nx_loc, globs%n_groups_x))
   allocate(v(SZ, globs%nx_loc, globs%n_groups_x))
   allocate(w(SZ, globs%nx_loc, globs%n_groups_x))

   ! GPU only
   !allocate(cuda_allocator_t :: allocator)
   !allocate(cuda_backend_t :: backend)

   !cuda_backend = cuda_backend_t(allocator, xdirps, ydirps, zdirps)

   !allocator = cuda_allocator_t([SZ, 512, 512*512/SZ])
   !backend = cuda_backend_t(allocator, xdirps, ydirps, zdirps)

   time_integrator = time_intg_t(allocator=allocator, backend=backend)
   print*, 'time integrator instantiated'
   solver = solver_t(backend, time_integrator, xdirps, ydirps, zdirps)
   print*, 'solver instantiated'

   call cpu_time(t_start)

   call solver%run(100, u, v, w)

   call cpu_time(t_end)

   print*, 'Time: ', t_end - t_start

   print*, 'norms', norm2(u), norm2(v), norm2(w)

end program xcompact
