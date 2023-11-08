program xcompact
   use m_base_backend
   use m_cuda_backend
   use m_allocator
   use m_cuda_allocator
   use m_common, only: pi, globs_t
   use m_cuda_common
   use m_time_integrator, only: time_intg_t
   use m_tdsops, only: tdsops_t
   use m_cuda_tdsops, only: cuda_tdsops_t

   implicit none

   type(globs_t) :: globs
   class(base_backend_t), pointer :: backend
   class(allocator_t), pointer :: allocator
   type(time_intg_t) :: time_integrator
   type(dirps_t) :: xdirps, ydirps, zdirps

   type(cuda_backend_t), target :: cuda_backend

   type(cuda_allocator_t), target :: cuda_allocator

   cuda_allocator = cuda_allocator_t([SZ, 512, 512*512/SZ])
   allocator => cuda_allocator
   print*, 'allocator done'

   ! read L_x/y/z from the input file
   globs%Lx = 2*pi; globs%Ly = 2*pi; globs%Lz = 2*pi

   ! read ns from the input file
   globs%nx = 512; globs%ny = 512; globs%nz = 512

   ! set nprocs based on run time arguments
   globs%nproc_x = 1; globs%nproc_y = 1; globs%nproc_z = 1

   ! lets assume simple cases for now
   globs%nx_loc = globs%nx/globs%nproc_x
   globs%ny_loc = globs%ny/globs%nproc_y
   globs%nz_loc = globs%nz/globs%nproc_z

   globs%dx = globs%Lx/globs%nx
   globs%dy = globs%Ly/globs%ny
   globs%dz = globs%Lz/globs%nz

   xdirps%n = 512
   ydirps%n = 512
   zdirps%n = 512

   cuda_backend = cuda_backend_t(globs, allocator, xdirps, ydirps, zdirps)
   backend => cuda_backend
   print*, 'backend done'
   ! GPU only
   !allocate(cuda_allocator_t :: allocator)
   !allocate(cuda_backend_t :: backend)

   !cuda_backend = cuda_backend_t(allocator, xdirps, ydirps, zdirps)

   !allocator = cuda_allocator_t([SZ, 512, 512*512/SZ])
   !backend = cuda_backend_t(allocator, xdirps, ydirps, zdirps)
   time_integrator = time_intg_t(allocator=allocator, &
                                 backend=backend)

   print*, 'time integrator done'
   call time_integrator%run(1)
   print*, 'end'

end program xcompact
