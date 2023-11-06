program xcompact
   use m_base_backend
   use m_cuda_backend
   use m_allocator
   use m_cuda_allocator
   use m_common
   use m_cuda_common
   use m_time_integrator, only: time_intg_t

   implicit none

   !class(base_backend_t), allocatable :: backend
   !class(allocator_t), allocatable :: allocator
   type(time_intg_t) :: time_integrator
   type(dirps_t) :: xdirps, ydirps, zdirps

   type(cuda_backend_t), target :: cuda_backend
   class(base_backend_t), pointer :: backend

   type(cuda_allocator_t), target :: cuda_allocator
   class(allocator_t), pointer :: allocator

   cuda_allocator = cuda_allocator_t([SZ, 512, 512*512/SZ])
   allocator => cuda_allocator
   print*, 'allocator done'

   xdirps%n = 512
   ydirps%n = 512
   zdirps%n = 512

   cuda_backend = cuda_backend_t(allocator, xdirps, ydirps, zdirps)
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
   call time_integrator%run(1000)
   print*, 'end'

end program xcompact
