module m_slab_cuda
   use m_allocator, only: cudaallocator_t, cudamemblock_t
   use m_slab, only: slab
   use m_diffengine, only: diffengine

   type, extends(slab) :: slab_cuda_t
      integer :: dims(3)
      type(cudaallocator_t), pointer :: allocator
      type(cudamemblock_t), pointer :: u1, u2, u3
   contains
      procedure, public :: transport
      procedure, public :: div
      procedure, public :: u => get_component_ptr
   end type slab_cuda_t

   interface slab_cuda_t
      module procedure make_slab_cuda
   end interface slab_cuda_t

contains

   function make_slab_cuda(allocator, diffeng, diffeng2) result(self)
      type(cudaallocator_t), pointer, intent(in) :: allocator
      type(diffengine), intent(in) :: diffeng, diffeng2
      type(slab_cuda_t) :: self

      call mpi_comm_rank(MPI_COMM_WORLD, self%rankid, errcode)
      call mpi_comm_size(MPI_COMM_WORLD, self%nranks, errcode)

      self%allocator => allocator
      self%diffeng = diffeng
      self%diffeng2 = diffeng2

      blockptr1 => self%allocator%get_block()
      blockptr2 => self%allocator%get_block()
      blockptr3 => self%allocator%get_block()
   end function make_slab_cuda

   function get_component_ptr(self, i) result(ptr)
      class(slab_cuda_t), intent(in) :: self
      integer, intent(in) :: i
      real, pointer :: ptr(:, :, :)

      select case (i)
      case (1)
         ptr => u1%data
      case (2)
         ptr => u2%data
      case (3)
         ptr => u3%data
      end select
   end function get_component_ptr

   pure function transport(self)
      class(slab_cuda_t), intent(in) :: self
      class(slab_cuda_t), allocatable :: transport

      allocate (slab_cuda_t :: transport)
      transport = cuda_simd(&
           & self%allocator, self%diffeng, self%diffeng2 &
           & )
      tranport%set_storage()
      du => self%allocator%get_block()
      d2u => self%allocator%get_block()
      dud => self%allocator%get_block()
      select type (tranport)
      type is (slab_cuda_t)
         call transport_dir <  <  < grid, tBlock >  >  > ( &
              & self%u(1), self%u(1), transport%u(1))
         call transport_dir <  <  < grid, tBlock >  >  > ( &
              & self%u(2), self%u(2), transport%u(2))
         call transport_dir <  <  < grid, tBlock >  >  > ( &
              & self%u(3), self%u(3), transport%u(3))
      class default
         error stop
      end select
   end function transport

   attributes(global) subroutine transport_dir(u, u_dir, rslt)
      real, device, intent(in) :: u(:, :, :)
      real, device, intent(in) :: u_dir(:, :, :)
      real, device, intent(out) :: rslt

      integer, parameter :: n = dim(u, 2)
      real, dimension(n) :: du, d2u, dud

      layers: do k = 1, size(u, 3)
         call diffeng%diff(u(:, :, k), du%data)
         call diffeng2%diff(u(:, :, k), d2u%data)
         do j = 1, size(u, 2)
            !$omp simd
            do i = 1, size(u, 1)
               u2%data(i, j) = u(i, j, k) * u_dir(i, j, k)
            end do
            !$omp end simd
         end do
         call diffeng%diff(du2%data, u2%data)
         reshape(transport_dir, shape(u))
         do j = 1, size(u, 2)
            !$omp simd
            do i = 1, size(u, 1)
               rslt(i, j, k) = -0.5 * &
                    (u(i, j, k) * du%data(i, j) + du2%data(i, j)) &
                    & + xnu * d2u%data(i, j)
               !$omp end simd
            end do
         end do
      end do layers

      call self%allocator%release_block(du)
      call self%allocator%release_block(d2u)
      call self%allocator%release_block(u2)
      call self%allocator%release_block(du2)
      end function transport_dir

      subroutine div(self, rslt)
         class(slab_cuda_t), intent(in) :: self
         class(slab_cuda_t), intent(inout) :: rslt

         select type (rslt)
         type is (slab_cuda_t)
            rslt%u = self%u + 1.
            rslt%v = self%v + 1.
            rslt%w = self%w + 1.
         class default
            error stop
         end select
      end subroutine div
   end module m_slab_cuda
