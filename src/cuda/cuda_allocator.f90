module m_cudaallocator
   use m_allocator, only: allocator_t, field_t
   implicit none

   type, extends(allocator_t) :: cudaallocator_t
   contains
      procedure :: create_block => create_cuda_block
   end type cudaallocator_t

   type, extends(field_t) :: cudafield_t
      real, allocatable, device :: data_d(:, :, :)
   end type cudafield_t

   interface cudafield_t
      module procedure cuda_field_constructor
   end interface cudafield_t

contains

   function cuda_field_constructor(dims, next, id) result(m)
      integer, intent(in) :: dims(3), id
      type(cudafield_t), pointer, intent(in) :: next
      type(cudafield_t) :: m

      allocate (m%data_d(dims(1), dims(2), dims(3)))
      m%refcount = 0
      m%next => next
      m%id = id
   end function cuda_field_constructor

   function create_cuda_block(self, next) result(ptr)
      class(cudaallocator_t), intent(inout) :: self
      type(cudafield_t), pointer, intent(in) :: next
      type(cudafield_t), pointer :: newblock
      class(field_t), pointer :: ptr
      allocate (newblock)
      self%id = self%id + 1
      newblock = cudafield_t(self%dims, next, id=self%id)
      ptr => newblock
   end function create_cuda_block

end module m_cudaallocator
