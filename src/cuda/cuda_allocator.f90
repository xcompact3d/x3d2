module m_cuda_allocator
   use m_allocator, only: allocator_t, field_t
   implicit none

   type, extends(allocator_t) :: cuda_allocator_t
   contains
      procedure :: create_block => create_cuda_block
   end type cuda_allocator_t

   type, extends(field_t) :: cuda_field_t
      real, allocatable, device :: data_d(:, :, :)
   end type cuda_field_t

   interface cuda_field_t
      module procedure cuda_field_constructor
   end interface cuda_field_t

contains

   function cuda_field_constructor(dims, next, id) result(m)
      integer, intent(in) :: dims(3), id
      type(cuda_field_t), pointer, intent(in) :: next
      type(cuda_field_t) :: m

      allocate (m%data_d(dims(1), dims(2), dims(3)))
      m%refcount = 0
      m%next => next
      m%id = id
   end function cuda_field_constructor

   function create_cuda_block(self, next) result(ptr)
      class(cuda_allocator_t), intent(inout) :: self
      type(cuda_field_t), pointer, intent(in) :: next
      type(cuda_field_t), pointer :: newblock
      class(field_t), pointer :: ptr
      allocate (newblock)
      self%id = self%id + 1
      newblock = cuda_field_t(self%dims, next, id=self%id)
      ptr => newblock
   end function create_cuda_block

end module m_cuda_allocator
