module m_cuda_allocator
   use m_allocator, only: allocator_t, field_t
   use m_common, only: dp

   implicit none

   type, extends(allocator_t) :: cuda_allocator_t
   contains
      procedure :: create_block => create_cuda_block
   end type cuda_allocator_t

   interface cuda_allocator_t
      module procedure cuda_allocator_init
   end interface cuda_allocator_t

   type, extends(field_t) :: cuda_field_t
      real(dp), device, pointer, private :: p_data_d(:)
      real(dp), device, pointer, contiguous :: data_d(:, :, :)
   end type cuda_field_t

   interface cuda_field_t
      module procedure cuda_field_init
   end interface cuda_field_t

contains

   function cuda_field_init(nx, ny, nz, sz, next, id) result(f)
      integer, intent(in) :: nx, ny, nz, sz, id
      type(cuda_field_t), pointer, intent(in) :: next
      type(cuda_field_t) :: f

      allocate (f%p_data_d(nx*ny*nz))
      ! will be removed, bounds remapping will be carried out by get_block.
      f%data_d(1:sz, 1:nx, 1:ny*nz/sz) => f%p_data_d
      f%refcount = 0
      f%next => next
      f%id = id
   end function cuda_field_init

   function cuda_allocator_init(nx, ny, nz, sz) result(allocator)
      integer, intent(in) :: nx, ny, nz, sz
      type(cuda_allocator_t) :: allocator

      allocator%allocator_t = allocator_t(nx, ny, nz, sz)
   end function cuda_allocator_init

   function create_cuda_block(self, next) result(ptr)
      class(cuda_allocator_t), intent(inout) :: self
      type(cuda_field_t), pointer, intent(in) :: next
      type(cuda_field_t), pointer :: newblock
      class(field_t), pointer :: ptr
      allocate (newblock)
      self%next_id = self%next_id + 1
      newblock = cuda_field_t(self%nx_padded, self%ny_padded, self%nz_padded, &
                              self%sz, next, id=self%next_id)
      ptr => newblock
   end function create_cuda_block

end module m_cuda_allocator
