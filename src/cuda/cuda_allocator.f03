module m_cudaallocator
  use m_allocator, only: allocator_t, memblock_t
  implicit none

  type, extends(allocator_t) :: cudaallocator_t
   contains
     procedure :: create_block => create_cuda_block
  end type cudaallocator_t

  type, extends(memblock_t) :: cudamemblock_t
     real, allocatable, device :: data_d(:, :, :)
  end type cudamemblock_t

  interface cudamemblock_t
     module procedure cuda_memblock_constructor
  end interface cudamemblock_t

contains

  function cuda_memblock_constructor(dims, next, id) result(m)
    integer, intent(in) :: dims(3), id
    type(cudamemblock_t), pointer, intent(in) :: next
    type(cudamemblock_t) :: m

    allocate(m%data_d(dims(1), dims(2), dims(3)))
    m%refcount = 0
    m%next => next
    m%id = id
  end function cuda_memblock_constructor

  function create_cuda_block(self, next) result(ptr)
    class(cudaallocator_t), intent(inout) :: self
    type(cudamemblock_t), pointer, intent(in) :: next
    type(cudamemblock_t), pointer :: newblock
    class(memblock_t), pointer :: ptr
    allocate(newblock)
    self%id = self%id + 1
    newblock = cudamemblock_t(self%dims, next, id=self%id + 1)
    ptr => newblock
  end function create_cuda_block

end module m_cudaallocator
