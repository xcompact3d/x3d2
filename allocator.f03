module allocator
  implicit none

  type :: allocator_t
     integer :: SZ, n, m
     class(memblock_t), pointer :: first => null()
   contains
     procedure :: get_block
     procedure :: release_block
  end type allocator_t

  type :: memblock_t
     real, allocatable :: data(:, :, :)
     type(memblock_t), pointer :: next
     integer :: refcount
     integer :: id
  end type memblock_t

contains

  function memblock_constructor(dims, id, next) result(m)
    integer, intent(in) :: dims(3), id
    type(memory_block_t), pointer, intent(in) :: next
    type(memory_block_t) :: m

    allocate(m%data(dims(1), dims(2), dims(3)))
    m%refcount = 0
    m%next => next
    m%id = id
  end function memblock_constructor

  function get_block(self) result(handle)
    !> Returns a pointer to the first available memory block, i.e. the
    !> head of the free memory block list.
    !> Example
    !>     f%data => get_memory_block()
    class(allocator_t), intent(inout) :: self
    type(memory_block_t), pointer :: handle
    ! If we're about to allocate the last block, extend the pool
    ! first.
    if(.not. associated(self%first)) then
       allocate(newblock)
       !> Construct a memory_block_t. This effectively allocates
       !> storage space.
       newblock = memory_block_t([self%SZ, self%n, self%m], first, id=id)
       first => newblock
    end if
    handle => first
    first => first%next
    handle%next => null()
  end function get_block

  subroutine release_block(self, handle)
    !> Release memory block pointed to by handle to the pool.  It is
    !> pushed to the front of the free memory block list.
    class(allocator_t), intent(inout) :: self
    type(memblock_t), pointer :: handle
    handle%next => self%first
    self%first => handle
  end subroutine release_block

  subroutine destroy(self)
    !> Destroy each memory block in the free memory block list, from
    !> head to tail.  Deallocation of a memblock_t triggers all
    !> final (see def of memblock_t type) procedures and therefore
    !> deallocation of the block's storage space.
    class(allocator_t), intent(inout) :: self
    type(memblock_t), pointer :: current
    do
       if(.not. associated(self%first)) exit
       current => first
       first => first%next
       deallocate(current)
    end do
  end subroutine destroy

  subroutine deallocate_memblocks_data(self)
    type(memblock_t), intent(inout) :: self
    deallocate(self%data)
  end subroutine deallocate_memblocks_data

  function get_block_ids(self)
    class(allocator_t), intent(inout) :: self
    integer, allocatable :: get_block_ids(:)
    type(memblock_t), pointer :: current

    current => first
    get_block_ids = [current%id]
    do
       if(.not. associated(current%next)) exit
       get_block_ids = [get_block_ids, current%next%id]
       current => current%next
    end do
  end function get_block_ids

end module allocator
