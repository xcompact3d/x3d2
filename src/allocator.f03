module m_allocator
  implicit none

  type :: allocator_t
     integer :: dims(3)
     integer :: id = 0
     class(memblock_t), pointer :: first => null()
   contains
     procedure :: get_block
     procedure :: release_block
     procedure :: create_block
     procedure :: get_block_ids
     procedure :: destroy
  end type allocator_t

  type :: memblock_t
     type(memblock_t), pointer :: next
     real, allocatable :: data(:, :, :)
     integer :: refcount = 0
     integer :: id
  end type memblock_t

  interface memblock_t
     module procedure memblock_constructor
  end interface memblock_t

contains

  function memblock_constructor(dims, next, id) result(m)
    integer, intent(in) :: dims(3), id
    type(memblock_t), pointer, intent(in) :: next
    type(memblock_t) :: m

    allocate(m%data(dims(1), dims(2), dims(3)))
    m%refcount = 0
    m%next => next
    m%id = id
  end function memblock_constructor

  function create_block(self, next) result(ptr)
    class(allocator_t), intent(inout) :: self
    type(memblock_t), pointer, intent(in) :: next
    type(memblock_t), pointer :: newblock
    class(memblock_t), pointer :: ptr
    self%id = self%id + 1
    allocate(newblock)
    newblock = memblock_t([8, 8, 8], next, id=self%id)
    ptr => newblock
  end function create_block

  function get_block(self) result(handle)
    !> Returns a pointer to the first available memory block, i.e. the
    !> head of the free memory block list.
    !> Example
    !>     f%data => get_memory_block()
    class(allocator_t), intent(inout) :: self
    class(memblock_t), pointer :: handle
    class(memblock_t), allocatable, target :: ptr
    ! If we're about to allocate the last block, extend the pool
    ! first.
    if(.not. associated(self%first)) then
       !> Construct a memblock_t. This effectively allocates
       !> storage space.
       self%first => self%create_block(next=self%first)
    end if
    handle => self%first
    self%first => self%first%next
    handle%next => null()
  end function get_block

  subroutine release_block(self, handle)
    !> Release memory block pointed to by handle to the pool.  It is
    !> pushed to the front of the free memory block list.
    class(allocator_t), intent(inout) :: self
    class(memblock_t), pointer :: handle
    handle%next => self%first
    self%first => handle
  end subroutine release_block

  subroutine destroy(self)
    !> Destroy each memory block in the free memory block list, from
    !> head to tail.  Deallocation of a memblock_t triggers all
    !> final (see def of memblock_t type) procedures and therefore
    !> deallocation of the block's storage space.
    class(allocator_t), intent(inout) :: self
    class(memblock_t), pointer :: current
    do
       if(.not. associated(self%first)) exit
       current => self%first
       self%first => self%first%next
       deallocate(current)
       self%id = self%id - 1
    end do
  end subroutine destroy

  function get_block_ids(self)
    class(allocator_t), intent(inout) :: self
    integer, allocatable :: get_block_ids(:)
    class(memblock_t), pointer :: current
    integer :: i

    current => self%first
    if(.not. associated(current)) then
       get_block_ids = [0]
    else
       i = current%id
       get_block_ids = [current%id]
       do
          if(.not. associated(current%next)) exit
          i = current%next%id
          get_block_ids = [get_block_ids, current%next%id]
          current => current%next
       end do
    end if
  end function get_block_ids

end module m_allocator
