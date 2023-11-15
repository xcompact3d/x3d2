module m_allocator
   use m_common, only: dp

   implicit none

   type :: allocator_t
     !! An instance of type allocator_t is responsible for the
     !! maintenance of a linked list of instances of equal size
     !! [[m_allocator(module):field_t(type)]] objects:
     !!
     !! ```
     !!       ---- ---- ----     ---- ---- ----
     !! ...-->|id=1|data|next|-->|id=0|data|next|-->null()
     !!       ---- ---- ----     ---- ---- ----
     !! ```
     !!
     !! the last block's `next` pointer being non associated.
     !!
     !! User code can request access to a memory block by using the
     !! type bound procedure
     !! [[m_allocator(module):get_block(function)]].  If the list is
     !! not empty, a pointer to the first block on the list is
     !! returned and the block is detached from the list.  If the list
     !! is empty (i.e. all initially allocated blocks are currently
     !! referenced to) then a new block is allocated before a pointer
     !! to it is returned.
     !!
     !! In order to reuse memory it is important that user code
     !! release blocks when they are not needed anymore.  This is done
     !! by calling the type bound procedure
     !! [[m_allocator(module):release_block(subroutine)]].  The
     !! released block is then pushed in front of the block list.

      integer :: dims(3)
      !> The id for the next allocated block.  This counter is
      !> incremented each time a new block is allocated.
      integer :: next_id = 0
      !> The pointer to the first block on the list.  Non associated if
      !> the list is empty
      ! TODO: Rename first to head
      class(field_t), pointer :: first => null()
   contains
      procedure :: get_block
      procedure :: release_block
      procedure :: create_block
      procedure :: get_block_ids
      procedure :: destroy
   end type allocator_t

   type :: field_t
     !! Memory block type holding both a 3D data field and a pointer
     !! to the next block.  The `field_t` type also holds a integer
     !! `refcount` that counts the number of references to this
     !! field.  User code is currently responsible for incrementing
     !! the reference count.
      type(field_t), pointer :: next
      real(dp), allocatable :: data(:, :, :)
      integer :: refcount = 0
      integer :: id !! An integer identifying the memory block.
   end type field_t

   interface field_t
      module procedure field_constructor
   end interface field_t

contains

   function field_constructor(dims, next, id) result(m)
      integer, intent(in) :: dims(3), id
      type(field_t), pointer, intent(in) :: next
      type(field_t) :: m

      allocate (m%data(dims(1), dims(2), dims(3)))
      m%refcount = 0
      m%next => next
      m%id = id
   end function field_constructor

   function create_block(self, next) result(ptr)
    !! Allocate memory for a new block and return a pointer to a new
    !! [[m_allocator(module):field_t(type)]] object.
      class(allocator_t), intent(inout) :: self
      type(field_t), pointer, intent(in) :: next
      type(field_t), pointer :: newblock
      class(field_t), pointer :: ptr
      self%next_id = self%next_id + 1
      allocate (newblock)
      newblock = field_t(self%dims, next, id=self%next_id)
      ptr => newblock
   end function create_block

   function get_block(self) result(handle)
    !! Return a pointer to the first available memory block, i.e. the
    !! current head of the block list.  If the list is empty, allocate
    !! a new block with [[m_allocator(module):create_block(function)]]
    !! first.
    !!
    !! Example
    !! ```
    !! f%data => get_block()
    !! ```
      class(allocator_t), intent(inout) :: self
      class(field_t), pointer :: handle
      class(field_t), allocatable, target :: ptr
      ! If the list is empty, allocate a new block before returning a
      ! pointer to it.
      if (.not. associated(self%first)) then
         ! Construct a field_t. This effectively allocates
         ! storage space.
         self%first => self%create_block(next=self%first)
      end if
      handle => self%first
      self%first => self%first%next ! 2nd block becomes head block
      handle%next => null() ! Detach ex-head block from the block list
   end function get_block

   subroutine release_block(self, handle)
    !! Release memory block pointed to by HANDLE to the block list.
    !! It is pushed to the front of the block list, in other words it
    !! is made the head block.
      class(allocator_t), intent(inout) :: self
      class(field_t), pointer :: handle
      handle%next => self%first
      self%first => handle
   end subroutine release_block

   subroutine destroy(self)
    !! Go through the block list from head to tail, deallocating each
    !! memory block in turn.  Deallocation of a
    !! [[m_allocator(module):field_t(type)]] object automatically
    !! deallocates its internal allocatable
    !! [[field_t(type):data(variable)]] array.
      class(allocator_t), intent(inout) :: self
      type(field_t), pointer :: current
      do
         if (.not. associated(self%first)) exit
         current => self%first
         self%first => self%first%next
         deallocate (current)
         self%next_id = self%next_id - 1
      end do
   end subroutine destroy

   function get_block_ids(self)
    !! Utility function that returns a array made of the `id` of the
    !! block currently in the block list.  Return the array [0] if
    !! block list is empty.
      ! TODO: Block indices should start at 1 or return [-1] in case of
      ! empty block list.
      class(allocator_t), intent(inout) :: self
      integer, allocatable :: get_block_ids(:)
      class(field_t), pointer :: current
      integer :: i

      current => self%first
      if (.not. associated(current)) then
         get_block_ids = [0]
      else
         i = current%id
         get_block_ids = [current%id]
         do
            if (.not. associated(current%next)) exit
            i = current%next%id
            get_block_ids = [get_block_ids, current%next%id]
            current => current%next
         end do
      end if
   end function get_block_ids

end module m_allocator
