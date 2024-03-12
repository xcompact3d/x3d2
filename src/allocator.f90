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

      integer :: nx_padded, ny_padded, nz_padded, sz
      !> The id for the next allocated block.  This counter is
      !> incremented each time a new block is allocated.
      integer :: next_id = 0
      !> Padded dimensions for x, y, and z oriented fields
      integer :: xdims(3), ydims(3), zdims(3)
      !> Padded dimensions for natural Cartesian ordering
      integer :: cdims(3)
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

   interface allocator_t
      module procedure allocator_init
   end interface allocator_t

   type :: field_t
     !! Memory block type holding both a data field and a pointer
     !! to the next block.  The `field_t` type also holds a integer
     !! `refcount` that counts the number of references to this
     !! field.  User code is currently responsible for incrementing
     !! the reference count.
      class(field_t), pointer :: next
      real(dp), pointer, private :: p_data(:)
      real(dp), pointer, contiguous :: data(:, :, :)
      integer :: refcount = 0
      integer :: id !! An integer identifying the memory block.
   contains
      procedure :: set_shape
   end type field_t

   interface field_t
      module procedure field_init
   end interface field_t

   type :: flist_t
      class(field_t), pointer :: ptr
   end type flist_t

contains

   function field_init(nx, ny, nz, sz, next, id) result(f)
      integer, intent(in) :: nx, ny, nz, sz, id
      type(field_t), pointer, intent(in) :: next
      type(field_t) :: f

      allocate (f%p_data(nx*ny*nz))
      f%refcount = 0
      f%next => next
      f%id = id
   end function field_init

   subroutine set_shape(self, dims)
      implicit none

      class(field_t) :: self
      integer, intent(in) :: dims(3)

      self%data(1:dims(1), 1:dims(2), 1:dims(3)) => self%p_data

   end subroutine set_shape

   function allocator_init(nx, ny, nz, sz) result(allocator)
      integer, intent(in) :: nx, ny, nz, sz
      type(allocator_t) :: allocator

      integer :: nx_padded, ny_padded, nz_padded

      ! Apply padding based on sz
      nx_padded = nx
      ny_padded = ny
      nz_padded = nz

      allocator%nx_padded = nx_padded
      allocator%ny_padded = ny_padded
      allocator%nz_padded = nz_padded
      allocator%sz = sz

      allocator%xdims = [sz, nx_padded, ny_padded*nz_padded/sz]
      allocator%ydims = [sz, ny_padded, nx_padded*nz_padded/sz]
      allocator%zdims = [sz, nz_padded, nx_padded*ny_padded/sz]
      allocator%cdims = [nx_padded, ny_padded, nz_padded]
   end function allocator_init

   function create_block(self, next) result(ptr)
    !! Allocate memory for a new block and return a pointer to a new
    !! [[m_allocator(module):field_t(type)]] object.
      class(allocator_t), intent(inout) :: self
      type(field_t), pointer, intent(in) :: next
      type(field_t), pointer :: newblock
      class(field_t), pointer :: ptr
      self%next_id = self%next_id + 1
      allocate (newblock)
      newblock = field_t(self%nx_padded, self%ny_padded, self%nz_padded, &
                         self%sz, next, id=self%next_id)
      ptr => newblock
   end function create_block

   function get_block(self, dir) result(handle)
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
      integer, optional, intent(in) :: dir
      integer :: dims(3)
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

      ! set dims based on dir
      call handle%set_shape(self%xdims)
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
