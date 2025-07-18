module m_allocator
  use iso_fortran_env, only: stderr => error_unit

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, NULL_LOC
  use m_field, only: field_t

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

    integer :: ngrid, sz
    !> The id for the next allocated block.  This counter is
    !> incremented each time a new block is allocated.
    integer :: next_id = 0
    !> padded dimensions and n_groups in all 'dir's
    integer, private :: dims_padded_dir(3, 4)
    integer, private :: n_groups_dir(3)
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
    procedure :: get_padded_dims
    procedure :: get_n_groups
  end type allocator_t

  interface allocator_t
    module procedure allocator_init
  end interface allocator_t

contains

  function allocator_init(dims, sz) result(allocator)
    integer, intent(in) :: dims(3), sz
    type(allocator_t) :: allocator

    integer :: nx, ny, nz, nx_padded, ny_padded, nz_padded

    nx = dims(1); ny = dims(2); nz = dims(3)

    ! Apply padding based on sz
    nx_padded = nx - 1 + mod(-(nx - 1), sz) + sz
    ny_padded = ny - 1 + mod(-(ny - 1), sz) + sz
    ! Current reorder functions do not require a padding in z-direction.
    nz_padded = nz

    allocator%ngrid = nx_padded*ny_padded*nz_padded
    allocator%sz = sz

    allocator%n_groups_dir(1:3) = [ny_padded*nz_padded/sz, &
                                   nx_padded*nz_padded/sz, &
                                   nx_padded*ny_padded/sz]

    allocator%dims_padded_dir(:, 1) = [sz, nx_padded, &
                                       allocator%n_groups_dir(1)]
    allocator%dims_padded_dir(:, 2) = [sz, ny_padded, &
                                       allocator%n_groups_dir(2)]
    allocator%dims_padded_dir(:, 3) = [sz, nz_padded, &
                                       allocator%n_groups_dir(3)]
    allocator%dims_padded_dir(:, 4) = [nx_padded, ny_padded, nz_padded]

  end function allocator_init

  function create_block(self, next) result(ptr)
    !! Allocate memory for a new block and return a pointer to a new
    !! [[m_allocator(module):field_t(type)]] object.
    class(allocator_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: next
    type(field_t), pointer :: newblock
    class(field_t), pointer :: ptr
    self%next_id = self%next_id + 1
    allocate (newblock)
    select type (next)
    type is (field_t)
      newblock = field_t(self%ngrid, next, id=self%next_id)
    class default
      error stop "Incorrect overloading for create_block"
    end select
    ptr => newblock
  end function create_block

  function get_block(self, dir, data_loc) result(handle)
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
    integer, intent(in) :: dir
    integer, intent(in), optional :: data_loc
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

    ! Store direction info in the field type.
    handle%dir = dir
    if (present(data_loc)) then
      handle%data_loc = data_loc
    else
      handle%data_loc = NULL_LOC
    end if

    ! Set dims based on direction
    dims = self%dims_padded_dir(1:3, dir)

    ! Apply bounds remapping based on requested direction
    call handle%set_shape(dims)
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

  function get_padded_dims(self, dir) result(dims)
    implicit none

    class(allocator_t), intent(inout) :: self
    integer, intent(in) :: dir
    integer :: dims(3)

    dims = self%dims_padded_dir(1:3, dir)
  end function get_padded_dims

  function get_n_groups(self, dir) result(n_groups)
    implicit none

    class(allocator_t), intent(inout) :: self
    integer, intent(in) :: dir
    integer :: n_groups

    n_groups = self%n_groups_dir(dir)
  end function get_n_groups

end module m_allocator
