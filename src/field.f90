module m_field

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C

  type :: field_t
     !! Memory block type holding both a data field and a pointer
     !! to the next block.  The `field_t` type also holds a integer
     !! `refcount` that counts the number of references to this
     !! field.  User code is currently responsible for incrementing
     !! the reference count.
    class(field_t), pointer :: next
    real(dp), pointer, private :: p_data(:)
    real(dp), pointer, contiguous :: data(:, :, :)
    integer :: dir
    integer :: data_loc
    integer :: refcount = 0
    integer :: id !! An integer identifying the memory block.
  contains
    procedure :: fill
    procedure :: get_shape
    procedure :: set_shape
    procedure :: set_data_loc
  end type field_t

  interface field_t
    module procedure field_init
  end interface field_t

  type :: flist_t
    !! Use for creating a list of field pointers
    class(field_t), pointer :: ptr
  end type flist_t

contains

  function field_init(ngrid, next, id) result(f)
    integer, intent(in) :: ngrid, id
    type(field_t), pointer, intent(in) :: next
    type(field_t) :: f

    allocate (f%p_data(ngrid))
    f%refcount = 0
    f%next => next
    f%id = id
  end function field_init

  subroutine fill(self, c)
    implicit none

    class(field_t) :: self
    real(dp), intent(in) :: c

    self%p_data(:) = c

  end subroutine fill

  subroutine set_data_loc(self, data_loc)
    class(field_t) :: self
    integer, intent(in) :: data_loc

    self%data_loc = data_loc

  end subroutine

  function get_shape(self) result(dims)
    implicit none

    class(field_t) :: self
    integer :: dims(3)

    dims = shape(self%data)

  end function get_shape

  subroutine set_shape(self, dims)
    implicit none

    class(field_t) :: self
    integer, intent(in) :: dims(3)

    self%data(1:dims(1), 1:dims(2), 1:dims(3)) => self%p_data

  end subroutine set_shape

end module m_field
