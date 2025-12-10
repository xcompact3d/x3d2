module m_field
  !! Field data structure module for managing computational grid data.
  !!
  !! This module provides the field_t type for storing 3D scalar fields
  !! on the computational grid. Fields can be organised in linked lists
  !! for memory management and support different data orientations
  !! (x-pencil, y-pencil, z-pencil).

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C

  type :: field_t
    !! Memory block type holding a 3D scalar field with metadata.
    !!
    !! The field_t type stores both a data field and a pointer to the next
    !! block, enabling linked list structures for memory management. The type
    !! tracks a reference count (currently managed by user code), data
    !! orientation (x-, y-, or z-pencil), and data location on the staggered grid.
    class(field_t), pointer :: next             !! Pointer to next field in linked list
    real(dp), pointer, private :: p_data(:)     !! 1D array storage for data
    real(dp), pointer, contiguous :: data(:, :, :)  !! 3D view of data array
    integer :: dir                              !! Data direction (DIR_X, DIR_Y, DIR_Z, or DIR_C)
    integer :: data_loc                         !! Data location flag (VERT, CELL, etc.)
    integer :: refcount = 0                     !! Reference count for memory management
    integer :: id                               !! Unique identifier for this memory block
  contains
    procedure :: fill          !! Fill field with a constant value
    procedure :: get_shape     !! Get 3D dimensions of data array
    procedure :: set_shape     !! Set 3D dimensions by reshaping p_data
    procedure :: set_data_loc  !! Set data location flag
  end type field_t

  interface field_t
    module procedure field_init
  end interface field_t

  type :: flist_t
    !! Wrapper type for creating arrays of field pointers.
    !!
    !! This type is used to create lists or arrays of field pointers,
    !! useful for managing multiple fields such as velocity components
    !! or transported scalar species.
    class(field_t), pointer :: ptr  !! Pointer to a field
  end type flist_t

contains

  function field_init(ngrid, next, id) result(f)
    !! Initialise a new field with allocated memory.
    !!
    !! Creates a new field_t instance with allocated storage for ngrid points.
    !! The field is linked to the next field in the list and assigned a unique ID.
    integer, intent(in) :: ngrid  !! Total number of grid points to allocate
    type(field_t), pointer, intent(in) :: next  !! Pointer to next field in linked list
    integer, intent(in) :: id     !! Unique identifier for this field
    type(field_t) :: f            !! Initialised field

    allocate (f%p_data(ngrid))
    f%refcount = 0
    f%next => next
    f%id = id
  end function field_init

  subroutine fill(self, c)
    !! Fill the entire field with a constant value.
    !!
    !! Sets all grid points in the field to the specified constant value.
    implicit none

    class(field_t) :: self        !! Field to fill
    real(dp), intent(in) :: c     !! Constant value to fill with

    self%p_data(:) = c

  end subroutine fill

  subroutine set_data_loc(self, data_loc)
    !! Set the data location flag for this field.
    !!
    !! The data location specifies where on the staggered grid the data
    !! is located (e.g., VERT, CELL, X_FACE, etc.).
    class(field_t) :: self           !! Field to modify
    integer, intent(in) :: data_loc  !! Data location flag

    self%data_loc = data_loc

  end subroutine

  function get_shape(self) result(dims)
    !! Get the 3D dimensions of the field data.
    !!
    !! Returns the current shape of the 3D data array.
    implicit none

    class(field_t) :: self  !! Field to query
    integer :: dims(3)      !! Array dimensions [nx, ny, nz]

    dims = shape(self%data)

  end function get_shape

  subroutine set_shape(self, dims)
    !! Reshape the field data to specified 3D dimensions.
    !!
    !! Maps the 1D storage array (p_data) to a 3D view with the specified
    !! dimensions. The total size must match the allocated storage.
    implicit none

    class(field_t) :: self        !! Field to reshape
    integer, intent(in) :: dims(3)  !! Target dimensions [nx, ny, nz]

    self%data(1:dims(1), 1:dims(2), 1:dims(3)) => self%p_data

  end subroutine set_shape

end module m_field
