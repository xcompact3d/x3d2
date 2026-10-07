!!! backends/omp/target/allocator.f90
!!
!! Implements an allocator specialised to OMP target offloading

module m_omptgt_allocator

  use iso_c_binding, only: c_ptr, c_f_pointer, &
                           c_sizeof, c_associated, c_null_ptr
  use omp_lib, only: omp_target_alloc, omp_target_free, &
                     omp_get_default_device, omp_get_num_devices, &
                     omp_get_initial_device

  use m_common, only: dp

  use m_allocator, only: allocator_t
  use m_mesh, only: mesh_t
  use m_field, only: field_t

  implicit none

  private
  public :: omptgt_allocator_t
  public :: omptgt_field_t

  type, extends(allocator_t) :: omptgt_allocator_t
  contains
    procedure :: create_block => create_block_omptgt
  end type omptgt_allocator_t

  interface omptgt_allocator_t
    module procedure omptgt_allocator_init
  end interface omptgt_allocator_t

  type, extends(field_t) :: omptgt_field_t
    !! A device-resident field.
    !!
    integer, private :: n
    integer, dimension(3), private :: dims
    integer, private :: dev_id
    type(c_ptr), private :: dev_ptr = c_null_ptr
  contains
    procedure :: fill => fill_omptgt
    procedure :: get_shape => get_shape_omptgt
    procedure :: set_shape => set_shape_omptgt
    procedure :: get_dev_ptr
    final :: omptgt_field_destroy
  end type omptgt_field_t

contains

  type(omptgt_allocator_t) function omptgt_allocator_init(dims, sz) result(a)
    !! Constructs the allocator, leaving the sizing of the memory blocks it
    !! hands out to the base allocator.
    integer, intent(in) :: dims(3)
    integer, intent(in) :: sz

    a%allocator_t = allocator_t(dims, sz)
  end function omptgt_allocator_init

  function create_block_omptgt(self, next) result(ptr)
    !! Creates a new device-resident block and links it in front of `next` in
    !! the allocator's list of blocks.
    class(omptgt_allocator_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: next
    class(field_t), pointer :: ptr

    self%next_id = self%next_id + 1
    allocate (omptgt_field_t :: ptr)
    select type (ptr)
    type is (omptgt_field_t)
      call omptgt_field_init(self%ngrid, next, id=self%next_id, f=ptr)
    class default
      error stop "Allocation failed to set type"
    end select

  end function create_block_omptgt

  subroutine omptgt_field_init(ngrid, next, id, f)
    !! Initialises a field and allocates its `ngrid` points of storage on the
    !! default target device. No device means offloading was not enabled at
    !! build time, so this errors out rather than silently running on the host.
    integer, intent(in) :: ngrid
    class(field_t), pointer, intent(in) :: next
    integer, intent(in) :: id
    type(omptgt_field_t) :: f

    f%refcount = 0
    f%next => next
    f%id = id

    if (omp_get_num_devices() < 1) then
      error stop "No OpenMP target device available, was offloading enabled?"
    end if
    f%dev_id = omp_get_default_device()
    if (f%dev_id == omp_get_initial_device()) then
      error stop "Device ID is HOST"
    end if

    f%n = ngrid
    f%dims = -1
    f%dev_ptr = omp_target_alloc(f%n*c_sizeof(0.0_dp), f%dev_id)
    if (.not. c_associated(f%dev_ptr)) then
      error stop "omp_target_alloc failed"
    end if

  end subroutine omptgt_field_init

  subroutine omptgt_field_destroy(self)
    !! Frees the device allocation when the field is finalised.
    type(omptgt_field_t) :: self

    if (c_associated(self%dev_ptr)) then
      call omp_target_free(self%dev_ptr, self%dev_id)
      self%dev_ptr = c_null_ptr
    end if
  end subroutine

  subroutine fill_omptgt(self, c)
    !! Sets every point of the field, padding included, to the constant `c`.
    class(omptgt_field_t) :: self
    real(dp), intent(in) :: c

    call fill_omptgt_(self%dev_ptr, c, self%n)

  end subroutine fill_omptgt

  subroutine fill_omptgt_(dev_ptr, c, n)
    !! Offloaded kernel behind `fill_omptgt`. Takes the device pointer rather
    !! than the field so the target region needs no mapping.
    type(c_ptr), intent(in) :: dev_ptr
    real(dp), intent(in) :: c
    integer, intent(in) :: n

    real(dp), dimension(:), pointer :: p_data_tgt
    integer :: i

    call c_f_pointer(dev_ptr, p_data_tgt, shape=[n])
    !$omp target is_device_ptr(dev_ptr)
    !$omp teams loop
    do i = 1, n
      p_data_tgt(i) = c
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  function get_shape_omptgt(self) result(dims)
    !! Returns the shape the field is currently viewed with, or -1s if
    !! `set_shape` has not been called yet.
    class(omptgt_field_t) :: self
    integer :: dims(3)

    dims = self%dims
  end function

  subroutine set_shape_omptgt(self, dims)
    !! Sets the shape the field is viewed with. Unlike the host and CUDA
    !! fields there is no array pointer to reshape here, so this only records
    !! the extents later used to cast the device pointer into an array.
    class(omptgt_field_t) :: self
    integer, intent(in) :: dims(3)

    if (product(dims) > self%n) then
      error stop "Trying to set shape of field greater than its capacity"
    end if
    self%dims = dims

  end subroutine

  type(c_ptr) function get_dev_ptr(self) result(ptr)
    !! Returns the device pointer to the field's storage, for use in the
    !! `is_device_ptr` clause of an offloaded region.
    class(omptgt_field_t) :: self

    ptr = self%dev_ptr
  end function

end module m_omptgt_allocator
