!!! backends/omp/target/allocator.f90
!!
!! Implements an allocator specialised to OMP target offloading

module m_omptgt_allocator

  use iso_c_binding, only: c_ptr, c_f_pointer, c_sizeof
  use omp_lib, only: omp_target_alloc, omp_target_free, omp_get_default_device

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
    ! A device-resident field
    integer, private :: dev_id
    type(c_ptr), private :: dev_ptr
    real(dp), pointer, private :: p_data_tgt(:) => null()
    real(dp), pointer, contiguous :: data_tgt(:, :, :) => null()
  contains
    procedure :: destroy => omptgt_field_destroy
    procedure :: fill => fill_omptgt
    procedure :: get_shape => get_shape_omptgt
    procedure :: set_shape => set_shape_omptgt
  end type omptgt_field_t
  
  interface omptgt_field_t
    module procedure omptgt_field_init
  end interface omptgt_field_t
  
contains

  ! Constructor for the OMP target offload allocator
  type(omptgt_allocator_t) function omptgt_allocator_init(dims, sz) result(a)
    integer, intent(in) :: dims(3)
    integer, intent(in) :: sz

    a%allocator_t = allocator_t(dims, sz)
  end function omptgt_allocator_init

  ! Allocates a device-resident block
  function create_block_omptgt(self, next) result(ptr)
    class(omptgt_allocator_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: next
    type(omptgt_field_t), pointer :: newblock_tgt
    class(field_t), pointer :: ptr

    self%next_id = self%next_id + 1
    allocate(newblock_tgt)
    newblock_tgt = omptgt_field_t(self%ngrid, next, id=self%next_id)
    ptr => newblock_tgt

  end function create_block_omptgt

  ! Constructs a device-resident field
  type(omptgt_field_t) function omptgt_field_init(ngrid, next, id) result(f)
    integer, intent(in) :: ngrid
    class(field_t), pointer, intent(in) :: next
    integer, intent(in) :: id

    f%refcount = 0
    f%next => next
    f%id = id

    f%dev_id = omp_get_default_device()
    f%dev_ptr = omp_target_alloc(ngrid * c_sizeof(0.0_dp), f%dev_id)
    call c_f_pointer(f%dev_ptr, f%p_data_tgt, shape=[ngrid])
    
  end function omptgt_field_init

  subroutine omptgt_field_destroy(self)
    class(omptgt_field_t) :: self

    nullify(self%data_tgt)
    nullify(self%p_data_tgt)
    call omp_target_free(self%dev_ptr, self%dev_id)
  end subroutine

  ! Deallocates device-resident memory before deallocating the base type
  subroutine destroy(self)
    class(omptgt_allocator_t) :: self

    class(field_t), pointer :: ptr

    ptr => self%first
    do
      if (.not. associated(ptr)) then
        exit
      end if

      select type(ptr)
      type is(omptgt_field_t)
        call ptr%destroy()
      end select

      ptr => ptr%next
    end do

    call self%allocator_t%destroy()
  end subroutine

  subroutine fill_omptgt(self, c)
    class(omptgt_field_t) :: self
    real(dp), intent(in) :: c

    call fill_omptgt_(self%p_data_tgt, c, size(self%p_data_tgt))

  end subroutine fill_omptgt

  subroutine fill_omptgt_(p_data_tgt, c, n)
    real(dp), dimension(:), intent(inout) :: p_data_tgt
    real(dp), intent(in) :: c
    integer, intent(in) :: n

    integer :: i
    
    !$omp target teams distribute parallel do has_device_addr(p_data_tgt)
    do i = 1, n
      p_data_tgt(i) = c
    end do
    !$omp end target teams distribute parallel do

  end subroutine

  function get_shape_omptgt(self) result(dims)
    class(omptgt_field_t) :: self
    integer :: dims(3)

    dims = shape(self%data_tgt)
  end function

  subroutine set_shape_omptgt(self, dims)
    class(omptgt_field_t) :: self
    integer, intent(in) :: dims(3)

    call c_f_pointer(self%dev_ptr, self%data_tgt, shape=dims)

  end subroutine

end module m_omptgt_allocator
  
