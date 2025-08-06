!!! backends/omp/target/allocator.f90
!!
!! Implements an allocator specialised to OMP target offloading

module m_omptgt_allocator

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
    real(dp), pointer, private :: p_data_tgt(:) => null()
    real(dp), pointer, contiguous :: data_tgt(:, :, :) => null()
  contains
    procedure :: fill => fill_omptgt
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

    print *, "Create OMPTGT block"

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

    ! allocate(f%p_data_tgt(ngrid))
    f%refcount = 0
    f%next => next
    f%id = id
    !$omp target enter data map(alloc:f%p_data_tgt(ngrid)) map(to:f%refcount) map(to:f%id) map(to:f%data_tgt)
    
  end function omptgt_field_init

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
        !$omp target exit data map(delete:ptr%p_data_tgt) map(delete:ptr%refcount) map(delete:ptr%id) map(delete:ptr%data_tgt)
      end select

      ptr => ptr%next
    end do

    call self%allocator_t%destroy()
  end subroutine

  subroutine fill_omptgt(self, c)
    class(omptgt_field_t) :: self
    real(dp), intent(in) :: c

    integer :: i
    
    !$omp target teams distribute parallel do
    do i = 1, size(self%p_data_tgt)
      self%p_data_tgt(i) = c
    end do
    !$omp end target teams distribute parallel do

  end subroutine fill_omptgt

end module m_omptgt_allocator
  
