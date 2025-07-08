module m_cuda_allocator
  use m_allocator, only: allocator_t
  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type, extends(allocator_t) :: cuda_allocator_t
  contains
    procedure :: create_block => create_cuda_block
  end type cuda_allocator_t

  interface cuda_allocator_t
    module procedure cuda_allocator_init
  end interface cuda_allocator_t

  type, extends(field_t) :: cuda_field_t
    real(dp), device, pointer, private :: p_data_d(:)
    real(dp), device, pointer, contiguous :: data_d(:, :, :)
  contains
    procedure :: fill => fill_cuda
    procedure :: get_shape => get_shape_cuda
    procedure :: set_shape => set_shape_cuda
  end type cuda_field_t

  interface cuda_field_t
    module procedure cuda_field_init
  end interface cuda_field_t

contains

  function cuda_field_init(ngrid, next, id) result(f)
    integer, intent(in) :: ngrid, id
    type(cuda_field_t), pointer, intent(in) :: next
    type(cuda_field_t) :: f

    allocate (f%p_data_d(ngrid))
    f%refcount = 0
    f%next => next
    f%id = id
  end function cuda_field_init

  subroutine fill_cuda(self, c)
    implicit none

    class(cuda_field_t) :: self
    real(dp), intent(in) :: c

    self%p_data_d = c

  end subroutine fill_cuda

  function get_shape_cuda(self) result(dims)
    implicit none

    class(cuda_field_t) :: self
    integer :: dims(3)

    dims = shape(self%data_d)

  end function get_shape_cuda

  subroutine set_shape_cuda(self, dims)
    implicit none

    class(cuda_field_t) :: self
    integer, intent(in) :: dims(3)

    self%data_d(1:dims(1), 1:dims(2), 1:dims(3)) => self%p_data_d

  end subroutine set_shape_cuda

  function cuda_allocator_init(nx, ny, nz, sz) result(allocator)
    integer, intent(in) :: nx, ny, nz, sz
    type(cuda_allocator_t) :: allocator

    allocator%allocator_t = allocator_t(nx, ny, nz, sz)
  end function cuda_allocator_init

  function create_cuda_block(self, next) result(ptr)
    class(cuda_allocator_t), intent(inout) :: self
    type(cuda_field_t), pointer, intent(in) :: next
    type(cuda_field_t), pointer :: newblock
    class(field_t), pointer :: ptr
    allocate (newblock)
    self%next_id = self%next_id + 1
    newblock = cuda_field_t(self%ngrid, next, id=self%next_id)
    ptr => newblock
  end function create_cuda_block

end module m_cuda_allocator
