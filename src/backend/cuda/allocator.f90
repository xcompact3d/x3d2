module m_cuda_allocator
  !! GPU memory allocator for CUDA backend.
  !!
  !! GPU memory (device memory) is physically separate from CPU memory (host).
  !! This allocator manages device-side storage, ensuring field data resides
  !! in GPU memory for kernel execution. Explicit device allocation avoids
  !! expensive implicit host-device transfers that would kill performance.
  !!
  !! **Design rationale:**
  !!
  !! - `cuda_field_t` extends `field_t` with device pointers (`p_data_d`, `data_d`)
  !! - Maintains both 1D and 3D views of same memory for flexibility
  !! - Reference counting prevents premature deallocation
  !! - Block-based allocation reduces allocation overhead
  !!
  use m_allocator, only: allocator_t
  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type, extends(allocator_t) :: cuda_allocator_t
    !! GPU memory allocator extending base allocator
  contains
    procedure :: create_block => create_cuda_block  !! Allocate GPU field block
  end type cuda_allocator_t

  interface cuda_allocator_t
    module procedure cuda_allocator_init
  end interface cuda_allocator_t

  type, extends(field_t) :: cuda_field_t
    !! Field residing in GPU device memory
    real(dp), device, pointer, private :: p_data_d(:)  !! 1D device memory pointer (raw allocation)
    real(dp), device, pointer, contiguous :: data_d(:, :, :)  !! 3D device view (for kernel access)
  contains
    procedure :: fill => fill_cuda               !! Fill with constant value
    procedure :: get_shape => get_shape_cuda     !! Query 3D dimensions
    procedure :: set_shape => set_shape_cuda     !! Reshape 3D view
  end type cuda_field_t

  interface cuda_field_t
    module procedure cuda_field_init
  end interface cuda_field_t

contains

  function cuda_field_init(ngrid, next, id) result(f)
    !! Initialise GPU field with device memory allocation.
    !!
    !! Device memory must be explicitly allocated before use. This constructor
    !! allocates the 1D device array and sets up metadata for later reshaping
    !! to 3D when dimensions are known.
    integer, intent(in) :: ngrid  !! Total number of grid points
    integer, intent(in) :: id     !! Unique field identifier
    type(cuda_field_t), pointer, intent(in) :: next  !! Next field in linked list
    type(cuda_field_t) :: f       !! Initialised field

    allocate (f%p_data_d(ngrid))
    f%refcount = 0
    f%next => next
    f%id = id
  end function cuda_field_init

  subroutine fill_cuda(self, c)
    !! Fill entire field with constant value on GPU.
    !!
    !! Initialising fields directly on GPU avoids transferring initialisation
    !! data from host. Single assignment to device array leverages GPU's
    !! memory controllers for efficient broadcast to all elements.
    implicit none

    class(cuda_field_t) :: self   !! Field to fill
    real(dp), intent(in) :: c     !! Constant value

    self%p_data_d = c

  end subroutine fill_cuda

  function get_shape_cuda(self) result(dims)
    !! Query current 3D dimensions of field.
    !!
    !! Fields are allocated with total size but reshaped dynamically based
    !! on decomposition. This query enables algorithms to adapt to actual
    !! current dimensions without hard-coding sizes.
    implicit none

    class(cuda_field_t) :: self   !! Field to query
    integer :: dims(3)            !! Current dimensions

    dims = shape(self%data_d)

  end function get_shape_cuda

  subroutine set_shape_cuda(self, dims)
    !! Reshape 3D view of device memory.
    !!
    !! Same 1D device allocation is reused for different pencil orientations
    !! (X-pencils, Y-pencils, Z-pencils). Reshaping avoids reallocating GPU
    !! memory, which is expensive. Fortran pointer remapping is essentially
    !! free, just changing metadata not data.
    implicit none

    class(cuda_field_t) :: self   !! Field to reshape
    integer, intent(in) :: dims(3)  !! New dimensions

    self%data_d(1:dims(1), 1:dims(2), 1:dims(3)) => self%p_data_d

  end subroutine set_shape_cuda

  function cuda_allocator_init(dims, sz) result(allocator)
    !! Initialise CUDA allocator with grid dimensions.
    !!
    !! Base allocator handles dimension calculations and block management
    !! logic. CUDA allocator only needs to override block creation to use
    !! device memory, avoiding code duplication.
    integer, intent(in) :: dims(3)  !! Grid dimensions
    integer, intent(in) :: sz       !! Pencil size (SZ)
    type(cuda_allocator_t) :: allocator  !! Initialised allocator

    allocator%allocator_t = allocator_t(dims, sz)
  end function cuda_allocator_init

  function create_cuda_block(self, next) result(ptr)
    !! Create new field block in GPU memory.
    !!
    !! Central allocation point ensures consistent initialisation and enables
    !! tracking (via IDs) for debugging memory issues. Returning base class
    !! pointer maintains polymorphism for generic algorithm code.
    class(cuda_allocator_t), intent(inout) :: self  !! Allocator instance
    type(cuda_field_t), pointer, intent(in) :: next  !! Next in linked list
    type(cuda_field_t), pointer :: newblock  !! Newly allocated block
    class(field_t), pointer :: ptr           !! Polymorphic return pointer
    allocate (newblock)
    self%next_id = self%next_id + 1
    newblock = cuda_field_t(self%ngrid, next, id=self%next_id)
    ptr => newblock
  end function create_cuda_block

end module m_cuda_allocator
