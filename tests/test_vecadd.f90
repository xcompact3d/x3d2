!!! test_vecadd
!!
!!  Tests the vecadd backend operation

program test_vecadd

  use MPI

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
#ifdef CUDA
  use m_cuda_common, only: SZ
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
#else
  use m_omp_common, only: SZ
  use m_omp_backend, only: omp_backend_t
#endif
  
  implicit none
  
  integer, dimension(4), parameter :: dirs = [ DIR_X, DIR_Y, DIR_Z, DIR_C ]
  character(len=5), dimension(4), parameter :: dirnames = [ "DIR_X", "DIR_Y", "DIR_Z", "DIR_C" ]

  class(mesh_t), allocatable :: mesh
  class(allocator_t), pointer :: allocator => null()
  class(base_backend_t), pointer :: backend => null()
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  class(field_t), pointer :: a => null()
  class(field_t), pointer :: b => null()
  class(field_t), pointer :: c => null()
  class(field_t), pointer :: z => null()
  class(field_t), pointer :: r => null()

  integer, dimension(3) :: dims_global
  integer, dimension(3) :: nproc_dir
  real(dp), dimension(3) :: L_global
  character(len=8), dimension(2) :: BC_x, BC_y, BC_z
  
  integer :: ierr

  logical :: test_pass
  integer :: d

  call MPI_Init(ierr)
  
  dims_global = [32, 32, 32]
  L_global = [1.0, 1.0, 1.0]
  nproc_dir = [1, 1, 1]
  BC_x = ["periodic", "periodic"]
  BC_y = BC_x
  BC_z = BC_x
  
  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)
#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh, SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  test_pass = .true.

  do d = 1, size(dirs)

    call initialise_data(a, b, c, z, r, d)
    
    ! a + 0 = a
    r%data = z%data
    call backend%vecadd(1.0_dp, a, 1.0_dp, r)
    if (any(r%data /= a%data)) then
      print *, dirnames(d), ": a + 0 = a FAILED"
      test_pass = .false.
    end if
    
    ! a + a = 2a
    r%data = a%data
    call backend%vecadd(1.0_dp, a, 1.0_dp, r)
    if (any(r%data /= 2 * a%data)) then
      print *, dirnames(d), ": a + a = 2a FAILED"
      test_pass = .false.
    end if

    ! a + (-a) = 0
    r%data = -a%data
    call backend%vecadd(1.0_dp, a, 1.0_dp, r)
    if (any(r%data /= 0.0_dp)) then
      print *, dirnames(d), ": a + (-a) = 0 FAILED"
      test_pass = .false.
    end if
    
    ! a + b = c
    r%data = b%data
    call backend%vecadd(1.0_dp, a, 1.0_dp, r)
    if (any(r%data /= c%data)) then
      print *, dirnames(d), ": a + b = c FAILED"
      test_pass = .false.
    end if

  end do

  call allocator%release_block(a)
  call allocator%release_block(b)
  call allocator%release_block(c)
  call allocator%release_block(z)
  call allocator%release_block(r)
  
  if (test_pass) then
    print *, "PASS"
  else
    error stop "FAIL"
  end if

  call MPI_Finalize(ierr)

contains

  subroutine initialise_data(a, b, c, z, r, d)

    integer, intent(in) :: d
    class(field_t), pointer, intent(inout) :: a, b, c, z, r

    if (associated(a)) then
      call allocator%release_block(a)
    end if
    if (associated(b)) then
      call allocator%release_block(b)
    end if
    if (associated(c)) then
      call allocator%release_block(c)
    end if
    if (associated(z)) then
      call allocator%release_block(z)
    end if
    if (associated(r)) then
      call allocator%release_block(r)
    end if

    a => allocator%get_block(dirs(d))
    b => allocator%get_block(dirs(d))
    c => allocator%get_block(dirs(d))
    z => allocator%get_block(dirs(d))
    r => allocator%get_block(dirs(d))

    ! Initialise values
    call random_number(a%data)
    call random_number(b%data)
    c%data = a%data + b%data
    z%data = 0.0_dp

  end subroutine initialise_data
  
end program test_vecadd
  
