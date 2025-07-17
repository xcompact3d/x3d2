!!! test_vecadd
!!
!!  Tests the vecadd backend operation

program test_vecadd

  use MPI

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT
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

  integer, dimension(3), parameter :: dirs = [DIR_X, DIR_Y, DIR_Z]
  character(len=5), dimension(3), parameter :: dirnames = &
                         ["DIR_X", "DIR_Y", "DIR_Z"]

  class(mesh_t), allocatable :: mesh
  class(allocator_t), pointer :: allocator => null()
  class(base_backend_t), pointer :: backend => null()
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(allocator_t), target :: host_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(allocator_t), pointer :: host_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  class(field_t), pointer :: a => null()
  class(field_t), pointer :: b => null()
  class(field_t), pointer :: c => null()
  class(field_t), pointer :: z => null()

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
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
  host_allocator = allocator_t(mesh%get_dims(VERT), SZ)

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  test_pass = .true.

  do d = 1, size(dirs)

    call initialise_data(a, b, c, z, d)

    call test_add_zero(d)
    call test_double(d)
    call test_subself(d)
    call test_addab(d)

  end do

  call allocator%release_block(a)
  call allocator%release_block(b)
  call allocator%release_block(c)
  call allocator%release_block(z)

  if (test_pass) then
    print *, "PASS"
  else
    error stop "FAIL"
  end if

  call MPI_Finalize(ierr)

contains

  subroutine test_add_zero(d)

    integer, intent(in) :: d

    class(field_t), pointer :: r
    class(field_t), pointer :: r_host, a_host

    print *, dirnames(d), " Testing: a + 0 = a"

    r => allocator%get_block(dirs(d))
    r_host => host_allocator%get_block(dirs(d))
    a_host => host_allocator%get_block(dirs(d))

    call backend%get_field_data(r_host%data, z, dirs(d))
    call backend%set_field_data(r, r_host%data, dirs(d))

    call backend%vecadd(1.0_dp, a, 1.0_dp, r)

    call backend%get_field_data(a_host%data, a, dirs(d))
    call backend%get_field_data(r_host%data, r, dirs(d))
    if (any(r_host%data /= a_host%data)) then
      print *, dirnames(d), ": a + 0 = a failed"
      test_pass = .false.
    end if

    call allocator%release_block(r)
    call host_allocator%release_block(r_host)
    call host_allocator%release_block(a_host)

  end subroutine

  subroutine test_double(d)

    integer, intent(in) :: d

    class(field_t), pointer :: r
    class(field_t), pointer :: r_host, a_host

    print *, dirnames(d), " Testing: a + a = 2a"

    r => allocator%get_block(dirs(d))
    r_host => host_allocator%get_block(dirs(d))
    a_host => host_allocator%get_block(dirs(d))

    call backend%get_field_data(a_host%data, a, dirs(d))
    r_host%data = a_host%data
    call backend%set_field_data(r, r_host%data, dirs(d))

    call backend%vecadd(1.0_dp, a, 1.0_dp, r)

    call backend%get_field_data(r_host%data, r, dirs(d))
    if (any(r_host%data /= 2*a_host%data)) then
      print *, dirnames(d), ": a + a = 2a failed"
      test_pass = .false.
    end if

    call allocator%release_block(r)
    call host_allocator%release_block(r_host)
    call host_allocator%release_block(a_host)

  end subroutine

  subroutine test_subself(d)

    integer, intent(in) :: d

    class(field_t), pointer :: r
    class(field_t), pointer :: r_host, a_host

    print *, dirnames(d), " Testing: a + (-a) = 0"

    r => allocator%get_block(dirs(d))
    r_host => host_allocator%get_block(dirs(d))
    a_host => host_allocator%get_block(dirs(d))

    call backend%get_field_data(a_host%data, a, dirs(d))
    r_host%data = -a_host%data
    call backend%set_field_data(r, r_host%data, dirs(d))

    call backend%vecadd(1.0_dp, a, 1.0_dp, r)

    call backend%get_field_data(r_host%data, r, dirs(d))
    if (any(r_host%data /= 0)) then
      print *, dirnames(d), ": a + (-a) = 0 failed"
      test_pass = .false.
    end if

    call allocator%release_block(r)
    call host_allocator%release_block(r_host)
    call host_allocator%release_block(a_host)

  end subroutine

  subroutine test_addab(d)

    integer, intent(in) :: d

    class(field_t), pointer :: r
    class(field_t), pointer :: r_host, b_host, c_host

    print *, dirnames(d), " Testing: a + b = c"

    r => allocator%get_block(dirs(d))
    r_host => host_allocator%get_block(dirs(d))
    b_host => host_allocator%get_block(dirs(d))
    c_host => host_allocator%get_block(dirs(d))

    call backend%get_field_data(b_host%data, b, dirs(d))
    r_host%data = b_host%data
    call backend%set_field_data(r, r_host%data, dirs(d))

    call backend%vecadd(1.0_dp, a, 1.0_dp, r)

    call backend%get_field_data(r_host%data, r, dirs(d))
    call backend%get_field_data(c_host%data, c, dirs(d))
    if (any(r_host%data /= c_host%data)) then
      print *, dirnames(d), ": a + b = c failed"
      test_pass = .false.
    end if

    call allocator%release_block(r)
    call host_allocator%release_block(r_host)
    call host_allocator%release_block(b_host)
    call host_allocator%release_block(c_host)

  end subroutine

  subroutine initialise_data(a, b, c, z, d)

    integer, intent(in) :: d
    class(field_t), pointer, intent(inout) :: a, b, c, z
    class(field_t), pointer :: a_host, b_host, c_host, z_host

    print *, dirnames(d), ": Initialising test data"

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

    a => allocator%get_block(dirs(d))
    b => allocator%get_block(dirs(d))
    c => allocator%get_block(dirs(d))
    z => allocator%get_block(dirs(d))

    a_host => host_allocator%get_block(dirs(d))
    b_host => host_allocator%get_block(dirs(d))
    c_host => host_allocator%get_block(dirs(d))
    z_host => host_allocator%get_block(dirs(d))

    ! Initialise values
    call random_number(a_host%data)
    call random_number(b_host%data)
    c_host%data = a_host%data + b_host%data
    z_host%data = 0.0_dp

    call backend%set_field_data(a, a_host%data, dirs(d))
    call backend%set_field_data(b, b_host%data, dirs(d))
    call backend%set_field_data(c, c_host%data, dirs(d))
    call backend%set_field_data(z, z_host%data, dirs(d))

    call host_allocator%release_block(a_host)
    call host_allocator%release_block(b_host)
    call host_allocator%release_block(c_host)
    call host_allocator%release_block(z_host)

  end subroutine initialise_data

end program test_vecadd

