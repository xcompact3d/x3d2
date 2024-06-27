program test_scalar_product
  !! Tests the implementation of the scalar product function.
  !
  !! Given two fields, a and b, computes s = a_i * b_i where repeated indices
  !! imply summation.

  use m_common, only: DIR_X, DIR_Y, DIR_Z, DIR_C

  use m_allocator
  use m_base_backend
#ifdef CUDA
#else
  use m_omp_backend
  use m_omp_common, only: SZ
#endif

  implicit none

  integer, parameter :: nx = 17, ny = 32, nz = 59
  real(dp), parameter :: lx = 1.618, ly = 3.141529, lz = 1.729

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
#ifdef CUDA
#else
  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
#endif

  type(mesh_t) :: mesh

  character(len=5), dimension(4), parameter :: test = &
    ["DIR_X", "DIR_Y", "DIR_Z", "DIR_C"]
  integer, dimension(4), parameter :: dir = [DIR_X, DIR_Y, DIR_Z, DIR_C]
  integer :: i

  integer :: nrank, nproc
  integer :: ierr

  logical :: test_pass = .true.

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  mesh = mesh_t([nx, ny, nz], &
                [1, 1, nproc], &
                [lx, ly, lz])

#ifdef CUDA
#else
  omp_allocator = allocator_t(mesh, SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  do i = 1, 4
    call runtest(test(i), dir(i))
  end do

  if (nrank == 0) then
    if (.not. test_pass) then
      error stop "Test failed"
    end if
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call MPI_Finalize(ierr)

contains

  subroutine runtest(test, dir)

    character(len=*), intent(in) :: test
    integer, intent(in) :: dir

    class(field_t), pointer :: a, b, c
    real(dp) :: s

    integer :: n
    integer :: expt
    logical :: check_pass

    if (nrank == 0) then
      print *, "Testing ", test
    end if
    
    a => backend%allocator%get_block(dir)
    b => backend%allocator%get_block(dir)

    call a%set_data_loc(VERT)
    call b%set_data_loc(VERT)

    if (nrank == 0) then
      print *, "Simplest check: dot(0, 0) = 0"
    end if
    a%data = 0
    b%data = 0
    s = backend%scalar_product(a, b)
    if (s /= 0) then
      check_pass = .false.
    else
      check_pass = .true.
    end if
    call MPI_Allreduce(MPI_IN_PLACE, check_pass, 1, &
      MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
      ierr)
    if (nrank == 0) then
      print *, "- Got: ", s
      print *, "- Expected: ", 0
      if (.not. check_pass) then
        print *, "- FAIL"
      else
        print *, "- PASS"
      end if
    end if
    test_pass = test_pass .and. check_pass
      
    if (nrank == 0) then
      print *, "Check: dot(nrank, nrank) = sum^{nrank-1}_i=0 sum_n(i) i**2"
    end if
    a%data = (nrank + 1)
    s = backend%scalar_product(a, a)

    ! Determine number of interior points, using a temporary DIR_C field
    c => backend%allocator%get_block(DIR_C)
    call c%set_data_loc(a%data_loc)
    n = product(mesh%get_field_dims(c))
    call backend%allocator%release_block(c)
    
    expt = n*(nrank + 1)**2
    call MPI_Allreduce(MPI_IN_PLACE, expt, 1, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, &
                       ierr)
    if (s /= real(expt, kind(s))) then
      check_pass = .false.
    else
      check_pass = .true.
    end if
    call MPI_Allreduce(MPI_IN_PLACE, check_pass, 1, &
      MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
      ierr)
    if (nrank == 0) then
      print *, "- Got: ", s
      print *, "- Expected: ", expt
      if (.not. check_pass) then
        print *, "- FAIL"
      else
        print *, "- PASS"
      end if
    end if
    test_pass = test_pass .and. check_pass

    call backend%allocator%release_block(a)
    call backend%allocator%release_block(b)

  end subroutine runtest

end program test_scalar_product
