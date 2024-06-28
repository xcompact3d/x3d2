program test_sum_intox
  !! Tests the implementation of summing a Y-oriented field into an X-oriented
  !! one.

  use m_common, only: DIR_X, DIR_Y, DIR_Z

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

  call runtest("YintoX", DIR_Y)
  call runtest("ZintoX", DIR_Z)

  if (nrank == 0) then
    if (.not. test_pass) then
      error stop "Test failed"
    end if
  end if

  call MPI_Finalize(ierr)

contains

  subroutine runtest(test, dir_from)

    use m_ordering, only: get_index_dir

    character(len=*), intent(in) :: test
    integer, intent(in) :: dir_from
    
    class(field_t), pointer :: a, b
    integer :: ctr
    integer :: i, j, k
    integer :: ii, jj, kk

    integer, dimension(3) :: dims
    logical :: check_pass

    if (nrank == 0) then
      print *, "Test ", test
    end if

    a => backend%allocator%get_block(DIR_X)
    b => backend%allocator%get_block(dir_from)

    dims = mesh%get_padded_dims(DIR_C)

    ! Initialise fields so that b = -a
    ctr = 0
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_dir(ii, jj, kk, i, j, k, DIR_X, SZ, &
                             dims(1), dims(2), dims(3))
          a%data(ii, jj, kk) = ctr
          call get_index_dir(ii, jj, kk, i, j, k, dir_from, SZ, &
                             dims(1), dims(2), dims(3))
          b%data(ii, jj, kk) = -ctr
          ctr = ctr + 1
        end do
      end do
    end do

    if (dir_from == DIR_Y) then
      call backend%sum_yintox(a, b)
    else
      call backend%sum_zintox(a, b)
    end if
    
    if ((minval(a%data) /= 0) .or. (maxval(a%data) /= 0)) then
      check_pass = .false.
    else
      check_pass = .true.
    end if
    call MPI_Allreduce(MPI_IN_PLACE, check_pass, 1, &
                       MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
                       ierr)
    test_pass = test_pass .and. check_pass

    if (nrank == 0) then
      if (check_pass) then
        print *, "- PASS"
      else
        print *, "- FAIL"
      end if
    end if

    call backend%allocator%release_block(a)
    call backend%allocator%release_block(b)

  end subroutine runtest

end program test_sum_intox
