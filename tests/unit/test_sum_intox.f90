program test_sum_intox
  !! Tests the implementation of summing a Y-oriented field into an X-oriented
  !! one.

  use m_common, only: DIR_X, DIR_Y, DIR_Z

  use m_allocator
  use m_base_backend
  use m_backend_runtime, only: backend_runtime_t, backend_sz
  use m_test_utils, only: initialise_mpi, finalise_test, global_all

  implicit none

  integer, parameter :: nx = 17, ny = 32, nz = 59
  real(dp), parameter :: lx = 1.618, ly = 3.141529, lz = 1.729

  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend

  type(mesh_t), target :: mesh

  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  integer :: nrank, nproc

  logical :: test_pass = .true.

  call initialise_mpi(nrank, nproc)

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t([nx, ny, nz], &
                [1, 1, nproc], &
                [lx, ly, lz], &
                BC_x, BC_y, BC_z)

  call runtime%init(mesh)
  backend => runtime%backend

  call runtest("YintoX", DIR_Y)
  call runtest("ZintoX", DIR_Z)

  call finalise_test(test_pass, nrank)

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

    dims = backend%allocator%get_padded_dims(DIR_C)

    ! Initialise fields so that b = -a
    ctr = 0
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_dir(ii, jj, kk, i, j, k, DIR_X, backend_sz, &
                             dims(1), dims(2), dims(3))
          a%data(ii, jj, kk) = ctr
          call get_index_dir(ii, jj, kk, i, j, k, dir_from, backend_sz, &
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

    check_pass = .not. ((minval(a%data) /= 0) .or. (maxval(a%data) /= 0))
    call global_all(check_pass)
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
