program test_hostside_io
  use iso_fortran_env, only: stderr => error_unit
  use adios2
  use mpi
  use adios_io

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, DIR_C
  use m_omp_common, only: SZ

  implicit none

  logical :: allpass !! flag indicating all tests pass

class(field_t), pointer :: arr_to_write
class(field_t), pointer :: arr_to_read

  !> MPI vars
  integer :: irank, nproc
  !> Error code, used for MPI and ADIOS
  integer :: ierr

  type(allocator_t), target :: omp_allocator
  type(adios_io_t) :: io

  integer(kind=8), dimension(3) :: ishape, istart, icount
  integer, parameter :: nx = 64, ny = 32, nz = 16

  integer :: i, j, k

  allpass = .true.

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  call io%init(MPI_COMM_WORLD)

  if (irank == 0) print*, 'Run with', nproc, 'ranks'

  !================ 
  ! SETUP TEST DATA
  !================

  icount = (/ nx, ny, nz/) ! local size
  ishape = (/ nproc*nx, ny, nz/) ! global size
  istart = (/ irank*nx, 0, 0/) ! local offset

  omp_allocator = allocator_t(nx, ny, nz, SZ)
  if (irank == 0) print*, 'OpenMP allocator instantiated'

  arr_to_write => omp_allocator%get_block(DIR_C)

  ! Initialise with a simple index to verify read later
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
  arr_to_write%data(i, j, k) = (istart(1) + i) + j*nx*nproc + k*nx*nproc*ny
  end do
  end do
  end do

  !================ 
  ! WRITE TEST DATA
  !================

  ! TODO check if adios can output filtered or coursened data
  ! TODO check if adios can output in a different precision
  ! TODO consider how to implement different mesh sizes for different variables (eg pressure & velocity)
  call io%write(arr_to_write, "test_omp_io.bp", "TestArr", icount, ishape, istart)
  call omp_allocator%release_block(arr_to_write)

  !================ 
  ! READ AND VERIFY TEST DATA
  !================

  arr_to_read => omp_allocator%get_block(DIR_C)
  call io%read(arr_to_read, "test_omp_io.bp", "TestArr", icount, ishape, istart)

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
  if(arr_to_read%data(i, j, k) /= (istart(1) + i) + j*nx*nproc + k*nx*nproc*ny) then
    write(stderr, '(a, f8.4, a, i8, a, i5, i5, i5)') &
      'Mismatch between read array(', arr_to_read%data(i, j, k), &
      ") and expected index (", (istart(1) + i) + j*nx*nproc + k*nx*nproc*ny, ") at (i,j,k) = ", i, j, k
    allpass = .false.
  end if
  end do
  end do
  end do

  call omp_allocator%release_block(arr_to_read)

  if (allpass) then
    if (irank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call io%deinit()
  call MPI_Finalize(ierr)

contains

  subroutine abort_test(msg)
    character(*), intent(in) :: msg
    integer :: ierr

    call io%deinit()
    call MPI_Finalize(ierr)

    error stop msg
  end subroutine
end program
