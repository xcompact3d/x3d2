
module adios_io
  use adios2
  use mpi
  use m_allocator, only: field_t

  implicit none

  type :: adios_io_t
    type(adios2_adios) :: adios_ctx
    integer :: irank, nproc
  contains
    procedure:: init 
    procedure:: deinit 
    procedure:: handle_fatal_error 
    procedure:: read 
    procedure:: write 
  end type

contains
  subroutine init(self)
  class(adios_io_t), intent(inout) :: self
    ! TODO pass in communicator?

    integer :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, self%irank, ierr)

    call adios2_init(self%adios_ctx, MPI_COMM_WORLD, ierr)
    if (.not.self%adios_ctx%valid) then
      print*, "ADIOS2 Error code:", ierr
      call self%handle_fatal_error("Cannot initialise ADIOS context.")
    endif
  end subroutine

  subroutine deinit(self)
  class(adios_io_t), intent(inout) :: self
    integer :: ierr

    if (self%adios_ctx%valid) then
      call adios2_finalize(self%adios_ctx, ierr)
    endif
  end subroutine

  subroutine write(self, in_arr, fpath, varname, icount, ishape, istart)
  class(adios_io_t), intent(inout) :: self
  class(field_t), pointer, intent(in) :: in_arr !! Field to be outputted
    character(*), intent(in) :: fpath !! Path to ouptut file
    character(*), intent(in) :: varname !! Name of variable in output file
    integer(kind=8), dimension(3), intent(in) :: icount !! Global size of in_arr
    integer(kind=8), dimension(3), intent(in) :: ishape !! Local size of in_arr
    integer(kind=8), dimension(3), intent(in) :: istart !! Local offset of in_arr

    type(adios2_io) :: io
    type(adios2_engine) :: writer
    type(adios2_variable) :: adios_var
    integer :: ierr

    call adios2_declare_io (io, self%adios_ctx, 'write', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO")
    endif

    call adios2_open(writer, io, fpath, adios2_mode_write, ierr)
    if (.not.writer%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 writer")
    endif

    ! TODO tie double precision variable to array
    call adios2_define_variable(adios_var, io, varname, adios2_type_dp, &
      3, ishape, istart, icount, adios2_constant_dims, ierr)

    call adios2_begin_step(writer, adios2_step_mode_append, ierr)
    call adios2_put(writer, adios_var, in_arr%data, ierr)
    call adios2_end_step(writer, ierr)

    if (writer%valid) then
      call adios2_close(writer, ierr)
    endif

  end subroutine

  subroutine read(self, out_arr, fpath, varname)
  class(adios_io_t), intent(inout) :: self
  class(field_t), pointer, intent(in) :: out_arr !! Field to be read from file
    character(*), intent(in) :: fpath !! Path to input file
    character(*), intent(in) :: varname !! Name of variable in input file

    integer(8), parameter :: initial_step = 0
    integer(8), parameter :: n_steps = 1

    type(adios2_io) :: io
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer :: ierr

    if (self%irank /= 0) then
      return
    endif

    call adios2_declare_io (io, self%adios_ctx, 'read', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO")
    endif

    call adios2_open(reader, io, fpath, adios2_mode_read, MPI_COMM_SELF, ierr)
    if (.not.reader%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 reader")
    endif

    call adios2_begin_step(reader, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(adios_var, io, varname, ierr)
    if (.not.adios_var%valid) then
      call self%handle_fatal_error("Cannot fetch ADIOS2 IO")
    endif

    call adios2_set_step_selection(adios_var, initial_step, n_steps, ierr)
    call adios2_get(reader, adios_var, out_arr%data, ierr)
    call adios2_end_step(reader, ierr)

    if (reader%valid) then
      call adios2_close(reader, ierr)
    endif
  end subroutine

  subroutine handle_fatal_error(self, msg)
  class(adios_io_t), intent(in) :: self
    character(*), intent(in) :: msg

    error stop msg
  end subroutine
end module


program test_omp_io
  use iso_fortran_env, only: stderr => error_unit
  use adios2
  use mpi
  use adios_io

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, DIR_C
  use m_omp_common, only: SZ

  implicit none

  logical :: allpass !! flag indicating all tests pass

class(field_t), pointer :: arr_to_write !! array to save & restore
class(field_t), pointer :: arr_to_read !! array to save & restore

  !> MPI vars
  integer :: irank, nproc
  !> Error code, used for MPI and ADIOS
  integer :: ierr

  type(allocator_t), target :: omp_allocator

  type(adios_io_t) :: io

  integer(kind=8), dimension(3) :: ishape, istart, icount
  integer :: n, n_block, n_glob
  integer, parameter :: nx = 64, ny = 32, nz = 16
  integer, parameter :: n_groups_x = 1

  integer :: i, j, k
  character(len=80)::fmt

  allpass = .true.

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  call io%init()

  ! if (nrank == 0 && nproc > 1) print*, 'Parallel run with', nproc, 'ranks'
  if (nproc > 1) call abort_test("Test does not support multiple MPI processes")

  !================ 
  ! SETUP TEST DATA
  !================

  icount = (/ nx, ny, nz/) ! global size
  ishape = (/ nx, ny, nz/) ! local size
  istart = (/ 0, 0, 0/) ! local offset

  omp_allocator = allocator_t(nx, ny, nz, SZ)
  print*, 'OpenMP allocator instantiated'

  arr_to_write => omp_allocator%get_block(DIR_C)

  ! Initialise with a simple index to verify read later
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
  arr_to_write%data(i, j, k) = i + j*nx + k*nx*ny
  end do
  end do
  end do

  call io%write(arr_to_write, "__FILE__.bp",  "TestArr", icount, ishape, istart)
  call omp_allocator%release_block(arr_to_write)

  arr_to_read => omp_allocator%get_block(DIR_C)
  call io%read(arr_to_read, "__FILE__.bp", "TestArr")

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
  if(arr_to_read%data(i, j, k) /= i + j*nx + k*nx*ny) then
    if (irank == 0) write(stderr, '(a, f8.4, a, f8.4, a, i5, i5, i5)') &
      'Mismatch between read array(', arr_to_read%data(i, j, k), &
      ") and expected index (", i + j*nx + k*nx*ny, "at (i,j,k) = ", i, j, k
    allpass = .false.
    return ! end test on first issue
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

