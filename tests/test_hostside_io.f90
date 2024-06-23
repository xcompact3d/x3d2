program test_hostside_io
  use iso_fortran_env, only: stderr => error_unit
  use adios2
  use mpi
  use m_adios_io

  implicit none

  !> MPI vars
  integer :: irank, nproc
  !> Error code, used for MPI and ADIOS
  integer :: ierr

  call run_tests()

contains

  subroutine run_tests()
    logical :: allpass = .true. !! flag indicating all tests pass

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (irank == 0) print*, 'Run with', nproc, 'ranks'

    if(.not. test_basic_io()) then
      write(stderr, '(a)') 'Basic IO test failed.'
      allpass = .false.
    endif
    if(.not. test_single_precision_io()) then
      write(stderr, '(a)') 'Single precision IO test failed.'
      allpass = .false.
    endif
    if(.not. test_strided_io()) then
      write(stderr, '(a)') 'Strided IO test failed.'
      allpass = .false.
    endif
    if(.not. test_strided_single_precision_io()) then
      write(stderr, '(a)') 'Strided + SP IO test failed.'
      allpass = .false.
    endif

    if (allpass) then
      if (irank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    end if

    call MPI_Finalize(ierr)

    if(.not. allpass) then
      write(stderr, '(a)') 'SOME TESTS FAILED.'
      stop
    endif
  end subroutine

  function test_basic_io() result(test_passed)
    use m_allocator, only: allocator_t, field_t
    use m_common, only: dp, DIR_C
    use m_omp_common, only: SZ

    class(field_t), pointer :: arr_to_write
    class(field_t), pointer :: arr_to_read

    !> Did the test pass?
    logical :: test_passed

    type(allocator_t), target :: omp_allocator
    type(adios_io_t) :: io

    integer(kind=8), dimension(3) :: ishape, istart, icount
    integer, parameter :: nx = 64, ny = 32, nz = 16

    integer :: i, j, k

    test_passed = .true.

    call io%init(MPI_COMM_WORLD)

    if(irank == 0) write(stderr, '(a)') 'Starting basic IO test.'

    !================ 
    ! SETUP TEST DATA
    !================

    icount = (/ nx, ny, nz/) ! local size
    ishape = (/ nproc*nx, ny, nz/) ! global size
    istart = (/ irank*nx, 0, 0/) ! local offset

    omp_allocator = allocator_t(nx, ny, nz, SZ)

    arr_to_write => omp_allocator%get_block(DIR_C)

    ! Initialise with a simple index to verify read later
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    arr_to_write%data(i, j, k) = (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny
    end do
    end do
    end do

    !================ 
    ! WRITE TEST DATA
    !================

    call io%write(arr_to_write, "test_basic_io.bp", "TestArr", icount, ishape, istart)
    call omp_allocator%release_block(arr_to_write)

    !================ 
    ! READ AND VERIFY TEST DATA
    !================

    arr_to_read => omp_allocator%get_block(DIR_C)
    call io%read(arr_to_read, "test_basic_io.bp", "TestArr", icount, ishape, istart)

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    if(arr_to_read%data(i, j, k) /= (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny) then
      write(stderr, '(a, f8.4, a, i8, a, i5, i5, i5)') &
        'Mismatch between read array(', arr_to_read%data(i, j, k), &
        ") and expected index (", (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny, ") at (i,j,k) = ", i, j, k
      test_passed = .false.
    end if
    end do
    end do
    end do

    call omp_allocator%release_block(arr_to_read)

    call io%deinit()

    if(irank == 0) write(stderr, '(a)') 'Finishing basic IO test.'
  end function

  function test_single_precision_io() result(test_passed)
    use m_allocator, only: allocator_t, field_t
    use m_common, only: dp, DIR_C
    use m_omp_common, only: SZ

    class(field_t), pointer :: arr_to_write

    !> Did the test pass?
    logical :: test_passed

    type(allocator_t), target :: omp_allocator
    type(adios_io_t) :: io

    type(adios2_io) :: aio
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer(kind=8) :: idump = 0

    integer(kind=8), dimension(3) :: ishape, istart, icount
    integer, parameter :: nx = 64, ny = 32, nz = 16

    real(kind=4) :: arr_to_read(nx, ny, nz)

    integer :: i, j, k

    test_passed = .true.

    call io%init(MPI_COMM_WORLD)

    if(irank == 0) write(stderr, '(a)') 'Starting single precision IO test.'

    !================ 
    ! SETUP TEST DATA
    !================

    icount = (/ nx, ny, nz/) ! local size
    ishape = (/ nproc*nx, ny, nz/) ! global size
    istart = (/ irank*nx, 0, 0/) ! local offset

    omp_allocator = allocator_t(nx, ny, nz, SZ)

    arr_to_write => omp_allocator%get_block(DIR_C)

    ! Initialise with a simple index to verify read later
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    arr_to_write%data(i, j, k) = (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny
    end do
    end do
    end do

    !================ 
    ! WRITE TEST DATA
    !================

    call io%write(arr_to_write, "test_single_precision_io_data.bp", "TestArr", icount, ishape, istart, .true.)
    call omp_allocator%release_block(arr_to_write)

    !================ 
    ! READ AND VERIFY TEST DATA
    !================

    call adios2_declare_io (aio, io%adios_ctx, 'read', ierr)
    if (.not. aio%valid) then
      write(stderr, *) "Cannot create ADIOS2 IO"
      test_passed = .false.
      return
    endif
    call adios2_open(reader, aio, "test_single_precision_io_data.bp", adios2_mode_read, MPI_COMM_SELF, ierr)
    if (.not. reader%valid) then
      write(stderr, *) "Cannot create ADIOS2 reader"
      test_passed = .false.
      return
    endif
    call adios2_begin_step(reader, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(adios_var, aio, "TestArr", ierr)
    if (.not. adios_var%valid) then
      write(stderr, *) "Cannot fetch ADIOS2 variable"
      test_passed = .false.
      return
    endif

    call adios2_set_step_selection(adios_var, 0_8, 1_8, ierr)
    call adios2_set_selection(adios_var, 3, istart, icount, ierr)
    call adios2_get(reader, adios_var, arr_to_read, ierr)
    call adios2_end_step(reader, ierr)

    if (reader%valid) then
      call adios2_close(reader, ierr)
    endif

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    if(arr_to_read(i, j, k) /= (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny) then
      write(stderr, '(a, f8.4, a, i8, a, i5, i5, i5)') &
        'Mismatch between read array(', arr_to_read(i, j, k), &
        ") and expected index (", (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny, ") at (i,j,k) = ", i, j, k
      test_passed = .false.
      stop
    end if
    end do
    end do
    end do

    call io%deinit()

    if(irank == 0) write(stderr, '(a)') 'Finishing single precision IO test.'
  end function

  function test_strided_io() result(test_passed)
    use m_allocator, only: allocator_t, field_t
    use m_common, only: dp, DIR_C
    use m_omp_common, only: SZ

    class(field_t), pointer :: arr_to_write
    real(kind=8), dimension(:,:,:), allocatable :: arr_to_read

    !> Did the test pass?
    logical :: test_passed

    type(allocator_t), target :: omp_allocator
    type(adios_io_t) :: io

    type(adios2_io) :: aio
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer(kind=8) :: idump = 0

    integer(kind=8), dimension(3) :: ishape, istart, icount, istride, icount_strided, ishape_strided
    integer, parameter :: nx = 32, ny = 16, nz = 64

    integer :: i, j, k, idx

    test_passed = .true.

    call io%init(MPI_COMM_WORLD)

    if(irank == 0) write(stderr, '(a)') 'Starting strided IO test.'

    !================ 
    ! SETUP TEST DATA
    !================

    icount = (/ nx, ny, nz/) ! local size
    ishape = (/ nproc*nx, ny, nz/) ! global size
    istart = (/ irank*nx, 0, 0/) ! local offset

    omp_allocator = allocator_t(nx, ny, nz, SZ)

    arr_to_write => omp_allocator%get_block(DIR_C)

    ! Initialise with a simple index to verify read later
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    arr_to_write%data(i, j, k) = (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny
    end do
    end do
    end do

    !================ 
    ! WRITE TEST DATA
    !================

    istride = [2_8, 3_8, 4_8]
    icount_strided = icount / istride

    call io%write(arr_to_write, "test_strided_io_data.bp", "TestArr", icount, ishape, istart, istride_in=istride)
    call omp_allocator%release_block(arr_to_write)

    !================ 
    ! READ AND VERIFY TEST DATA
    !================

    ! Read only from one rank; let's not overcomplicate this
    if(irank == 0) then
      ishape_strided = ishape/istride
      allocate(arr_to_read(ishape_strided(1), ishape_strided(2), ishape_strided(3)))

      call adios2_declare_io (aio, io%adios_ctx, 'read', ierr)
      if (.not. aio%valid) then
        write(stderr, *) "Cannot create ADIOS2 IO"
        test_passed = .false.
        return
      endif
      call adios2_open(reader, aio, "test_strided_io_data.bp", adios2_mode_read, MPI_COMM_SELF, ierr)
      if (.not. reader%valid) then
        write(stderr, *) "Cannot create ADIOS2 reader"
        test_passed = .false.
        return
      endif
      call adios2_begin_step(reader, adios2_step_mode_read, ierr)
      call adios2_inquire_variable(adios_var, aio, "TestArr", ierr)
      if (.not. adios_var%valid) then
        write(stderr, *) "Cannot fetch ADIOS2 variable"
        test_passed = .false.
        return
      endif

      call adios2_set_step_selection(adios_var, 0_8, 1_8, ierr)
      call adios2_set_selection(adios_var, 3, [0_8,0_8,0_8], ishape_strided, ierr)
      call adios2_get(reader, adios_var, arr_to_read, ierr)
      call adios2_end_step(reader, ierr)

      if (reader%valid) then
        call adios2_close(reader, ierr)
      endif

      do k = 1, ishape_strided(3)
      do j = 1, ishape_strided(2)
      do i = 1, ishape_strided(1)
      idx = (istart(1) + i-1)*istride(1) + (j-1)*istride(2)*nx*nproc + (k-1)*istride(3)*nx*nproc*ny
      if(arr_to_read(i, j, k) /= idx) then
        write(stderr, '(a, f8.4, a, i8, a, i5, i5, i5)') &
          'Mismatch between read array(', arr_to_read(i, j, k), &
          ") and expected index (", idx, ") at (i,j,k) = ", i, j, k
        test_passed = .false.
      end if
      end do
      end do
      end do

      deallocate(arr_to_read)
    endif

    call io%deinit()

    if(irank == 0) write(stderr, '(a)') 'Finishing strided IO test.'
  end function

  function test_strided_single_precision_io() result(test_passed)
    use m_allocator, only: allocator_t, field_t
    use m_common, only: dp, DIR_C
    use m_omp_common, only: SZ

    class(field_t), pointer :: arr_to_write
    real(kind=4), dimension(:,:,:), allocatable :: arr_to_read

    !> Did the test pass?
    logical :: test_passed

    type(allocator_t), target :: omp_allocator
    type(adios_io_t) :: io

    type(adios2_io) :: aio
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer(kind=8) :: idump = 0

    integer(kind=8), dimension(3) :: ishape, istart, icount, istride, icount_strided, ishape_strided
    integer, parameter :: nx = 32, ny = 16, nz = 64

    integer :: i, j, k, idx

    test_passed = .true.

    call io%init(MPI_COMM_WORLD)

    if(irank == 0) write(stderr, '(a)') 'Starting strided IO test.'

    !================ 
    ! SETUP TEST DATA
    !================

    icount = (/ nx, ny, nz/) ! local size
    ishape = (/ nproc*nx, ny, nz/) ! global size
    istart = (/ irank*nx, 0, 0/) ! local offset

    omp_allocator = allocator_t(nx, ny, nz, SZ)

    arr_to_write => omp_allocator%get_block(DIR_C)

    ! Initialise with a simple index to verify read later
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    arr_to_write%data(i, j, k) = (istart(1) + i-1) + (j-1)*nx*nproc + (k-1)*nx*nproc*ny
    end do
    end do
    end do

    !================ 
    ! WRITE TEST DATA
    !================

    istride = [2_8, 3_8, 4_8]
    icount_strided = icount / istride

    call io%write(arr_to_write, "test_strided_io_data.bp", "TestArr", icount, ishape, istart, convert_to_sp_in=.true., istride_in=istride)
    call omp_allocator%release_block(arr_to_write)

    !================ 
    ! READ AND VERIFY TEST DATA
    !================

    ! Read only from one rank; let's not overcomplicate this
    if(irank == 0) then
      ishape_strided = ishape/istride
      allocate(arr_to_read(ishape_strided(1), ishape_strided(2), ishape_strided(3)))

      call adios2_declare_io (aio, io%adios_ctx, 'read', ierr)
      if (.not. aio%valid) then
        write(stderr, *) "Cannot create ADIOS2 IO"
        test_passed = .false.
        return
      endif
      call adios2_open(reader, aio, "test_strided_io_data.bp", adios2_mode_read, MPI_COMM_SELF, ierr)
      if (.not. reader%valid) then
        write(stderr, *) "Cannot create ADIOS2 reader"
        test_passed = .false.
        return
      endif
      call adios2_begin_step(reader, adios2_step_mode_read, ierr)
      call adios2_inquire_variable(adios_var, aio, "TestArr", ierr)
      if (.not. adios_var%valid) then
        write(stderr, *) "Cannot fetch ADIOS2 variable"
        test_passed = .false.
        return
      endif

      call adios2_set_step_selection(adios_var, 0_8, 1_8, ierr)
      call adios2_set_selection(adios_var, 3, [0_8,0_8,0_8], ishape_strided, ierr)
      call adios2_get(reader, adios_var, arr_to_read, ierr)
      call adios2_end_step(reader, ierr)

      if (reader%valid) then
        call adios2_close(reader, ierr)
      endif

      do k = 1, ishape_strided(3)
      do j = 1, ishape_strided(2)
      do i = 1, ishape_strided(1)
      idx = (istart(1) + i-1)*istride(1) + (j-1)*istride(2)*nx*nproc + (k-1)*istride(3)*nx*nproc*ny
      if(arr_to_read(i, j, k) /= idx) then
        write(stderr, '(a, f8.4, a, i8, a, i5, i5, i5)') &
          'Mismatch between read array(', arr_to_read(i, j, k), &
          ") and expected index (", idx, ") at (i,j,k) = ", i, j, k
        test_passed = .false.
      end if
      end do
      end do
      end do

      deallocate(arr_to_read)
    endif

    call io%deinit()

    if(irank == 0) write(stderr, '(a)') 'Finishing strided IO test.'
  end function

end program
