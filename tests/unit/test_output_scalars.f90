program test_output_scalars
  !! Unit tests for the scalar CSV output helper used by post-processing.
  !!
  !! These checks validate the file-writing behavior shared by scalar time
  !! series outputs, including the infrastructure used by monitoring.csv.
  use iso_fortran_env, only: stderr => error_unit, iostat_end
  use mpi

  use m_common, only: dp
  use m_scalar_series, only: scalar_series_t

  implicit none

  integer :: ierr, irank

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)

  if (irank == 0) then
    call test_root_writes_named_columns()
    call test_append_adds_rows_without_header()
    call test_non_root_is_noop()
    write (stderr, '(a)') 'Output scalar tests passed.'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine test_root_writes_named_columns()
    !! Verify that a root writer creates a fresh CSV file with the requested
    !! column names. The test writes two rows and reads them back to check both
    !! the header formatting and numeric values.
    type(scalar_series_t) :: series
    character(len=*), parameter :: filename = 'test_output_scalars_write.csv'
    character(len=128) :: header
    integer :: unit, ios
    real(dp) :: t, a, b, c

    call remove_file(filename)

    call series%init(filename, &
                     [character(len=16) :: 'lift', 'drag', 'moment'], &
                     is_root=.true., append=.false.)
    call series%write_step(0.0_dp, [1.25_dp, -2.5_dp, 3.75_dp])
    call series%write_step(0.5_dp, [4.0_dp, 5.5_dp, -6.25_dp])
    call series%finalise()

    open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
    call assert_iostat_ok(ios, 'open root output')

    read (unit, '(A)', iostat=ios) header
    call assert_iostat_ok(ios, 'read root header')
    call assert_equal(trim(header), '# time, lift, drag, moment', &
                      'root header')

    read (unit, *, iostat=ios) t, a, b, c
    call assert_iostat_ok(ios, 'read first root row')
    call assert_close(t, 0.0_dp, 'first row time')
    call assert_close(a, 1.25_dp, 'first row lift')
    call assert_close(b, -2.5_dp, 'first row drag')
    call assert_close(c, 3.75_dp, 'first row moment')

    read (unit, *, iostat=ios) t, a, b, c
    call assert_iostat_ok(ios, 'read second root row')
    call assert_close(t, 0.5_dp, 'second row time')
    call assert_close(a, 4.0_dp, 'second row lift')
    call assert_close(b, 5.5_dp, 'second row drag')
    call assert_close(c, -6.25_dp, 'second row moment')

    read (unit, *, iostat=ios)
    call assert_eof(ios, 'root output')
    close (unit)

    call remove_file(filename)
  end subroutine test_root_writes_named_columns

  subroutine test_append_adds_rows_without_header()
    !! Verify append mode for an existing scalar output file. The first writer
    !! creates the header and one row, then the second writer appends a row
    !! without emitting a duplicate header.
    type(scalar_series_t) :: series
    character(len=*), parameter :: filename = 'test_output_scalars_append.csv'
    character(len=128) :: header
    integer :: unit, ios
    real(dp) :: t, a, b

    call remove_file(filename)

    call series%init(filename, [character(len=16) :: 'ke', 'enstrophy'], &
                     is_root=.true., append=.false.)
    call series%write_step(1.0_dp, [10.0_dp, 20.0_dp])
    call series%finalise()

    call series%init(filename, [character(len=16) :: 'ke', 'enstrophy'], &
                     is_root=.true., append=.true.)
    call series%write_step(2.0_dp, [30.0_dp, 40.0_dp])
    call series%finalise()

    open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
    call assert_iostat_ok(ios, 'open append output')

    read (unit, '(A)', iostat=ios) header
    call assert_iostat_ok(ios, 'read append header')
    call assert_equal(trim(header), '# time, ke, enstrophy', 'append header')

    read (unit, *, iostat=ios) t, a, b
    call assert_iostat_ok(ios, 'read first append row')
    call assert_close(t, 1.0_dp, 'first append row time')
    call assert_close(a, 10.0_dp, 'first append row ke')
    call assert_close(b, 20.0_dp, 'first append row enstrophy')

    read (unit, *, iostat=ios) t, a, b
    call assert_iostat_ok(ios, 'read second append row')
    call assert_close(t, 2.0_dp, 'second append row time')
    call assert_close(a, 30.0_dp, 'second append row ke')
    call assert_close(b, 40.0_dp, 'second append row enstrophy')

    read (unit, *, iostat=ios)
    call assert_eof(ios, 'append output')
    close (unit)

    call remove_file(filename)
  end subroutine test_append_adds_rows_without_header

  subroutine test_non_root_is_noop()
    !! Verify that non-root ranks do not create or write scalar output files.
    !! This mirrors the intended root-only behavior used during MPI runs.
    type(scalar_series_t) :: series
    character(len=*), parameter :: filename = 'test_output_scalars_nonroot.csv'
    logical :: exists

    call remove_file(filename)

    call series%init(filename, [character(len=16) :: 'ignored'], &
                     is_root=.false., append=.false.)
    call series%write_step(1.0_dp, [99.0_dp])
    call series%finalise()

    inquire (file=filename, exist=exists)
    if (exists) then
      call remove_file(filename)
      error stop 'non-root scalar series writer created a file'
    end if
  end subroutine test_non_root_is_noop

  subroutine remove_file(filename)
    !! Delete a test output file if it already exists. This keeps each test
    !! independent of stale files from earlier failed or interrupted runs.
    character(len=*), intent(in) :: filename

    integer :: unit, ios
    logical :: exists

    inquire (file=filename, exist=exists)
    if (.not. exists) return

    open (newunit=unit, file=filename, status='old', action='readwrite', &
          iostat=ios)
    call assert_iostat_ok(ios, 'open file for cleanup')
    close (unit, status='delete')
  end subroutine remove_file

  subroutine assert_equal(actual, expected, label)
    !! Compare two strings and stop the test with a useful diagnostic if they
    !! differ. The label identifies which output field or line failed.
    character(len=*), intent(in) :: actual, expected, label

    if (actual /= expected) then
      write (stderr, '(a,": expected [",a,"], got [",a,"]")') &
        trim(label), expected, actual
      error stop 'scalar series string assertion failed'
    end if
  end subroutine assert_equal

  subroutine assert_close(actual, expected, label)
    !! Compare two floating-point values using a tight tolerance appropriate
    !! for values written and read in this unit test. The label makes failures
    !! point to the specific CSV row or column being checked.
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label

    real(dp), parameter :: tol = 1.0e-12_dp

    if (abs(actual - expected) > tol) then
      write (stderr, '(a,": expected ",ES20.12,", got ",ES20.12)') &
        trim(label), expected, actual
      error stop 'scalar series numeric assertion failed'
    end if
  end subroutine assert_close

  subroutine assert_iostat_ok(ios, label)
    !! Check that a file operation completed successfully. Non-zero iostat
    !! values are reported with context before stopping the test.
    integer, intent(in) :: ios
    character(len=*), intent(in) :: label

    if (ios /= 0) then
      write (stderr, '(a,": iostat=",i0)') trim(label), ios
      error stop 'scalar series file operation failed'
    end if
  end subroutine assert_iostat_ok

  subroutine assert_eof(ios, label)
    !! Check that reading past the expected rows reaches end-of-file. This
    !! catches accidental extra output, such as duplicate headers in append mode.
    integer, intent(in) :: ios
    character(len=*), intent(in) :: label

    if (ios /= iostat_end) then
      write (stderr, '(a,": expected EOF, iostat=",i0)') trim(label), ios
      error stop 'scalar series output had unexpected extra data'
    end if
  end subroutine assert_eof

end program test_output_scalars
