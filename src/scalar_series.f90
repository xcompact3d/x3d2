module m_scalar_series
  !! Root-only writer for named scalar time series.
  use m_common, only: dp

  implicit none

  private
  public :: scalar_series_t

  type :: scalar_series_t
    integer, private :: file_unit = -1
    integer, private :: n_columns = 0
    logical, private :: is_root = .false.
  contains
    procedure :: init
    procedure :: write_step
    procedure :: finalise
  end type scalar_series_t

contains

  subroutine init(self, filename, column_names, is_root, append)
    class(scalar_series_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: column_names(:)
    logical, intent(in) :: is_root, append

    integer :: i
    logical :: exists

    self%is_root = is_root
    self%n_columns = size(column_names)
    if (.not. is_root) return

    inquire (file=filename, exist=exists)
    if (append .and. exists) then
      open (newunit=self%file_unit, file=filename, status='old', &
            position='append', action='write')
      return
    end if

    open (newunit=self%file_unit, file=filename, status='replace', &
          action='write')
    write (self%file_unit, '(A)', advance='no') '# time'
    do i = 1, size(column_names)
      write (self%file_unit, '(A)', advance='no') ', '//trim(column_names(i))
    end do
    write (self%file_unit, *)
  end subroutine init

  subroutine write_step(self, t, values)
    class(scalar_series_t), intent(inout) :: self
    real(dp), intent(in) :: t
    real(dp), intent(in) :: values(:)

    integer :: i

    if (.not. self%is_root) return
    if (size(values) /= self%n_columns) then
      error stop 'Scalar series value count does not match its header.'
    end if
    write (self%file_unit, '(ES20.12)', advance='no') t
    do i = 1, size(values)
      write (self%file_unit, '(",",ES20.12)', advance='no') values(i)
    end do
    write (self%file_unit, *)
    flush (self%file_unit)
  end subroutine write_step

  subroutine finalise(self)
    class(scalar_series_t), intent(inout) :: self

    if (self%is_root .and. self%file_unit /= -1) close (self%file_unit)
    self%file_unit = -1
    self%n_columns = 0
  end subroutine finalise

end module m_scalar_series
