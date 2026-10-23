module m_alm_airfoil
  !! Static airfoil polars for the actuator line model (ALM).
  !!
  !! Ported from Incompact3d's `airfoils.f90` (itself taken from CACTUS). An
  !! airfoil holds one table of lift, drag and quarter-chord moment against
  !! angle of attack for each tabulated Reynolds number. The coefficients at a
  !! given (Re, alpha) are found by linear interpolation in alpha within the
  !! two bracketing Reynolds-number tables, then linearly in Re. Outside the
  !! tabulated Re range the nearest table is used.
  !!
  !! File format (`.air`): a title line containing ':', then the lines
  !! `Thickness to Chord Ratio: <t/c>`, `Zero Lift AOA (deg): <alpha_0>` and
  !! `Reverse Camber Direction: <0|1>`, followed by one block per Reynolds
  !! number: `Reynolds Number: <Re>`, a column-header line, and rows of
  !! `AOA[deg] CL CD CM25` terminated by a blank line or the end of file. Every
  !! table must span -180 to 180 degrees.
  use m_common, only: dp
  use m_logger, only: abort_run, log_warn

  implicit none

  private
  public :: airfoil_t, read_airfoil, eval_static_coeffs

  integer, parameter :: max_aoa_vals = 1000 !! AOA rows read per Re table
  integer, parameter :: max_re_vals = 10    !! Re tables read per airfoil

  type :: airfoil_t
    character(len=256) :: filename = ''
    character(len=256) :: title = ''
    integer :: camb = 0                !! reverse camber direction flag
    real(dp) :: tc = 0._dp             !! thickness to chord ratio
    real(dp) :: alzer = 0._dp          !! zero-lift angle of attack [deg]
    integer :: n_re = 0                !! number of Reynolds-number tables
    real(dp), allocatable :: re(:)     !! table Reynolds numbers
    integer, allocatable :: n_aoa(:)   !! AOA rows in each table
    real(dp), allocatable :: aoa(:, :) !! table AOA values [deg]
    real(dp), allocatable :: cl(:, :)
    real(dp), allocatable :: cd(:, :)
    real(dp), allocatable :: cm(:, :)
  end type airfoil_t

contains

  subroutine read_airfoil(airfoil, filename)
    !! Read an `.air` polar file into `airfoil`.
    implicit none
    type(airfoil_t), intent(out) :: airfoil
    character(len=*), intent(in) :: filename

    real(dp), allocatable :: ta(:, :), tcl(:, :), tcd(:, :), tcm(:, :)
    real(dp) :: tre(max_re_vals), row(4, max_aoa_vals)
    integer :: ntbl(max_re_vals)
    character(len=512) :: line
    integer :: unit, ios, ci, i, ii, jj, n_max
    logical :: eof

    airfoil%filename = filename
    open (newunit=unit, file=trim(filename), status='old', action='read', &
          iostat=ios)
    if (ios /= 0) call abort_run('ALM: cannot open airfoil file: ' &
                                 //trim(filename))

    ! Title block: the first line that contains ':'
    do
      read (unit, '(A)', iostat=ios) line
      if (ios /= 0) call abort_run('ALM: no title in airfoil file: ' &
                                   //trim(filename))
      ci = index(line, ':')
      if (ci > 0) exit
    end do
    if (len_trim(line) > ci) then
      airfoil%title = adjustl(line(ci + 1:len_trim(line)))
    else
      airfoil%title = 'No Title'
    end if
    call read_header_value(unit, filename, airfoil%tc)
    call read_header_value(unit, filename, airfoil%alzer)
    read (unit, '(A)', iostat=ios) line
    if (ios == 0) read (line(index(line, ':') + 1:), *, iostat=ios) &
      airfoil%camb
    if (ios /= 0) call abort_run('ALM: bad header in airfoil file: ' &
                                 //trim(filename))
    if (airfoil%camb == 1) airfoil%alzer = -airfoil%alzer

    allocate (ta(max_aoa_vals, max_re_vals), tcl(max_aoa_vals, max_re_vals), &
              tcd(max_aoa_vals, max_re_vals), tcm(max_aoa_vals, max_re_vals))

    ! First Reynolds-number block
    eof = .false.
    call find_next_block(unit, line, eof)

    i = 0
    do while (.not. eof .and. i < max_re_vals)
      i = i + 1
      read (line(index(line, ':') + 1:), *, iostat=ios) tre(i)
      if (ios /= 0) call abort_run('ALM: bad Reynolds number in airfoil &
                                   &file: '//trim(filename))
      read (unit, '(A)', iostat=ios) line ! column header
      ii = 0
      do
        read (unit, '(A)', iostat=ios) line
        if (ios /= 0) then
          eof = .true.
          exit
        end if
        if (is_blank(line)) exit
        if (ii == max_aoa_vals) then
          call log_warn('ALM: max. allowed AOA values exceeded in airfoil &
                        &data file: '//trim(airfoil%title))
          exit
        end if
        ii = ii + 1
        read (line, *, iostat=ios) ta(ii, i), tcl(ii, i), tcd(ii, i), &
          tcm(ii, i)
        if (ios /= 0) call abort_run('ALM: bad polar row in airfoil file: ' &
                                     //trim(filename))
      end do
      ntbl(i) = ii

      if (ii < 2) call abort_run('ALM: empty polar table in airfoil file: ' &
                                 //trim(filename))
      if (ta(1, i) > -180._dp .or. ta(ii, i) < 180._dp) then
        call abort_run('ALM: AOA data needs to be +/-180 deg in airfoil &
                       &data file: '//trim(airfoil%title))
      end if

      if (airfoil%camb == 1) then
        row(1, 1:ii) = ta(1:ii, i)
        row(2, 1:ii) = tcl(1:ii, i)
        row(3, 1:ii) = tcd(1:ii, i)
        row(4, 1:ii) = tcm(1:ii, i)
        do ii = 1, ntbl(i)
          jj = ntbl(i) - (ii - 1)
          ta(ii, i) = -row(1, jj)
          tcl(ii, i) = -row(2, jj)
          tcd(ii, i) = row(3, jj)
          tcm(ii, i) = -row(4, jj)
        end do
      end if

      if (.not. eof) call find_next_block(unit, line, eof)
    end do
    close (unit)

    if (i == 0) call abort_run('ALM: error reading airfoil data file: ' &
                               //trim(filename))
    if (.not. eof) call log_warn('ALM: max. allowed Re values exceeded in &
                                 &airfoil data file: '//trim(airfoil%title))

    airfoil%n_re = i
    n_max = maxval(ntbl(1:i))
    allocate (airfoil%re(i), airfoil%n_aoa(i))
    allocate (airfoil%aoa(n_max, i), airfoil%cl(n_max, i), &
              airfoil%cd(n_max, i), airfoil%cm(n_max, i))
    airfoil%re = tre(1:i)
    airfoil%n_aoa = ntbl(1:i)
    airfoil%aoa = ta(1:n_max, 1:i)
    airfoil%cl = tcl(1:n_max, 1:i)
    airfoil%cd = tcd(1:n_max, 1:i)
    airfoil%cm = tcm(1:n_max, 1:i)
  end subroutine read_airfoil

  pure subroutine eval_static_coeffs(airfoil, re, alpha, cl, cd, cm25)
    !! Interpolate the static polars at Reynolds number `re` and angle of
    !! attack `alpha` [deg]. Port of CACTUS' EvalStaticCoeff.
    implicit none
    type(airfoil_t), intent(in) :: airfoil
    real(dp), intent(in) :: re, alpha
    real(dp), intent(out) :: cl, cd, cm25

    real(dp) :: cla(2), cda(2), cma(2), xre, xa
    integer :: i_lb, i_ub, j, n, u1, l1, x1, m

    ! Bracketing Re tables
    xre = 0._dp
    if (re >= airfoil%re(1) .and. airfoil%n_re > 1) then
      i_ub = 2
      do
        if (re <= airfoil%re(i_ub)) then
          if (re == airfoil%re(i_ub)) then
            i_lb = i_ub
          else
            i_lb = i_ub - 1
            xre = (re - airfoil%re(i_lb))/(airfoil%re(i_ub) - airfoil%re(i_lb))
          end if
          exit
        else if (i_ub == airfoil%n_re) then
          ! Above the highest tabulated Re: use the last table
          i_lb = i_ub
          exit
        else
          i_ub = i_ub + 1
        end if
      end do
    else
      ! Below the lowest tabulated Re (or a single table): use the first
      i_lb = 1
      i_ub = 1
    end if

    cla = 0._dp
    cda = 0._dp
    cma = 0._dp
    m = 1
    do j = i_lb, i_ub
      n = airfoil%n_aoa(j)
      ! Interval halving on alpha
      u1 = n
      l1 = 1
      x1 = n/2
      do
        if (alpha < airfoil%aoa(x1, j)) then
          u1 = x1
        else
          l1 = x1
        end if
        if (u1 - l1 == 1) exit
        x1 = l1 + (u1 - l1)/2
      end do
      xa = (alpha - airfoil%aoa(l1, j)) &
           /(airfoil%aoa(u1, j) - airfoil%aoa(l1, j))
      cla(m) = airfoil%cl(l1, j) + xa*(airfoil%cl(u1, j) - airfoil%cl(l1, j))
      cda(m) = airfoil%cd(l1, j) + xa*(airfoil%cd(u1, j) - airfoil%cd(l1, j))
      cma(m) = airfoil%cm(l1, j) + xa*(airfoil%cm(u1, j) - airfoil%cm(l1, j))
      m = m + 1
    end do

    cl = cla(1) + xre*(cla(2) - cla(1))
    cd = cda(1) + xre*(cda(2) - cda(1))
    cm25 = cma(1) + xre*(cma(2) - cma(1))
  end subroutine eval_static_coeffs

  subroutine read_header_value(unit, filename, val)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: val

    character(len=512) :: line
    integer :: ios

    read (unit, '(A)', iostat=ios) line
    if (ios == 0) read (line(index(line, ':') + 1:), *, iostat=ios) val
    if (ios /= 0) call abort_run('ALM: bad header in airfoil file: ' &
                                 //trim(filename))
  end subroutine read_header_value

  subroutine find_next_block(unit, line, eof)
    !! Advance to the next line containing ':' (a Reynolds-number block).
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    logical, intent(out) :: eof

    integer :: ios

    eof = .false.
    do
      read (unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        eof = .true.
        return
      end if
      if (index(line, ':') > 0) return
    end do
  end subroutine find_next_block

  pure logical function is_blank(line)
    !! True for an empty line, also when it holds only a carriage return.
    character(len=*), intent(in) :: line

    is_blank = len_trim(line) == 0
    if (len_trim(line) == 1) is_blank = iachar(line(1:1)) == 13
  end function is_blank

end module m_alm_airfoil
