module m_alm_output
  !! Per-element load output of the actuator line model (ALM).
  !!
  !! Incompact3d writes one file per element; here each actuator line (blade
  !! or tower) gets one root-written CSV file, `<line name>.load`, with a row
  !! per element at every output step. Rotor performance goes through
  !! m_scalar_series in m_alm.
  use m_common, only: dp
  use m_logger, only: abort_run
  use m_alm_element, only: actuator_line_t, deg2rad

  implicit none

  private
  public :: init_loads_file, write_line_loads

contains

  subroutine init_loads_file(line, is_root, append)
    !! Start `<line name>.load` with its header, unless appending to an
    !! existing file on restart.
    implicit none
    type(actuator_line_t), intent(in) :: line
    logical, intent(in) :: is_root, append

    character(len=:), allocatable :: filename
    integer :: unit
    logical :: exists

    if (.not. is_root) return
    filename = trim(line%name)//'.load'
    inquire (file=filename, exist=exists)
    if (append .and. exists) return
    open (newunit=unit, file=filename, status='replace', action='write')
    write (unit, '(A)') '# iteration, time, element, x, y, z, r/R, AOA[deg], &
      &Re, Ur, CL, CD, CM25, Cn, Ct, Fn, Ft, F_end'
    close (unit)
  end subroutine init_loads_file

  subroutine write_line_loads(line, points, iter, time, is_root)
    !! Append the element loads of `line` at iteration `iter`. `points` are
    !! the element positions the loads were computed at, before the rotor
    !! advanced.
    implicit none
    type(actuator_line_t), intent(in) :: line
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: iter
    real(dp), intent(in) :: time
    logical, intent(in) :: is_root

    integer :: unit, ie, ios

    if (.not. is_root) return
    open (newunit=unit, file=trim(line%name)//'.load', status='old', &
          position='append', action='write', iostat=ios)
    if (ios /= 0) call abort_run('ALM: cannot open loads file for ' &
                                 //trim(line%name))
    do ie = 1, line%n_elem
      write (unit, '(I0,",",ES20.12,",",I0,15(",",ES20.12))') &
        iter, time, ie, points(:, ie), line%erdist(ie)/line%length, &
        line%eaoa(ie)/deg2rad, line%ere(ie), line%eur(ie), line%ecl(ie), &
        line%ecd(ie), line%ecm(ie), line%ecn(ie), line%ect(ie), &
        line%efn(ie), line%eft(ie), line%end_factor(ie)
    end do
    close (unit)
  end subroutine write_line_loads

end module m_alm_output
