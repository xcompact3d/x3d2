module m_alm_element
  !! Actuator lines for the actuator line model (ALM): blade-element geometry,
  !! rigid-body kinematics and aerodynamic loads.
  !!
  !! Ported from Incompact3d's `acl_elem.f90` and `acl_utils.f90`. An actuator
  !! line is defined by its stations (element ends) on the quarter-chord line;
  !! element centres, unit vectors (spanwise s, chordwise t, normal n), areas
  !! and chords are derived from them by make_line_geometry. Loads come from
  !! the static airfoil polars,
  !! \[
  !! \alpha = \operatorname{atan2}(u_n, u_t),\qquad
  !! \mathbf{F} = \tfrac12\rho\,U_r^2 A\,F_{end}
  !!   \left(C_N\mathbf{n} + C_T\mathbf{t}\right),
  !! \]
  !! with \(C_N = C_L\cos\alpha + C_D\sin\alpha\),
  !! \(C_T = -C_L\sin\alpha + C_D\cos\alpha\) and \(u_n, u_t\) the normal and
  !! chordwise components of the fluid velocity relative to the element.
  !!
  !! The dynamic-stall, added-mass, random-walk and constant-circulation
  !! options of Incompact3d are not ported.
  use m_common, only: dp, pi
  use m_alm_airfoil, only: airfoil_t, eval_static_coeffs

  implicit none

  private
  public :: actuator_line_t, allocate_line, make_line_geometry, &
            rotate_line, pitch_line, compute_line_forces, &
            compute_tower_forces, populate_line_airfoils, &
            read_line_geometry, quat_rot, cross

  real(dp), parameter, public :: deg2rad = pi/180._dp

  type :: actuator_line_t
    character(len=256) :: name = ''
    integer :: n_elem = 0
    logical :: flip_n = .false.        !! flip the normal (clockwise rotor)
    real(dp) :: cor(3) = 0._dp         !! centre of rotation
    real(dp) :: length = 0._dp         !! reference length (blade radius)

    ! Station (element end) values
    real(dp), allocatable :: qc(:, :)   !! (3, n_elem + 1) quarter-chord points
    real(dp), allocatable :: t(:, :)    !! (3, n_elem + 1) chordwise unit vec
    real(dp), allocatable :: chord(:)
    real(dp), allocatable :: thick(:)   !! thickness to chord ratio
    real(dp), allocatable :: pitch(:)   !! station pitch [rad]

    ! Element values
    real(dp), allocatable :: pe(:, :)   !! (3, n_elem) element centres
    real(dp), allocatable :: te(:, :)   !! chordwise unit vector
    real(dp), allocatable :: ne(:, :)   !! normal unit vector
    real(dp), allocatable :: se(:, :)   !! spanwise unit vector
    real(dp), allocatable :: ec(:)      !! chord
    real(dp), allocatable :: eds(:)     !! spanwise length
    real(dp), allocatable :: earea(:)   !! planform area
    real(dp), allocatable :: erdist(:)  !! distance from cor
    real(dp), allocatable :: ettoc(:)   !! thickness to chord ratio
    real(dp), allocatable :: epitch(:)  !! pitch [rad]
    integer, allocatable :: eairfoil(:) !! index into the owner's airfoils

    ! Flow and body velocities at the elements
    real(dp), allocatable :: evel(:, :) !! (3, n_elem) sampled fluid velocity
    real(dp), allocatable :: evb(:, :)  !! (3, n_elem) body velocity

    ! Element aerodynamics
    real(dp), allocatable :: eaoa(:)    !! angle of attack [rad]
    real(dp), allocatable :: ere(:)     !! chord Reynolds number
    real(dp), allocatable :: eur(:)     !! relative speed
    real(dp), allocatable :: ecl(:), ecd(:), ecm(:), ecn(:), ect(:)
    real(dp), allocatable :: efn(:), eft(:), ems(:)
    real(dp), allocatable :: ef(:, :)   !! (3, n_elem) force [N]
    real(dp), allocatable :: end_factor(:) !! tip/root loss factor
  end type actuator_line_t

contains

  subroutine allocate_line(line, n_stations)
    !! (Re)allocate the station and element arrays of `line`.
    implicit none
    type(actuator_line_t), intent(inout) :: line
    integer, intent(in) :: n_stations

    integer :: ne

    ne = n_stations - 1
    line%n_elem = ne
    call realloc2(line%qc, n_stations)
    call realloc2(line%t, n_stations)
    call realloc1(line%chord, n_stations)
    call realloc1(line%thick, n_stations)
    call realloc1(line%pitch, n_stations)
    call realloc2(line%pe, ne)
    call realloc2(line%te, ne)
    call realloc2(line%ne, ne)
    call realloc2(line%se, ne)
    call realloc1(line%ec, ne)
    call realloc1(line%eds, ne)
    call realloc1(line%earea, ne)
    call realloc1(line%erdist, ne)
    call realloc1(line%ettoc, ne)
    call realloc1(line%epitch, ne)
    call realloc2(line%evel, ne)
    call realloc2(line%evb, ne)
    call realloc1(line%eaoa, ne)
    call realloc1(line%ere, ne)
    call realloc1(line%eur, ne)
    call realloc1(line%ecl, ne)
    call realloc1(line%ecd, ne)
    call realloc1(line%ecm, ne)
    call realloc1(line%ecn, ne)
    call realloc1(line%ect, ne)
    call realloc1(line%efn, ne)
    call realloc1(line%eft, ne)
    call realloc1(line%ems, ne)
    call realloc2(line%ef, ne)
    call realloc1(line%end_factor, ne)
    if (allocated(line%eairfoil)) deallocate (line%eairfoil)
    allocate (line%eairfoil(ne))
    line%eairfoil = 1
    line%end_factor = 1._dp

  contains

    subroutine realloc1(a, n)
      real(dp), allocatable, intent(inout) :: a(:)
      integer, intent(in) :: n
      if (allocated(a)) deallocate (a)
      allocate (a(n))
      a = 0._dp
    end subroutine realloc1

    subroutine realloc2(a, n)
      real(dp), allocatable, intent(inout) :: a(:, :)
      integer, intent(in) :: n
      if (allocated(a)) deallocate (a)
      allocate (a(3, n))
      a = 0._dp
    end subroutine realloc2
  end subroutine allocate_line

  subroutine read_line_geometry(filename, length, span_vec, r_over_l, &
                                c_over_l, pitch, thick)
    !! Read an actuator line geometry (`.al`) file:
    !!   `<label>: <length>`, `<label>: <sx> <sy> <sz>`, `<label>: <n>`,
    !! then n rows of `r/L  c/L  pitch[deg]  t/c`.
    use m_logger, only: abort_run
    implicit none
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: length, span_vec(3)
    real(dp), allocatable, intent(out) :: r_over_l(:), c_over_l(:), &
                                          pitch(:), thick(:)

    character(len=1024) :: line
    integer :: unit, ios, n, i

    open (newunit=unit, file=trim(filename), status='old', action='read', &
          iostat=ios)
    if (ios /= 0) call abort_run('ALM: cannot open geometry file: ' &
                                 //trim(filename))
    read (unit, '(A)', iostat=ios) line
    if (ios == 0) read (line(index(line, ':') + 1:), *, iostat=ios) length
    if (ios == 0) read (unit, '(A)', iostat=ios) line
    if (ios == 0) read (line(index(line, ':') + 1:), *, iostat=ios) span_vec
    if (ios == 0) read (unit, '(A)', iostat=ios) line
    if (ios == 0) read (line(index(line, ':') + 1:), *, iostat=ios) n
    if (ios /= 0) call abort_run('ALM: bad header in geometry file: ' &
                                 //trim(filename))
    if (n < 2) call abort_run('ALM: geometry file needs at least two &
                              &stations: '//trim(filename))

    allocate (r_over_l(n), c_over_l(n), pitch(n), thick(n))
    do i = 1, n
      read (unit, *, iostat=ios) r_over_l(i), c_over_l(i), pitch(i), thick(i)
      if (ios /= 0) call abort_run('ALM: bad station row in geometry file: ' &
                                   //trim(filename))
    end do
    close (unit)
  end subroutine read_line_geometry

  subroutine make_line_geometry(line)
    !! Element centres, unit vectors, areas and chords from the stations.
    implicit none
    type(actuator_line_t), intent(inout) :: line

    real(dp) :: s_e(3), t_e(3), n_e(3), p1(3), p2(3), p3(3), p4(3)
    real(dp) :: a1(3), a2(3), s_mag
    integer :: j

    do j = 1, line%n_elem
      line%pe(:, j) = (line%qc(:, j + 1) + line%qc(:, j))/2._dp
      line%erdist(j) = norm2(line%pe(:, j) - line%cor)

      ! Spanwise along the quarter-chord line
      s_e = line%qc(:, j + 1) - line%qc(:, j)
      s_mag = sqrt(dot_product(s_e, s_e))
      line%eds(j) = s_mag
      s_e = s_e/s_mag

      ! Mean chordwise vector, made normal to the spanwise one
      t_e = (line%t(:, j + 1) + line%t(:, j))/2._dp
      t_e = t_e - dot_product(t_e, s_e)*s_e
      t_e = t_e/sqrt(dot_product(t_e, t_e))

      n_e = cross(s_e, t_e)
      n_e = n_e/sqrt(dot_product(n_e, n_e))

      line%se(:, j) = s_e
      line%te(:, j) = t_e
      line%ne(:, j) = n_e
      if (line%flip_n) then
        line%ne(:, j) = -n_e
        line%se(:, j) = -s_e
      end if

      ! Element planform: quad from the leading and trailing edges of both
      ! stations, split into two triangles
      p1 = line%qc(:, j) - 0.25_dp*line%chord(j)*line%t(:, j)
      p2 = line%qc(:, j) + 0.75_dp*line%chord(j)*line%t(:, j)
      p3 = line%qc(:, j + 1) + 0.75_dp*line%chord(j + 1)*line%t(:, j + 1)
      p4 = line%qc(:, j + 1) - 0.25_dp*line%chord(j + 1)*line%t(:, j + 1)
      a1 = cross(p2 - p1, p3 - p2)/2._dp
      a2 = cross(p4 - p3, p1 - p4)/2._dp
      line%earea(j) = sqrt(dot_product(a1, a1)) + sqrt(dot_product(a2, a2))

      line%ec(j) = line%earea(j)/s_mag
      line%ettoc(j) = 0.5_dp*(line%thick(j + 1) + line%thick(j))
      line%epitch(j) = 0.5_dp*(line%pitch(j + 1) + line%pitch(j))
    end do
  end subroutine make_line_geometry

  subroutine rotate_line(line, origin, axis, theta)
    !! Rotate the stations of `line` by `theta` about `axis` through `origin`
    !! (chordwise vectors about `axis` through zero). Element geometry is not
    !! rebuilt; call make_line_geometry afterwards.
    implicit none
    type(actuator_line_t), intent(inout) :: line
    real(dp), intent(in) :: origin(3), axis(3), theta

    integer :: i

    do i = 1, line%n_elem + 1
      line%qc(:, i) = quat_rot(line%qc(:, i), theta, axis, origin)
      line%t(:, i) = quat_rot(line%t(:, i), theta, axis, [0._dp, 0._dp, 0._dp])
      line%t(:, i) = line%t(:, i)/norm2(line%t(:, i))
    end do
  end subroutine rotate_line

  subroutine pitch_line(line, pitch_angle)
    !! Pitch every station by `pitch_angle` [rad] about the line's root-to-tip
    !! direction, then rebuild the element geometry.
    implicit none
    type(actuator_line_t), intent(inout) :: line
    real(dp), intent(in) :: pitch_angle

    real(dp) :: s(3)
    integer :: i, n_st

    n_st = line%n_elem + 1
    s = line%qc(:, n_st) - line%cor
    s = s/norm2(s)
    do i = 1, n_st
      line%t(:, i) = quat_rot(line%t(:, i), pitch_angle, s, &
                              [0._dp, 0._dp, 0._dp])
      line%t(:, i) = line%t(:, i)/norm2(line%t(:, i))
      line%pitch(i) = line%pitch(i) - pitch_angle
    end do
    call make_line_geometry(line)
  end subroutine pitch_line

  subroutine populate_line_airfoils(line, airfoils)
    !! Give every element the airfoil whose thickness ratio is closest to its
    !! own, clamped to the thinnest and thickest available.
    implicit none
    type(actuator_line_t), intent(inout) :: line
    type(airfoil_t), intent(in) :: airfoils(:)

    real(dp) :: thicks(size(airfoils))
    integer :: ie, imin, imax, a

    do a = 1, size(airfoils)
      thicks(a) = airfoils(a)%tc
    end do
    imin = minloc(thicks, 1)
    imax = maxloc(thicks, 1)
    do ie = 1, line%n_elem
      if (line%ettoc(ie) >= thicks(imax)) then
        line%eairfoil(ie) = imax
      else if (line%ettoc(ie) <= thicks(imin)) then
        line%eairfoil(ie) = imin
      else
        line%eairfoil(ie) = minloc(abs(thicks - line%ettoc(ie)), 1)
      end if
    end do
  end subroutine populate_line_airfoils

  subroutine compute_line_forces(line, airfoils, rho_air, visc)
    !! Blade-element loads of every element from the sampled fluid velocity,
    !! the body velocity and the static polars.
    implicit none
    type(actuator_line_t), intent(inout) :: line
    type(airfoil_t), intent(in) :: airfoils(:)
    real(dp), intent(in) :: rho_air, visc

    real(dp) :: urel(3), urdn, urdc, ur, alpha, cl, cd, cm25, cn, ct
    real(dp) :: fn, ft, ms
    integer :: ie

    do ie = 1, line%n_elem
      urel = line%evel(:, ie) - line%evb(:, ie)
      urdn = dot_product(line%ne(:, ie), urel)
      urdc = dot_product(line%te(:, ie), urel)
      ur = sqrt(urdn**2 + urdc**2)
      alpha = atan2(urdn, urdc)

      line%eur(ie) = ur
      line%eaoa(ie) = alpha
      line%ere(ie) = ur*line%ec(ie)/visc

      call eval_static_coeffs(airfoils(line%eairfoil(ie)), line%ere(ie), &
                              alpha/deg2rad, cl, cd, cm25)

      cn = cl*cos(alpha) + cd*sin(alpha)
      ct = -cl*sin(alpha) + cd*cos(alpha)

      fn = 0.5_dp*rho_air*cn*line%earea(ie)*ur**2
      ft = 0.5_dp*rho_air*ct*line%earea(ie)*ur**2
      ms = 0.5_dp*rho_air*cm25*line%ec(ie)*line%earea(ie)*ur**2

      fn = fn*line%end_factor(ie)
      ft = ft*line%end_factor(ie)

      line%ecl(ie) = cl
      line%ecd(ie) = cd
      line%ecm(ie) = cm25
      line%ecn(ie) = cn
      line%ect(ie) = ct
      line%efn(ie) = fn
      line%eft(ie) = ft
      line%ems(ie) = ms
      line%ef(:, ie) = fn*line%ne(:, ie) + ft*line%te(:, ie)
    end do
  end subroutine compute_line_forces

  subroutine compute_tower_forces(tower, rho_air, visc, time, cl, cd, &
                                  strouhal, noise)
    !! Tower loads: constant drag and a lift oscillating at the shedding
    !! frequency St*U/D, perturbed by `noise` (uniform in [0, 1) per element).
    implicit none
    type(actuator_line_t), intent(inout) :: tower
    real(dp), intent(in) :: rho_air, visc, time, cl, cd, strouhal
    real(dp), intent(in) :: noise(:)

    real(dp) :: urdn, urdc, ur, alpha, freq, cn, ct, fn, ft
    integer :: ie

    do ie = 1, tower%n_elem
      urdn = dot_product(tower%ne(:, ie), tower%evel(:, ie))
      urdc = dot_product(tower%te(:, ie), tower%evel(:, ie))
      ur = sqrt(urdn**2 + urdc**2)
      alpha = atan2(urdn, urdc)
      tower%eur(ie) = ur
      tower%eaoa(ie) = alpha
      tower%ere(ie) = ur*tower%ec(ie)/visc

      freq = strouhal*ur/max(tower%ec(ie), 0.0001_dp)
      tower%ecl(ie) = cl*sin(2._dp*freq*pi*time) &
                      *(1._dp + 0.05_dp*(-1._dp + 2._dp*noise(ie)))
      tower%ecd(ie) = cd

      cn = tower%ecl(ie)*cos(alpha) + tower%ecd(ie)*sin(alpha)
      ct = -tower%ecl(ie)*sin(alpha) + tower%ecd(ie)*cos(alpha)
      fn = 0.5_dp*rho_air*cn*tower%earea(ie)*ur**2
      ft = 0.5_dp*rho_air*ct*tower%earea(ie)*ur**2

      tower%ecn(ie) = cn
      tower%ect(ie) = ct
      tower%efn(ie) = fn
      tower%eft(ie) = ft
      tower%ems(ie) = 0._dp
      tower%ef(:, ie) = fn*tower%ne(:, ie) + ft*tower%te(:, ie)
    end do
  end subroutine compute_tower_forces

  pure function quat_rot(v, theta, axis, origin) result(v_rot)
    !! Rotate point `v` by `theta` [rad] about `axis` through `origin`
    !! (right-hand rule). Equivalent to Incompact3d's quaternion QuatRot,
    !! written in Rodrigues' form; `axis` need not be normalised.
    implicit none
    real(dp), intent(in) :: v(3), theta, axis(3), origin(3)
    real(dp) :: v_rot(3)

    real(dp) :: n(3), p(3)

    n = axis/norm2(axis)
    p = v - origin
    v_rot = origin + p*cos(theta) + cross(n, p)*sin(theta) &
            + n*dot_product(n, p)*(1._dp - cos(theta))
  end function quat_rot

  pure function cross(a, b) result(c)
    implicit none
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end module m_alm_element
