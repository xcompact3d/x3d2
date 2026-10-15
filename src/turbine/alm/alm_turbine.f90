module m_alm_turbine
  !! Horizontal-axis wind turbines for the actuator line model (ALM).
  !!
  !! Ported from Incompact3d's `acl_turb.f90` and the turbine half of
  !! `acl_model.f90`. A turbine owns its blades and optional tower (each an
  !! actuator_line_t), its airfoil polars, and its rotor state (azimuth,
  !! angular velocity, collective pitch, controller). This module is free of
  !! fields and MPI: the flow velocity at the actuator points is sampled and
  !! the forces are projected by m_alm.
  !!
  !! A turbine is described by a `&TurbineSpecs` namelist file (`.turb`), in
  !! the same format as Incompact3d. Supported options:
  !!   - TypeFlag = 1 (horizontal axis)
  !!   - OperFlag = 1 (constant rotation, omega = uref*tsr/R),
  !!     2 (NREL 5 MW torque and pitch controller) or
  !!     4 (omega = U_upstream*tsr/R from an upstream probe)
  !!   - RotFlag = 1 (clockwise) or 2 (counterclockwise)
  !!   - towerFlag = 0/1, EndEffectsFlag = 0 (none), 1 (Glauert) or
  !!     2 (Shen et al. 2005) with TipCorr/RootCorr
  !! OperFlag = 3 (list controller), AddedMassFlag, DynStallFlag,
  !! RandomWalkForcingFlag, ConstantCirculationFlag and AeroElastFlag are not
  !! ported and stop the run when set.
  use m_common, only: dp, pi
  use m_logger, only: abort_run, log_info, log_warn, str
  use m_alm_airfoil, only: airfoil_t, read_airfoil
  use m_alm_controller, only: controller_t, init_controller, &
                              operate_controller, max_ctrl_blades
  use m_alm_element, only: actuator_line_t, allocate_line, &
                           make_line_geometry, rotate_line, pitch_line, &
                           populate_line_airfoils, read_line_geometry, &
                           compute_line_forces, compute_tower_forces, &
                           quat_rot, cross, deg2rad

  implicit none

  private
  public :: turbine_t, read_turbine_file, set_turbine_geometry, &
            compute_turbine_forces, advance_turbine, upstream_probe

  integer, parameter, public :: OPER_CONSTANT = 1, OPER_NREL = 2, &
                                OPER_UPSTREAM = 4

  integer, parameter :: max_n_airfoils = 80

  type :: turbine_t
    character(len=100) :: name = ''
    character(len=256) :: blade_geom_file = ''
    character(len=256) :: tower_geom_file = ''
    integer :: n_blades = 0
    real(dp) :: origin(3) = 0._dp
    real(dp) :: rot_n(3) = [1._dp, 0._dp, 0._dp] !! rotation axis
    real(dp) :: shaft_tilt_angle = 0._dp     !! [deg]
    real(dp) :: blade_cone_angle = 0._dp     !! [deg]
    real(dp) :: yaw_angle = 0._dp            !! [deg]
    real(dp) :: rmax = 0._dp                 !! rotor radius
    real(dp) :: area = 0._dp                 !! swept area
    real(dp) :: blade_inertia = 0._dp
    real(dp) :: i_rotor = 0._dp
    integer :: oper = OPER_CONSTANT
    logical :: clockwise = .true.
    real(dp) :: tsr = 0._dp
    real(dp) :: uref = 0._dp

    ! Rotor state
    real(dp) :: azim_angle = 0._dp
    real(dp) :: angular_vel = 0._dp
    real(dp) :: delta_omega = 0._dp
    real(dp) :: cbp = 0._dp, cbp_old = 0._dp !! collective blade pitch [rad]
    real(dp) :: u_upstream(3) = 0._dp
    type(controller_t) :: controller

    ! End effects
    logical :: end_effects = .false.
    logical :: end_glauert = .false.
    logical :: end_shen = .false.
    logical :: tip_correction = .false.
    logical :: root_correction = .false.
    real(dp) :: shen_c1 = 0.125_dp, shen_c2 = 21._dp

    ! Tower
    logical :: has_tower = .false.
    real(dp) :: tower_height = 0._dp
    real(dp) :: tower_offset = 0._dp
    real(dp) :: tower_drag = 1._dp
    real(dp) :: tower_lift = 0.3_dp
    real(dp) :: tower_strouhal = 0.21_dp

    type(airfoil_t), allocatable :: airfoils(:)
    type(actuator_line_t), allocatable :: blade(:)
    type(actuator_line_t) :: tower

    ! Performance
    real(dp) :: thrust = 0._dp, torque = 0._dp, power = 0._dp
    real(dp) :: ct = 0._dp, cp = 0._dp, ctr = 0._dp
    real(dp) :: cfx = 0._dp, cfy = 0._dp, cfz = 0._dp
  end type turbine_t

contains

  subroutine read_turbine_file(turb, filename, u_inflow)
    !! Read a `&TurbineSpecs` namelist file and the airfoils it lists.
    !! `u_inflow` is the initial reference speed of the NREL controller
    !! (Incompact3d uses 0.5*(u1 + u2)).
    implicit none
    type(turbine_t), intent(inout) :: turb
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: u_inflow

    character(len=100) :: name
    character(len=256) :: blade_geom, tower_geom, dynstall_param_file, &
                          AeroElastInputFile, AeroElastSolverFile, &
                          list_controller_file
    character(len=256) :: afname(max_n_airfoils)
    real(dp) :: origin(3), toweroffset, tower_drag, tower_lift, &
                tower_strouhal, uref, tsr, ShenC1, ShenC2, GammaCirc, &
                BladeInertia, GeneratorInertia, GBRatio, GBEfficiency, &
                RatedGenSpeed, RatedLimitGenTorque, CutInGenSpeed, &
                Region2StartGenSpeed, Region2EndGenSpeed, Kgen, RatedPower, &
                MaximumTorque, yaw_angle, shaft_tilt_angle, blade_cone_angle
    integer :: numblades, numfoil, towerFlag, TypeFlag, OperFlag, RotFlag, &
               AddedMassFlag, DynStallFlag, EndEffectsFlag, TipCorr, &
               RootCorr, RandomWalkForcingFlag, AeroElastFlag, &
               AeroElastModel, ConstantCirculationFlag
    integer :: unit, ios, k

    namelist /TurbineSpecs/ name, origin, numblades, blade_geom, numfoil, &
      afname, towerFlag, towerOffset, tower_geom, tower_drag, tower_lift, &
      tower_strouhal, TypeFlag, OperFlag, tsr, uref, RotFlag, &
      AddedMassFlag, RandomWalkForcingFlag, DynStallFlag, &
      dynstall_param_file, EndEffectsFlag, TipCorr, RootCorr, ShenC1, &
      ShenC2, ConstantCirculationFlag, GammaCirc, yaw_angle, &
      shaft_tilt_angle, blade_cone_angle, AeroElastFlag, AeroElastModel, &
      AeroElastInputFile, AeroElastSolverFile, BladeInertia, &
      GeneratorInertia, GBRatio, GBEfficiency, RatedGenSpeed, &
      RatedLimitGenTorque, CutInGenSpeed, Region2StartGenSpeed, &
      Region2EndGenSpeed, Kgen, RatedPower, MaximumTorque, &
      list_controller_file

    ! Defaults (as Incompact3d)
    name = 'turbine'
    origin = 0._dp
    numblades = 3
    blade_geom = ''
    numfoil = 0
    afname = ''
    towerFlag = 0
    toweroffset = 0._dp
    tower_geom = ''
    tower_lift = 0.3_dp
    tower_drag = 1._dp
    tower_strouhal = 0.21_dp
    TypeFlag = 1
    OperFlag = 1
    tsr = 0._dp
    uref = 0._dp
    RotFlag = 1
    AddedMassFlag = 0
    RandomWalkForcingFlag = 0
    DynStallFlag = 0
    EndEffectsFlag = 0
    ConstantCirculationFlag = 0
    TipCorr = 0
    RootCorr = 0
    ShenC1 = 0.125_dp
    ShenC2 = 21._dp
    GammaCirc = 0._dp
    yaw_angle = 0._dp
    shaft_tilt_angle = 0._dp
    blade_cone_angle = 0._dp
    AeroElastFlag = 0
    AeroElastModel = 0
    BladeInertia = 0._dp
    GeneratorInertia = 0._dp
    GBRatio = 1._dp
    GBEfficiency = 1._dp
    RatedGenSpeed = 0._dp
    RatedLimitGenTorque = 0._dp
    CutInGenSpeed = 0._dp
    Region2StartGenSpeed = 0._dp
    Region2EndGenSpeed = 0._dp
    Kgen = 0._dp
    RatedPower = 0._dp
    MaximumTorque = 0._dp

    open (newunit=unit, file=trim(filename), status='old', action='read', &
          iostat=ios)
    if (ios /= 0) call abort_run('ALM: cannot open turbine file: ' &
                                 //trim(filename), collective=.true.)
    read (unit, nml=TurbineSpecs, iostat=ios)
    close (unit)
    if (ios /= 0) call abort_run('ALM: error reading &TurbineSpecs in: ' &
                                 //trim(filename), collective=.true.)

    turb%name = name
    turb%origin = origin
    turb%n_blades = numblades
    turb%blade_geom_file = blade_geom
    turb%blade_inertia = BladeInertia
    if (numblades < 1) call abort_run('ALM: numblades must be positive in ' &
                                      //trim(filename), collective=.true.)

    if (numfoil < 1 .or. numfoil > max_n_airfoils) &
      call abort_run('ALM: provide between 1 and '//str(max_n_airfoils) &
                     //' airfoils (numfoil) in '//trim(filename), &
                     collective=.true.)
    allocate (turb%airfoils(numfoil))
    do k = 1, numfoil
      call read_airfoil(turb%airfoils(k), afname(k))
      call log_info('ALM: '//trim(name)//' airfoil '//str(k)//' is ' &
                    //trim(afname(k)))
    end do

    if (towerFlag == 1) then
      turb%has_tower = .true.
      turb%tower_offset = toweroffset
      turb%tower_geom_file = tower_geom
      turb%tower_drag = tower_drag
      turb%tower_lift = tower_lift
      turb%tower_strouhal = tower_strouhal
    end if

    if (TypeFlag /= 1) call abort_run('ALM: only horizontal-axis turbines &
                                      &(TypeFlag = 1) are supported.', &
                                      collective=.true.)
    turb%shaft_tilt_angle = shaft_tilt_angle
    turb%yaw_angle = yaw_angle
    turb%blade_cone_angle = blade_cone_angle

    select case (OperFlag)
    case (OPER_CONSTANT, OPER_UPSTREAM)
      turb%uref = uref
      turb%tsr = tsr
    case (OPER_NREL)
      turb%uref = u_inflow
      turb%tsr = tsr
      if (numblades > max_ctrl_blades) &
        call abort_run('ALM: the NREL controller supports at most ' &
                       //str(max_ctrl_blades)//' blades.', collective=.true.)
      call init_controller(turb%controller, GeneratorInertia, GBRatio, &
                           RatedGenSpeed, RatedLimitGenTorque, &
                           CutInGenSpeed, Region2StartGenSpeed, Kgen, &
                           RatedPower, MaximumTorque)
    case (3)
      call abort_run('ALM: the list controller (OperFlag = 3) is not &
                     &supported yet.', collective=.true.)
    case default
      call abort_run('ALM: unknown OperFlag in '//trim(filename), &
                     collective=.true.)
    end select
    turb%oper = OperFlag

    select case (RotFlag)
    case (1)
      turb%clockwise = .true.
    case (2)
      turb%clockwise = .false.
    case default
      call abort_run('ALM: RotFlag must be 1 (clockwise) or 2 &
                     &(counterclockwise).', collective=.true.)
    end select

    if (AddedMassFlag /= 0 .or. DynStallFlag /= 0 &
        .or. RandomWalkForcingFlag /= 0 .or. ConstantCirculationFlag /= 0 &
        .or. AeroElastFlag /= 0) then
      call abort_run('ALM: added mass, dynamic stall, random-walk forcing, &
                     &constant circulation and aeroelasticity are not &
                     &supported yet; set their flags to 0.', &
                     collective=.true.)
    end if

    select case (EndEffectsFlag)
    case (0)
    case (1, 2)
      turb%end_effects = .true.
      turb%end_glauert = EndEffectsFlag == 1
      turb%end_shen = EndEffectsFlag == 2
      turb%shen_c1 = ShenC1
      turb%shen_c2 = ShenC2
      turb%tip_correction = TipCorr == 1
      turb%root_correction = RootCorr == 1
    case default
      call abort_run('ALM: EndEffectsFlag must be 0, 1 or 2.', &
                     collective=.true.)
    end select
  end subroutine read_turbine_file

  subroutine set_turbine_geometry(turb, restart)
    !! Build the blades (and tower) at the current azimuth and collective
    !! pitch. Unless `restart`, the angular velocity is set to uref*tsr/R.
    implicit none
    type(turbine_t), intent(inout) :: turb
    logical, intent(in) :: restart

    real(dp), allocatable :: r_over_l(:), c_over_l(:), pitch(:), thick(:)
    real(dp) :: span_vec(3), theta, length, yaw_axis(3), tilt_axis(3)
    integer :: n_st, ib, ist

    yaw_axis = [0._dp, 1._dp, 0._dp]
    tilt_axis = [0._dp, 0._dp, 1._dp]

    ! Unyawed, untilted rotation axis; set again so a rebuild on restart
    ! does not yaw and tilt it twice
    turb%rot_n = [1._dp, 0._dp, 0._dp]
    if (.not. turb%clockwise) turb%rot_n = -turb%rot_n

    call read_line_geometry(turb%blade_geom_file, turb%rmax, span_vec, &
                            r_over_l, c_over_l, pitch, thick)
    n_st = size(r_over_l)

    ! Blades are built along z through the origin, coned, translated to the
    ! turbine origin and spread in azimuth about the rotation axis
    if (allocated(turb%blade)) deallocate (turb%blade)
    allocate (turb%blade(turb%n_blades))
    theta = 2._dp*pi/turb%n_blades
    do ib = 1, turb%n_blades
      associate (bl => turb%blade(ib))
        call allocate_line(bl, n_st)
        write (bl%name, '(A,"_blade_",I0)') trim(turb%name), ib
        bl%cor = turb%origin
        bl%length = turb%rmax
        bl%flip_n = turb%clockwise
        do ist = 1, n_st
          bl%qc(:, ist) = [0._dp, 0._dp, r_over_l(ist)*turb%rmax]
          if (turb%clockwise) then
            bl%t(:, ist) = [sin(pitch(ist)*deg2rad), &
                            cos(pitch(ist)*deg2rad), 0._dp]
          else
            bl%t(:, ist) = [sin(pitch(ist)*deg2rad), &
                            -cos(pitch(ist)*deg2rad), 0._dp]
          end if
          bl%chord(ist) = c_over_l(ist)*turb%rmax
          bl%thick(ist) = thick(ist)
          bl%pitch(ist) = pitch(ist)*deg2rad

          ! Blade cone angle (about y)
          bl%qc(:, ist) = quat_rot(bl%qc(:, ist), &
                                   turb%blade_cone_angle*deg2rad, &
                                   yaw_axis, [0._dp, 0._dp, 0._dp])
          bl%t(:, ist) = quat_rot(bl%t(:, ist), &
                                  turb%blade_cone_angle*deg2rad, &
                                  yaw_axis, [0._dp, 0._dp, 0._dp])
          bl%qc(:, ist) = bl%qc(:, ist) + bl%cor
        end do

        call rotate_line(bl, bl%cor, turb%rot_n, &
                         (ib - 1)*theta + turb%azim_angle)
        ! pitch_line also builds the element geometry
        call pitch_line(bl, turb%cbp_old)
        call populate_line_airfoils(bl, turb%airfoils)
      end associate
    end do

    ! Yaw and tilt the rotor, and its rotation axis
    call rotate_turbine(turb, yaw_axis, turb%yaw_angle*deg2rad)
    call rotate_turbine(turb, tilt_axis, -turb%shaft_tilt_angle*deg2rad)
    turb%rot_n = quat_rot(turb%rot_n, turb%yaw_angle*deg2rad, yaw_axis, &
                          [0._dp, 0._dp, 0._dp])
    turb%rot_n = quat_rot(turb%rot_n, -turb%shaft_tilt_angle*deg2rad, &
                          tilt_axis, [0._dp, 0._dp, 0._dp])

    if (turb%has_tower) then
      call read_line_geometry(turb%tower_geom_file, length, span_vec, &
                              r_over_l, c_over_l, pitch, thick)
      turb%tower_height = length
      n_st = size(r_over_l)
      associate (tw => turb%tower)
        call allocate_line(tw, n_st)
        tw%name = trim(turb%name)//'_tower'
        tw%cor = turb%origin
        tw%length = length
        tw%flip_n = .false.
        ! The tower stands on y = 0, offset in x from the rotor centre
        do ist = 1, n_st
          tw%qc(:, ist) = [tw%cor(1) + turb%tower_offset, &
                           r_over_l(ist)*length, tw%cor(3)]
          tw%t(:, ist) = [1._dp, 0._dp, 0._dp]
          tw%chord(ist) = c_over_l(ist)*length
          tw%thick(ist) = thick(ist)
          tw%pitch(ist) = pitch(ist)*deg2rad
        end do
        call make_line_geometry(tw)
        tw%evb = 0._dp
      end associate
    end if

    if (.not. restart) turb%angular_vel = turb%uref*turb%tsr/turb%rmax
    turb%area = pi*turb%rmax**2
    turb%i_rotor = turb%n_blades*turb%blade_inertia
    call compute_rot_vel(turb)

    call log_info('ALM: turbine '//trim(turb%name)//': '// &
                  str(turb%n_blades)//' blades, R = '// &
                  str(turb%rmax, '(es12.5)')//', origin = ('// &
                  str(turb%origin(1), '(es12.5)')//', '// &
                  str(turb%origin(2), '(es12.5)')//', '// &
                  str(turb%origin(3), '(es12.5)')//'), axis = ('// &
                  str(turb%rot_n(1), '(f8.5)')//', '// &
                  str(turb%rot_n(2), '(f8.5)')//', '// &
                  str(turb%rot_n(3), '(f8.5)')//')')
  end subroutine set_turbine_geometry

  subroutine compute_turbine_forces(turb, rho_air, visc, time, tower_noise)
    !! Loads on the blades (with end effects) and tower from the element
    !! velocities, then the rotor performance.
    implicit none
    type(turbine_t), intent(inout) :: turb
    real(dp), intent(in) :: rho_air, visc, time
    real(dp), intent(in) :: tower_noise(:)

    integer :: ib

    if (turb%end_effects) call compute_end_effects(turb)
    do ib = 1, turb%n_blades
      call compute_line_forces(turb%blade(ib), turb%airfoils, rho_air, visc)
    end do
    call compute_performance(turb, rho_air)
    if (turb%has_tower) then
      call compute_tower_forces(turb%tower, rho_air, visc, time, &
                                turb%tower_lift, turb%tower_drag, &
                                turb%tower_strouhal, tower_noise)
    end if
  end subroutine compute_turbine_forces

  subroutine advance_turbine(turb, time, dt)
    !! Advance the rotor by one time step: update the angular velocity (and
    !! pitch) from the operating mode, then rotate the blades. For the
    !! controlled modes `u_upstream` must already hold the probe velocity.
    implicit none
    type(turbine_t), intent(inout) :: turb
    real(dp), intent(in) :: time, dt

    real(dp) :: theta, delta_pitch
    integer :: ib

    select case (turb%oper)
    case (OPER_CONSTANT)
      theta = turb%angular_vel*dt
      turb%azim_angle = turb%azim_angle + theta
      call rotate_turbine(turb, turb%rot_n, theta)
      call compute_rot_vel(turb)

    case (OPER_NREL)
      turb%uref = norm2(turb%u_upstream)
      call operate_controller(turb%controller, time, turb%n_blades, &
                              turb%angular_vel)
      turb%delta_omega = (turb%torque - turb%controller%gear_box_ratio &
                          *turb%controller%gen_trq) &
                         /(turb%i_rotor + turb%controller%gear_box_ratio**2 &
                           *turb%controller%i_generator)*dt
      turb%angular_vel = turb%angular_vel + turb%delta_omega
      theta = turb%angular_vel*dt
      turb%azim_angle = turb%azim_angle + theta

      ! Collective pitch (the controller acts on clockwise rotors only)
      if (.not. turb%clockwise) &
        call abort_run('ALM: the NREL controller needs a clockwise rotor.')
      do ib = 1, turb%n_blades
        turb%cbp = -turb%controller%pit_com(ib)
        delta_pitch = turb%cbp - turb%cbp_old
        call pitch_line(turb%blade(ib), delta_pitch)
      end do
      turb%cbp_old = turb%cbp

      call rotate_turbine(turb, turb%rot_n, theta)
      call compute_rot_vel(turb)
      turb%controller%status = turb%controller%status + 1

    case (OPER_UPSTREAM)
      turb%uref = norm2(turb%u_upstream)
      turb%angular_vel = turb%uref*turb%tsr/turb%rmax
      theta = turb%angular_vel*dt
      turb%azim_angle = turb%azim_angle + theta
      call rotate_turbine(turb, turb%rot_n, theta)
      call compute_rot_vel(turb)
    end select
  end subroutine advance_turbine

  pure function upstream_probe(turb) result(point)
    !! Probe point for the controlled modes, two radii along -rot_n from the
    !! rotor centre (as Incompact3d; for a counterclockwise rotor rot_n
    !! points upstream, so the probe lies downstream).
    implicit none
    type(turbine_t), intent(in) :: turb
    real(dp) :: point(3)

    point = turb%origin - turb%rot_n*2._dp*turb%rmax
  end function upstream_probe

  subroutine compute_performance(turb, rho_air)
    implicit none
    type(turbine_t), intent(inout) :: turb
    real(dp), intent(in) :: rho_air

    real(dp) :: f_tot(3), torq_tot, tr(3), q_ref
    integer :: ib, ie

    f_tot = 0._dp
    torq_tot = 0._dp
    do ib = 1, turb%n_blades
      associate (bl => turb%blade(ib))
        do ie = 1, bl%n_elem
          tr = cross(bl%pe(:, ie) - bl%cor, bl%ef(:, ie))
          torq_tot = torq_tot + dot_product(tr, turb%rot_n) &
                     + bl%ems(ie)*dot_product(bl%se(:, ie), turb%rot_n)
          f_tot = f_tot + bl%ef(:, ie)
        end do
      end associate
    end do

    q_ref = 0.5_dp*rho_air*turb%area*turb%uref**2
    turb%cfx = f_tot(1)/q_ref
    turb%cfy = f_tot(2)/q_ref
    turb%cfz = f_tot(3)/q_ref
    turb%thrust = norm2(f_tot)
    turb%ct = sqrt(turb%cfx**2 + turb%cfy**2 + turb%cfz**2)
    turb%torque = torq_tot
    turb%ctr = torq_tot/(q_ref*turb%rmax)
    turb%power = abs(torq_tot)*turb%angular_vel
    turb%cp = abs(turb%ctr)*turb%tsr
  end subroutine compute_performance

  subroutine compute_end_effects(turb)
    !! Glauert or Shen et al. (2005) tip and root loss factors.
    implicit none
    type(turbine_t), intent(inout) :: turb

    real(dp) :: g1, phi, rroot, rtip, f_tip, f_root, e
    integer :: ib, ie

    if (turb%end_glauert) then
      g1 = 1._dp
    else
      g1 = exp(-turb%shen_c1*(turb%n_blades*turb%tsr - turb%shen_c2)) + 0.1_dp
    end if

    do ib = 1, turb%n_blades
      associate (bl => turb%blade(ib))
        do ie = 1, bl%n_elem
          phi = bl%eaoa(ie) + bl%epitch(ie)
          rroot = bl%erdist(ie)/turb%rmax
          rtip = (turb%rmax - bl%erdist(ie))/turb%rmax
          f_tip = 1._dp
          f_root = 1._dp
          ! As Incompact3d, the "tip" factor uses r/R and the "root" factor
          ! (R - r)/R
          if (turb%tip_correction) then
            e = exp(-g1*turb%n_blades*0.5_dp*(1._dp/rroot - 1._dp)/sin(phi))
            if (abs(e) > 1._dp) then
              call log_warn('ALM: tip correction failed, phi = ' &
                            //str(phi, '(es12.5)'))
            else
              f_tip = 2._dp/pi*acos(e)
            end if
          end if
          if (turb%root_correction) then
            e = exp(-g1*turb%n_blades*0.5_dp*(1._dp/rtip - 1._dp)/sin(phi))
            if (abs(e) > 1._dp) then
              call log_warn('ALM: root correction failed, phi = ' &
                            //str(phi, '(es12.5)'))
            else
              f_root = 2._dp/pi*acos(e)
            end if
          end if
          bl%end_factor(ie) = f_tip*f_root
        end do
      end associate
    end do
  end subroutine compute_end_effects

  subroutine compute_rot_vel(turb)
    !! Rigid-body velocity omega x r of every blade element.
    implicit none
    type(turbine_t), intent(inout) :: turb

    real(dp) :: w_rot(3)
    integer :: ib, ie

    w_rot = turb%angular_vel*turb%rot_n
    do ib = 1, turb%n_blades
      associate (bl => turb%blade(ib))
        do ie = 1, bl%n_elem
          bl%evb(:, ie) = cross(w_rot, bl%pe(:, ie) - bl%cor)
        end do
      end associate
    end do
  end subroutine compute_rot_vel

  subroutine rotate_turbine(turb, axis, theta)
    !! Rotate all blades by `theta` about `axis` through the turbine origin.
    implicit none
    type(turbine_t), intent(inout) :: turb
    real(dp), intent(in) :: axis(3), theta

    integer :: ib

    do ib = 1, turb%n_blades
      call rotate_line(turb%blade(ib), turb%origin, axis, theta)
      call make_line_geometry(turb%blade(ib))
    end do
  end subroutine rotate_turbine

end module m_alm_turbine
