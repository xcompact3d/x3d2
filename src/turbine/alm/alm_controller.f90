module m_alm_controller
  !! Variable-speed generator-torque controller and PI collective blade-pitch
  !! controller of the NREL offshore 5 MW baseline turbine, as written by
  !! J. Jonkman (NREL/NWTC) for the IEA Annex XXIII OC3 studies (DISCON).
  !!
  !! Ported from Incompact3d's `acl_controller.f90`. The torque law has five
  !! regions (1, 1 1/2, 2, 2 1/2 and 3) on the low-pass filtered generator
  !! speed; the pitch controller is a gain-scheduled PI on the speed error,
  !! with pitch-angle and pitch-rate limits. The pitch gains and limits are
  !! those of the NREL 5 MW turbine; the torque parameters come from the
  !! turbine input file.
  use m_common, only: dp
  use m_logger, only: abort_run, log_info

  implicit none

  private
  public :: controller_t, init_controller, operate_controller

  integer, parameter, public :: max_ctrl_blades = 3

  type :: controller_t
    ! Input parameters
    real(dp) :: corner_freq = 0._dp !! speed filter corner frequency [rad/s]
    real(dp) :: pc_dt = 0._dp       !! pitch controller interval [s]
    real(dp) :: pc_ki = 0._dp       !! integral gain at rated pitch [-]
    real(dp) :: pc_kk = 0._dp       !! pitch where dP/dpitch doubles [rad]
    real(dp) :: pc_kp = 0._dp       !! proportional gain at rated pitch [s]
    real(dp) :: pc_max_pit = 0._dp  !! maximum pitch [rad]
    real(dp) :: pc_max_rat = 0._dp  !! maximum pitch rate [rad/s]
    real(dp) :: pc_min_pit = 0._dp  !! minimum pitch [rad]
    real(dp) :: pc_ref_spd = 0._dp  !! reference generator speed [rad/s]
    real(dp) :: vs_ct_in_sp = 0._dp !! cut-in generator speed [rad/s]
    real(dp) :: vs_dt = 0._dp       !! torque controller interval [s]
    real(dp) :: vs_max_rat = 0._dp  !! maximum torque rate [N m/s]
    real(dp) :: vs_max_tq = 0._dp   !! maximum generator torque [N m]
    real(dp) :: vs_rgn2_k = 0._dp   !! region 2 torque constant
    real(dp) :: vs_rgn2_sp = 0._dp  !! region 1 1/2 to 2 speed [rad/s]
    real(dp) :: vs_rgn3_mp = 0._dp  !! pitch above which region 3 applies
    real(dp) :: vs_rt_gn_sp = 0._dp !! rated generator speed [rad/s]
    real(dp) :: vs_rt_pwr = 0._dp   !! rated power [W]
    real(dp) :: vs_sl_pc = 0._dp    !! region 2 1/2 slip [%]
    real(dp) :: gear_box_ratio = 0._dp
    real(dp) :: i_generator = 0._dp !! generator moment of inertia

    ! State
    real(dp) :: bl_pitch(max_ctrl_blades) = 0._dp
    real(dp) :: gen_speed = 0._dp
    real(dp) :: gen_speed_f = 0._dp
    real(dp) :: gen_trq = 0._dp
    real(dp) :: gk = 0._dp
    real(dp) :: int_spd_err = 0._dp
    real(dp) :: last_gen_trq = 0._dp
    real(dp) :: last_time = 0._dp
    real(dp) :: last_time_pc = 0._dp
    real(dp) :: last_time_vs = 0._dp
    real(dp) :: pit_com(max_ctrl_blades) = 0._dp
    real(dp) :: vs_slope15 = 0._dp
    real(dp) :: vs_slope25 = 0._dp
    real(dp) :: vs_sy_sp = 0._dp
    real(dp) :: vs_tr_gn_sp = 0._dp
    integer :: status = 0 !! 0 on the first call, then > 0
  end type controller_t

contains

  subroutine init_controller(ctrl, generator_inertia, gb_ratio, &
                             rated_gen_speed, rated_limit_gen_torque, &
                             cut_in_gen_speed, region2_start_gen_speed, &
                             kgen, rated_power, maximum_torque)
    implicit none
    type(controller_t), intent(inout) :: ctrl
    real(dp), intent(in) :: generator_inertia, gb_ratio, rated_gen_speed, &
                            rated_limit_gen_torque, cut_in_gen_speed, &
                            region2_start_gen_speed, kgen, rated_power, &
                            maximum_torque

    ctrl%gear_box_ratio = gb_ratio
    ctrl%i_generator = generator_inertia
    ctrl%vs_ct_in_sp = cut_in_gen_speed
    ctrl%vs_rt_gn_sp = rated_gen_speed
    ctrl%vs_max_rat = rated_limit_gen_torque
    ctrl%vs_rgn2_k = kgen
    ctrl%vs_rgn2_sp = region2_start_gen_speed
    ctrl%vs_rt_pwr = rated_power
    ctrl%vs_max_tq = maximum_torque

    ctrl%corner_freq = 1.570796_dp
    ctrl%pc_dt = 0.00125_dp
    ctrl%pc_ki = 0.008068634_dp
    ctrl%pc_kk = 0.1099965_dp
    ctrl%pc_kp = 0.01882681_dp
    ctrl%pc_max_pit = 1.570796_dp
    ctrl%pc_max_rat = 0.1396263_dp
    ctrl%pc_min_pit = 0._dp
    ctrl%pc_ref_spd = 122.9096_dp
    ctrl%vs_dt = 0.00125_dp
    ctrl%vs_rgn3_mp = 0.01745329_dp
    ctrl%vs_sl_pc = 10._dp

    ! Derived torque-control parameters
    ctrl%vs_sy_sp = ctrl%vs_rt_gn_sp/(1._dp + 0.01_dp*ctrl%vs_sl_pc)
    ctrl%vs_slope15 = (ctrl%vs_rgn2_k*ctrl%vs_rgn2_sp*ctrl%vs_rgn2_sp) &
                      /(ctrl%vs_rgn2_sp - ctrl%vs_ct_in_sp)
    ctrl%vs_slope25 = (ctrl%vs_rt_pwr/ctrl%vs_rt_gn_sp) &
                      /(ctrl%vs_rt_gn_sp - ctrl%vs_sy_sp)
    if (ctrl%vs_rgn2_k == 0._dp) then
      ! Flat region 2 torque
      ctrl%vs_tr_gn_sp = ctrl%vs_sy_sp
    else
      ! Region 2 torque quadratic in speed
      ctrl%vs_tr_gn_sp = (ctrl%vs_slope25 &
                          - sqrt(ctrl%vs_slope25*(ctrl%vs_slope25 &
                                                  - 4._dp*ctrl%vs_rgn2_k &
                                                  *ctrl%vs_sy_sp))) &
                         /(2._dp*ctrl%vs_rgn2_k)
    end if

    call check(ctrl%vs_ct_in_sp >= 0._dp, 'VS_CtInSp must not be negative.')
    call check(ctrl%vs_rgn2_sp > ctrl%vs_ct_in_sp, &
               'VS_Rgn2Sp must be greater than VS_CtInSp.')
    call check(ctrl%vs_tr_gn_sp >= ctrl%vs_rgn2_sp, &
               'VS_TrGnSp must not be less than VS_Rgn2Sp.')
    call check(ctrl%vs_max_rat > 0._dp, 'VS_MaxRat must be greater than zero.')
    call check(ctrl%vs_rt_pwr >= 0._dp, 'VS_RtPwr must not be negative.')
    call check(ctrl%vs_rgn2_k >= 0._dp, 'VS_Rgn2K must not be negative.')
    call check(ctrl%vs_rgn2_k*ctrl%vs_rt_gn_sp*ctrl%vs_rt_gn_sp &
               <= ctrl%vs_rt_pwr/ctrl%vs_rt_gn_sp, &
               'VS_Rgn2K*VS_RtGnSp^2 must not be greater than &
               &VS_RtPwr/VS_RtGnSp.')
    call check(ctrl%vs_max_tq >= ctrl%vs_rt_pwr/ctrl%vs_rt_gn_sp, &
               'VS_RtPwr/VS_RtGnSp must not be greater than VS_MaxTq.')

    ctrl%status = 0
    call log_info('ALM: NREL 5MW torque and pitch controller enabled')

  contains

    subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) call abort_run('ALM controller: '//message, &
                                   collective=.true.)
    end subroutine check
  end subroutine init_controller

  subroutine operate_controller(ctrl, time, n_blades, rot_speed)
    !! Advance the controller to `time` given the rotor speed [rad/s]. The
    !! commanded generator torque is left in gen_trq and the commanded blade
    !! pitch angles in pit_com.
    implicit none
    type(controller_t), intent(inout) :: ctrl
    real(dp), intent(in) :: time, rot_speed
    integer, intent(in) :: n_blades

    ! Slightly above one, so a controller interval equal to the time step
    ! fires every step despite rounding
    real(dp), parameter :: one_plus_eps = 1._dp + epsilon(1._dp)
    real(dp) :: alpha, elap_time, trq_rate, spd_err, pit_com_p, pit_com_i, &
                pit_com_t, pit_rate
    integer :: k

    ctrl%gen_speed = rot_speed*ctrl%gear_box_ratio

    if (ctrl%status == 0) then
      ! First call: start the filter and integrator so the pitch stays
      ! unchanged for a zero speed error, and fire both controllers
      ctrl%gen_speed_f = ctrl%gen_speed
      ctrl%pit_com = ctrl%bl_pitch
      ctrl%gk = 1._dp/(1._dp + ctrl%pit_com(1)/ctrl%pc_kk)
      ctrl%int_spd_err = ctrl%pit_com(1)/(ctrl%gk*ctrl%pc_ki)
      ctrl%last_time = time
      ctrl%last_time_pc = time - ctrl%pc_dt
      ctrl%last_time_vs = time - ctrl%vs_dt
      return
    end if

    ! Single-pole low-pass filter on the generator speed
    alpha = exp((ctrl%last_time - time)*ctrl%corner_freq)
    ctrl%gen_speed_f = (1._dp - alpha)*ctrl%gen_speed &
                       + alpha*ctrl%gen_speed_f

    ! Variable-speed torque control
    elap_time = time - ctrl%last_time_vs
    if (time*one_plus_eps - ctrl%last_time_vs >= ctrl%vs_dt) then
      if (ctrl%gen_speed_f >= ctrl%vs_rt_gn_sp &
          .or. ctrl%pit_com(1) >= ctrl%vs_rgn3_mp) then
        ! Region 3: constant power
        ctrl%gen_trq = ctrl%vs_rt_pwr/ctrl%gen_speed_f
      else if (ctrl%gen_speed_f <= ctrl%vs_ct_in_sp) then
        ! Region 1: no torque
        ctrl%gen_trq = 0._dp
      else if (ctrl%gen_speed_f < ctrl%vs_rgn2_sp) then
        ! Region 1 1/2: linear ramp from zero to optimal
        ctrl%gen_trq = ctrl%vs_slope15*(ctrl%gen_speed_f - ctrl%vs_ct_in_sp)
      else if (ctrl%gen_speed_f < ctrl%vs_tr_gn_sp) then
        ! Region 2: optimal torque, quadratic in speed
        ctrl%gen_trq = ctrl%vs_rgn2_k*ctrl%gen_speed_f*ctrl%gen_speed_f
      else
        ! Region 2 1/2: induction-generator transition
        ctrl%gen_trq = ctrl%vs_slope25*(ctrl%gen_speed_f - ctrl%vs_sy_sp)
      end if

      ! Saturate with the torque and torque-rate limits
      ctrl%gen_trq = min(ctrl%gen_trq, ctrl%vs_max_tq)
      trq_rate = (ctrl%gen_trq - ctrl%last_gen_trq)/elap_time
      trq_rate = min(max(trq_rate, -ctrl%vs_max_rat), ctrl%vs_max_rat)
      ctrl%gen_trq = ctrl%last_gen_trq + trq_rate*elap_time
      ctrl%last_time_vs = time
      ctrl%last_gen_trq = ctrl%gen_trq
    end if

    ! Pitch control
    elap_time = time - ctrl%last_time_pc
    if (time*one_plus_eps - ctrl%last_time_pc >= ctrl%pc_dt) then
      ! Gain scheduling on the previously commanded pitch of blade 1
      ctrl%gk = 1._dp/(1._dp + ctrl%pit_com(1)/ctrl%pc_kk)

      ! Speed error and its integral, saturated with the pitch limits
      spd_err = ctrl%gen_speed_f - ctrl%pc_ref_spd
      ctrl%int_spd_err = ctrl%int_spd_err + spd_err*elap_time
      ctrl%int_spd_err = min(max(ctrl%int_spd_err, &
                                 ctrl%pc_min_pit/(ctrl%gk*ctrl%pc_ki)), &
                             ctrl%pc_max_pit/(ctrl%gk*ctrl%pc_ki))

      pit_com_p = ctrl%gk*ctrl%pc_kp*spd_err
      pit_com_i = ctrl%gk*ctrl%pc_ki*ctrl%int_spd_err
      pit_com_t = min(max(pit_com_p + pit_com_i, ctrl%pc_min_pit), &
                      ctrl%pc_max_pit)

      ! Saturate each blade's command with the pitch-rate limit
      do k = 1, min(n_blades, max_ctrl_blades)
        pit_rate = (pit_com_t - ctrl%bl_pitch(k))/elap_time
        pit_rate = min(max(pit_rate, -ctrl%pc_max_rat), ctrl%pc_max_rat)
        ctrl%pit_com(k) = ctrl%bl_pitch(k) + pit_rate*elap_time
        ctrl%bl_pitch(k) = ctrl%pit_com(k)
      end do
      ctrl%last_time_pc = time
    end if

    ctrl%last_time = time
  end subroutine operate_controller

end module m_alm_controller
