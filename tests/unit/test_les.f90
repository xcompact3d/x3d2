program test_les
  use m_common, only: dp
  use m_config, only: les_config_t
  use m_les, only: les_t, filter_width, strain_rate_magnitude, &
                   wall_damped_mixing_length

  implicit none

  type(les_config_t) :: config
  type(les_t) :: les
  real(dp) :: gradient(3, 3), spacing(3), expected, actual, length
  real(dp), parameter :: tol = 100._dp*epsilon(1._dp)
  logical :: all_pass

  all_pass = .true.
  spacing = [0.5_dp, 0.25_dp, 0.125_dp]

  call config%read(nml_string='&les_params model="smagorinsky", '// &
    'smagorinsky_constant=0.2, wall_damping=.false. /')
  les = les_t(config)

  ! Simple shear u=gamma*y has |S|=|gamma|.
  gradient = 0._dp
  gradient(1, 2) = 4._dp
  expected = (0.2_dp*filter_width(spacing))**2*4._dp
  actual = les%nut_from_gradient(gradient, spacing)
  call check_close('simple-shear nut', actual, expected, tol, all_pass)
  call check_close('simple-shear strain', strain_rate_magnitude(gradient), &
                   4._dp, tol, all_pass)

  ! Solid-body rotation has an antisymmetric gradient and no strain.
  gradient = 0._dp
  gradient(1, 2) = -3._dp
  gradient(2, 1) = 3._dp
  actual = les%nut_from_gradient(gradient, spacing)
  call check_close('solid-body rotation', actual, 0._dp, tol, all_pass)

  ! Mason-Thomson damping reaches zero at a smooth wall and asymptotes to Cs*Delta.
  length = wall_damped_mixing_length(1._dp, 0._dp, 0.2_dp, &
                                      0.4_dp, 3._dp, 0._dp)
  call check_close('wall damping at wall', length, 0._dp, tol, all_pass)
  length = wall_damped_mixing_length(1._dp, 1.0e6_dp, 0.2_dp, &
                                      0.4_dp, 3._dp, 0._dp)
  call check_close('wall damping far from wall', length, 0.2_dp, &
                   tol, all_pass)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'

contains

  subroutine check_close(label, value, reference, tolerance, pass)
    character(*), intent(in) :: label
    real(dp), intent(in) :: value, reference, tolerance
    logical, intent(inout) :: pass

    if (abs(value - reference) > tolerance) then
      print *, 'FAIL: ', label, value, reference
      pass = .false.
    else
      print *, 'PASS: ', label
    end if
  end subroutine check_close

end program test_les
