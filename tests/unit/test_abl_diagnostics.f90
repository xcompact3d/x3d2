program test_abl_diagnostics
  use m_abl_diagnostics, only: friction_velocity, neutral_log_law
  use m_common, only: dp

  implicit none

  real(dp), parameter :: tolerance = 100._dp*epsilon(1._dp)
  real(dp) :: expected, actual
  logical :: all_pass

  all_pass = .true.

  expected = 0.45_dp/0.41_dp*log(10._dp/0.1_dp)
  actual = neutral_log_law(10._dp, 0.45_dp, 0.41_dp, 0.1_dp)
  call check_close('neutral log law', actual, expected, tolerance, all_pass)

  actual = friction_velocity(-9._dp, 0._dp)
  call check_close('streamwise friction velocity', actual, 3._dp, &
                   tolerance, all_pass)

  actual = friction_velocity(-3._dp, -4._dp)
  call check_close('vector friction velocity', actual, sqrt(5._dp), &
                   tolerance, all_pass)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'

contains

  subroutine check_close(label, value, reference, tol, pass)
    character(*), intent(in) :: label
    real(dp), intent(in) :: value, reference, tol
    logical, intent(inout) :: pass

    if (abs(value - reference) > tol) then
      print *, 'FAIL: ', label, value, reference
      pass = .false.
    else
      print *, 'PASS: ', label
    end if
  end subroutine check_close

end program test_abl_diagnostics
