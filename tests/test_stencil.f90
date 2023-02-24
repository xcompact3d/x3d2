program test_stencil
  use iso_fortran_env, only: stderr => error_unit

  use m_stencil, only: stencil
  use m_stencil_definitions, only: get_stencil, get_boundary_stencils

  implicit none

  logical :: allpass
  real, parameter :: tol = 0.001

  type(stencil) :: test_stencil_order_1, test_stencil_order_2
  type(stencil) :: result_stencil, expected_stencil

  test_stencil_order_1 = stencil( &
       & order = 1, &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = [2., 1., 2., 1.], &
       & lower = 0., upper = 0. &
       & )
  test_stencil_order_2 = stencil( &
       & order = 2, &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = [2., 1., 2., 1.], &
       & lower = 0., upper = 0. &
       & )

  !! Test that stencil is flipped correctly (order 1)
  expected_stencil = stencil( &
       & order = 1, &
       & nodes = [+1, 0, -1, -2], &
       & coeffs = [2., 1., 2., 1.], &
       & lower = 0., upper = 0. &
       & )
  result_stencil = test_stencil_order_1%flip()
  if (.not. expected_stencil%is_equal(result_stencil, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is flipped correctly (order 1)... failed'
  else
     write(stderr, '(a)') 'Stencil is flipped correctly (order 1)... passed'
  end if

  !! Test that stencil is flipped correctly (order 2)
  expected_stencil = stencil( &
       & order = 2, &
       & nodes = [+1, 0, -1, -2], &
       & coeffs = [-2., -1., -2., -1.], &
       & lower = 0., upper = 0. &
       & )
  result_stencil = test_stencil_order_2%flip()
  if (.not. expected_stencil%is_equal(result_stencil, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is flipped correctly (order 2)... failed'
  else
     write(stderr, '(a)') 'Stencil is flipped correctly (order 2)... passed'
  end if
end program test_stencil
