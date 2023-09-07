program test_stencil
  use iso_fortran_env, only: stderr => error_unit

  use m_stencil, only: stencil
  use m_stencil_definitions, only: get_stencil, get_boundary_stencils

  implicit none

  logical :: allpass = .true.
  real, parameter :: tol = 0.001

  type(stencil) :: test_stencil_order_1, test_stencil_order_2
  type(stencil) :: result_stencil, expected_stencil
  type(stencil) :: result_bd_stencil(2), expected_bd_stencil(2)

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

  !! Test that we can get stencils properly
  expected_stencil = stencil( &
       & order = 2, &
       & nodes = [-2, -1, 0, 1, 2], &
       & coeffs = [3. / 44., 12. / 11., &
       & - 2. * (12. / 11. + 3. / 44.), &
       & 12. / 11., 3. / 44.], &
       & upper = 2. / 11., &
       & lower = 2. / 11. &
       & )
  result_stencil = get_stencil("compact6", 2)
  if (.not. expected_stencil%is_equal(result_stencil, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Getting bulk stencil... failed'
  else
     write(stderr, '(a)') 'Getting bulk stencil... passed'
  end if

  expected_stencil = stencil( &
       & order = 2, &
       & nodes = [-2, -1, 0, 1, 2], &
       & coeffs = [3. / 44., 12. / 11., &
       & - 2. * (12. / 11. + 3. / 44.), &
       & 12. / 11., 3. / 44.], &
       & upper = 2. / 11., &
       & lower = 2. / 11. &
       & )
  result_stencil = get_stencil("compact6", 2)
  if (.not. expected_stencil%is_equal(result_stencil, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Getting bulk stencil... failed'
  else
     write(stderr, '(a)') 'Getting bulk stencil... passed'
  end if

  expected_bd_stencil(1) = stencil( &
       & order = 1, &
       & nodes = [0, -1, -2, -3], &
       & coeffs = [-5. / 2., 2., 0.5, 0.], &
       & lower = 0., upper = 2. &
       & )
  expected_bd_stencil(2) = stencil( &
       & order = 1, &
       & nodes = [+1, 0, -1, -2], &
       & coeffs = [-3. / 4., 0., 3. / 4., 0.], &
       & lower = 1. / 4., upper = 1. / 4. &
       & )
  result_bd_stencil = get_boundary_stencils("dirichlet", 1, right=.true.)
  if (.not. expected_bd_stencil(1)%is_equal(result_bd_stencil(1), tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Getting boundary stencil... failed'
  else
     write(stderr, '(a)') 'Getting boundary stencil... passed'
  end if

  if (allpass) then
     write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
     error stop 'SOME TESTS FAILED.'
  end if
end program test_stencil
