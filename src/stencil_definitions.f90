module m_stencil_definitions
  use m_stencil, only: stencil
  implicit none

  type(stencil), private :: dirichlet(2, 2)
  type(stencil), private :: compact6(2)
  logical :: stencils_defined = .false.


contains

  subroutine define_stencils()
    compact6(1) = stencil( &
         & order = 1, &
         & nodes = [-2, -1, 1, 2], &
         & coeffs = [- 1. / 36., -7. / 9., + 7. / 9., + 1. / 36.], &
         & upper = 1. / 3., &
         & lower = 1. / 3. &
         & )

    compact6(2) = stencil( &
         & order = 2, &
         & nodes = [-2, -1, 0, 1, 2], &
         & coeffs = [3. / 44., 12. / 11., &
         & - 2. * (12. / 11. + 3. / 44.), &
         & 12. / 11., 3. / 44.], &
         & upper = 2. / 11., &
         & lower = 2. / 11. &
         & )

    dirichlet(1, 1) = stencil( &
         & order = 1, &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [-5. / 2., 2., 0.5, 0.], &
         & lower = 0., upper = 2. &
         & )
    dirichlet(2, 1) = stencil( &
         & order = 1, &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [-3. / 4., 0., 3. / 4., 0.], &
         & lower = 1. / 4., upper = 1. / 4. &
         & )

    dirichlet(1, 2) = stencil( &
         & order = 2, &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [13., -27., 15., -1.], &
         & lower = 0., upper = 11. &
         & )
    dirichlet(2, 2) = stencil( &
         & order = 2, &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [6. / 5., -12. / 5., 6. / 5., 0.], &
         & lower = 1. / 10., upper = 1. / 10. &
         & )
    stencils_defined = .true.
  end subroutine define_stencils

  function get_stencil(key, order) result(s)
    integer, intent(in) :: order
    character(*), intent(in) :: key
    type(stencil) :: s
    if (.not. stencils_defined) then
       call define_stencils()
    endif
    if (.not. any([1, 2] == order)) then
       error stop "order must be 1 or 2"
    end if

    if (key == "compact6") then
       s = compact6(order)
    else
       error stop "Unknown key"
    end if
  end function get_stencil

  !! Need another function because boundary stencils are made of two
  !! components, so the function returns a type(stencil) array of size
  !! 2.
  function get_boundary_stencils(key, order, right) result(s)
    character(*), intent(in) :: key
    integer, intent(in) :: order
    logical, optional, intent(in) :: right
    type(stencil) :: s(2)
    if (.not. stencils_defined) then
       call define_stencils()
    endif
    if (.not. any([1, 2] == order)) then
       error stop "order must be 1 or 2"
    end if
    if (key == "dirichlet") then
       s = dirichlet(:, order)
    else
       error stop "Unknown key"
    end if
    if (present(right) .and. right) then
       s = s%flip()
    end if
  end function get_boundary_stencils
end module m_stencil_definitions
