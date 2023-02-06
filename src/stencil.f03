module m_stencil
  implicit none
  type :: stencil
     integer, allocatable :: nodes(:)
     real, allocatable :: coeffs(:)
     real :: lower, upper
   contains
     procedure, private :: stencil_mul_real
     procedure, public :: apply => apply_stencil
     procedure, public :: apply_along => apply_stencil_along
     procedure, public :: get_upper, get_lower
     procedure, public :: flip, is_equal
     generic :: operator(*) => stencil_mul_real
  end type stencil

contains

  pure logical function is_equal(self, st, tol)
    class(stencil), intent(in) :: self
    type(stencil), intent(in) :: st
    real, intent(in) :: tol
    logical nodes_equal, coeffs_equal
    nodes_equal = all(self%nodes == st%nodes)
    coeffs_equal = all(abs(self%coeffs - st%coeffs) < tol)

    is_equal = nodes_equal .and. coeffs_equal
  end function is_equal

  pure elemental type(stencil) function stencil_mul_real(self, a)
    class(stencil), intent(in) :: self
    real, intent(in) :: a
    integer, allocatable :: nodes(:)
    real, allocatable :: coeffs(:)

    nodes = self%nodes
    coeffs = self%coeffs
    stencil_mul_real = stencil( &
         & nodes = nodes, &
         & coeffs = a * coeffs, &
         & lower = self%lower, upper = self%upper &
         & )
  end function stencil_mul_real

  pure elemental type(stencil) function flip(self)
    class(stencil), intent(in) :: self
    integer, allocatable :: nodes(:)
    real, allocatable :: coeffs(:)

    nodes = self%nodes
    coeffs = self%coeffs
    flip = stencil( &
         & nodes = - nodes, &
         & coeffs = coeffs, &
         & lower = self%lower, upper = self%upper &
         & )
  end function flip

  pure elemental real function get_upper(self)
    class(stencil), intent(in) :: self
    get_upper = self%upper
  end function get_upper

  pure elemental real function get_lower(self)
    class(stencil), intent(in) :: self
    get_lower = self%lower
  end function get_lower

  pure real function apply_stencil(self, f, ref)
    class(stencil), intent(in) :: self
    real, intent(in) :: f(:)
    integer, intent(in) :: ref
    real, allocatable :: eval(:)

    eval = f(self%nodes + ref)
    apply_stencil = dot_product(eval, self%coeffs)
  end function apply_stencil

  pure function apply_stencil_along(self, f)
    class(stencil), intent(in) :: self
    real, intent(in) :: f(:)
    real, allocatable :: apply_stencil_along(:), f_padded(:), eval(:)
    integer :: ref, i, lpad, rpad

    lpad = minval(self%nodes) + 1
    rpad = size(f) + maxval(self%nodes)
    allocate(f_padded(lpad:rpad))
    do i = lpad, 0
       f_padded(i) = f(size(f) + i - 1)
    end do
    do i = size(f) + 1, rpad
       f_padded(i) = f(i - size(f) + 1)
    end do
    f_padded(1:size(f)) = f
    allocate(apply_stencil_along, source = f)
    do ref = 1, size(f)
       ! Would be nice to reuse self%apply_stencil here but any
       ! negative indices are lost when passing f_padded
       eval = f_padded(self%nodes + ref)
       apply_stencil_along(ref) = dot_product(eval, self%coeffs)
    end do
  end function apply_stencil_along

end module m_stencil
