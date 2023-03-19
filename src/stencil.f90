module m_stencil
   implicit none
   type :: stencil
      integer :: order
      integer, allocatable :: nodes(:)
      real, allocatable :: coeffs(:)
      real :: lower, upper
   contains
      procedure, public :: get_upper, get_lower
      procedure, public :: flip, is_equal
   end type stencil

contains

   pure elemental type(stencil) function flip(self)
      class(stencil), intent(in) :: self
      integer, allocatable :: nodes(:)
      real, allocatable :: coeffs(:)

      nodes = self%nodes
      coeffs = self%coeffs
      flip = stencil( &
           & order=self%order, &
           & nodes=-nodes, &
           & coeffs=coeffs, &
           & lower=self%lower, upper=self%upper &
           & )
      if (self%order == 2) then
         flip%coeffs = -1.*flip%coeffs
      end if
   end function flip

   pure elemental real function get_upper(self)
      class(stencil), intent(in) :: self
      get_upper = self%upper
   end function get_upper

   pure elemental real function get_lower(self)
      class(stencil), intent(in) :: self
      get_lower = self%lower
   end function get_lower

   pure logical function is_equal(self, st, tol)
      class(stencil), intent(in) :: self
      type(stencil), intent(in) :: st
      real, intent(in) :: tol
      logical nodes_equal, coeffs_equal
      nodes_equal = all(self%nodes == st%nodes)
      coeffs_equal = all(abs(self%coeffs - st%coeffs) < tol)

      is_equal = nodes_equal .and. coeffs_equal
   end function is_equal

end module m_stencil
