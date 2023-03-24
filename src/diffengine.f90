module m_diffengine
  !! This module defines the diffengine_t derived type representing a
  !! diffentiation engine for a given set of spatial stencils (both at
  !! the boundaries and in the bulk of the domain).
   use m_tridiagsolv, only: tridiagsolv, periodic_tridiagsolv
   use m_stencil, only: stencil
   use m_stencil_definitions, only: get_boundary_stencils, get_stencil

   implicit none

   type :: diffengine_t
      !! Concrete derived type representing a diffenrentiation engine.
      !! The purpose of a differentiation engine is to make available
      !! a procedure `diff` performing differentiation along a data
      !! array, accoring to a set of boundary and bulk differentiation
      !! stencils:
      !!
      !! L - L - O - O - O - O - O - O - O - O - R - R
      !!
      !! where 'L' represent the two left side boundary points and 'R'
      !! the two right side boundary points.
      !!
      !! A differentiation engine is instanciated through the custom
      !! type constructor [[m_diffengine:diffengine_constructor]]:
      !!
      !! ```
      !! diffengine = diffengine_t(bulk_key, length, order, dx, &
      !!                  & left_key, right_key)
      ! ```
      !!
      !! where bulk_key, left_key, and right_key are character strings
      !! specifying the stencil type. Example are "compact6" or
      !! "dirichlet".
      !!
      !! Both the `left_key` and `right_key` parameters are optional.
      !! If they are omitted, the derivation domain is assumed to be
      !! periodic, i.e.
      !!
      !! O - O - O - O - O - O - O - O - O - O - O - O
      !!
      type(stencil) :: bulk_stencil
      type(stencil) :: left_stencils(2)
      type(stencil) :: right_stencils(2)
      class(tridiagsolv), allocatable :: tomsolv
      real :: dx
   contains
      procedure, public :: diff
   end type diffengine_t

   interface diffengine_t
      module procedure diffengine_constructor
   end interface diffengine_t

contains

   function diffengine_constructor(bulk_key, length, order, dx, &
        & left_key, right_key) result(diffengine)
      !> Stencil key for bulk points, can be one of `"compact6"`.
      character(*), intent(in) :: bulk_key
      !> Stencil key for left boundary points, can be one of
      !> `"dirichlet"`.
      character(*), optional :: left_key
      !> Stencil key for right boundary points, can be one of
      !> `"dirichlet"`.
      character(*), optional :: right_key
      !> Length of the derivation domain
      integer, intent(in) :: length
      !> Differentiation order, can be either `1` or `2`.
      integer, intent(in) :: order
      !> Mesh step size
      real, intent(in) :: dx
      type(diffengine_t) :: diffengine
      integer :: i
      logical :: is_periodic
      real, allocatable :: lower_diag(:), upper_diag(:)

      diffengine%dx = dx

      diffengine%bulk_stencil = get_stencil(bulk_key, order)

      if (present(left_key) .and. present(right_key)) then
         is_periodic = .false.
         diffengine%left_stencils = get_boundary_stencils(left_key, order)
         diffengine%right_stencils = get_boundary_stencils(right_key, &
              & order, right=.true.)
      else if ((.not. present(left_key)) .and. &
           &   (.not. present(right_key))) then
         is_periodic = .true.
         diffengine%left_stencils(:) = diffengine%bulk_stencil
         diffengine%right_stencils(:) = diffengine%bulk_stencil
      else
         error stop "Both left and right boundary types must be specified."
      end if

      ! In the periodic case the first element of the differentiation
      ! pencil (f(1)) is the same as the last (f(n)).  Alhtough we have
      ! n data, the effective domain size if n - 1.
      ! See periodic_tridiagsolv%solve_periodic (thomas.f90)
      if (is_periodic) then
         upper_diag = [(diffengine%bulk_stencil%get_upper(), i=1, length - 2)]
         lower_diag = [(diffengine%bulk_stencil%get_lower(), i=1, length - 2)]
         diffengine%tomsolv = periodic_tridiagsolv(lower_diag, upper_diag)
      else
         associate (left => diffengine%left_stencils, &
                    right => diffengine%right_stencils, &
                    bulk => diffengine%bulk_stencil)
            upper_diag = [ &
                 & left%get_upper(), &
                 & (bulk%get_upper(), i=size(left), length - size(right)), &
                 & reverse(right%get_upper()) &
                 ]
            lower_diag = [ &
                 & left%get_upper(), &
                 & (bulk%get_lower(), i=size(left), length - size(right)), &
                 & reverse(right%get_upper()) &
                 ]
         end associate
         if (size(upper_diag) /= length) then
            error stop "upper-diag has the wrong size"
         end if
         if (size(lower_diag) /= length) then
            error stop "lower-diag has the wrong size"
         end if
         upper_diag = upper_diag(1:(length - 1))
         lower_diag = lower_diag(2:length)
         diffengine%tomsolv = tridiagsolv(lower_diag, upper_diag)
      end if

   end function diffengine_constructor

   subroutine diff(self, f, df)
     !! Compute derivatives along the second dimension of an input 2D
     !! array using compact finite-differences stencils. The stencils
     !! are first applied before computing the solution by solving the
     !! associated tridiagonal systems.
     !!
     !!              -->
     !!               _
     !! O - O - O - | O | - O - O - O - O - O - O - O - O
     !! O - O - O - | O | - O - O - O - O - O - O - O - O
     !! O - O - O - | O | - O - O - O - O - O - O - O - O
     !! O - O - O - | O | - O - O - O - O - O - O - O - O
     !!               -
     !!
      class(diffengine_t), intent(in) :: self
      !> Input array.  Size along second dimension should be equal to
      !> the `allocator_t` instance `length` component.
      real, intent(in) :: f(:, :)
      !> Output array. Shape is identical to the input array's.
      real, intent(out) :: df(:, :)
      integer :: i ! Loop counter
      integer :: pos ! current position along differentiation axis
      integer :: n ! Size of the differentiation pencil
      integer, allocatable :: lnodes_1(:), lnodes_2(:)
      integer, allocatable :: rnodes_1(:), rnodes_2(:)

      n = size(f, 2)

      pos = 1
      lnodes_1 = modulo(self%left_stencils(1)%nodes + (pos - 1), n - 1) + 1
      pos = 2
      lnodes_2 = modulo(self%left_stencils(2)%nodes + (pos - 1), n - 1) + 1
      !$omp simd
      do i = 1, size(f, 1)
         df(i, 1) = dot_product(self%left_stencils(1)%coeffs, f(i, lnodes_1))
         df(i, 2) = dot_product(self%left_stencils(2)%coeffs, f(i, lnodes_2))
      end do

      !$omp end simd
      do pos = 3, n - 2
         !$omp simd
         do i = 1, size(f, 1)
            df(i, pos) = dot_product(self%bulk_stencil%coeffs, &
                 & f(i, self%bulk_stencil%nodes + pos))
         end do
         !$omp end simd
      end do

      pos = n - 1
      rnodes_1 = modulo(self%right_stencils(1)%nodes + (pos - 1), n - 1) + 1
      pos = n
      rnodes_2 = modulo(self%right_stencils(2)%nodes + (pos - 1), n - 1) + 1
      !$omp simd
      do i = 1, size(f, 1)
         df(i, n - 1) = dot_product( &
              & self%right_stencils(1)%coeffs, &
              & f(i, rnodes_1))
         df(i, n) = dot_product( &
              & self%right_stencils(2)%coeffs, &
              & f(i, rnodes_2))
      end do
      !$omp end simd

      df = df / self%dx

      call self%tomsolv%solve(df, df)
   end subroutine diff

   pure function reverse(x)
      real, intent(in) :: x(:)
      real, allocatable :: reverse(:)
      reverse = x(size(x):1:-1)
   end function reverse

end module m_diffengine
