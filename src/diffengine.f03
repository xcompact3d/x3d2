module m_diffengine
  use m_tridiagsolv
  use m_stencil

  implicit none

  type :: diffengine
     type(stencil) :: bulk_stencil
     type(stencil) :: left_stencils(2)
     type(stencil) :: right_stencils(2)
     class(tridiagsolv), allocatable :: tomsolv
  end type diffengine

  type(tri

contains

  function diffengine_constructor(bulk_stencil, left_stencils, right_stencils, n) &
       & result(diffengine)
    type(stencil), intent(in) :: bulk_stencil
    type(stencil), intent(in) :: left_stencils(2), &
         & right_stencils(2)
    integer, intent(in) :: n

    if (present(left_stencils) .and. present(right_stencils)) then
       upper_diag = [ &
            & left_stencils%get_upper(), &
            & (bulk_stencil%get_upper(), i= 3, n - 2), &
            & reverse(right_stencils%get_upper()) &
            ]
       lower_diag = [ &
            & left_stencils%get_upper(), &
            & (bulk_stencil%get_lower(), i = 3, n), &
            & reverse(right_stencils%get_upper()) &
            ]
       diffengine%tomsolv = tridiagsolv(lower_diag, upper_diag, n)
    else if ((.not. present(left_stencil)) .and. &
         &   (.not. present(right_stencil)) then
       diffengine%tomsolv = periodic_tridiagsolv( &
            & bulk_stencil%lower, bulk_stencil%upper, n)
    end if
  end function diffengine_constructor


  subroutine diff(self, f, df)
    class(diffengine), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    integer, parameter :: n = size(f, 2)
    type(stencil), pointer :: s
    !$omp simd
    do i = 1, size(f, 1)
       s => self%left_stencils(1)
       du(i, 1) = dot_product(s%coeffs, u(i, s%nodes + 1))
       s => self%stencils(2)
       du(i, 2) = dot_product(s%coeffs, u(i, s%nodes + 2))
    end do

    !$omp end simd
    s => self%bulk_stencil
    do j = 3, n - 2
       !$omp simd
       do i = 1, size(f, 1)
          du(i, j) = dot_product(s%coeffs, u(i, s%nodes + j))
       end do
       !$omp end simd
    end do
    !$omp simd
    do i = 1, size(f, 1)
       s => self%stencils(1)
       du(i, n - 1) = dot_product(s%coeffs, u(i, s%nodes + n - 1))
       s => self%stencils(2)
       du(i, n) = dot_product(s%coeffs, u(i, s%nodes + n))
    end do
    !$omp end simd

    call self%tomsolv%solve(du, df)
  end subroutine diff

end module m_diffengine
