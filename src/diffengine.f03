module m_diffengine
  use m_tridiagsolv
  use m_stencil
  use m_stencil_definitions, only: get_boundary_stencils, get_stencil

  implicit none

  type :: diffengine
     type(stencil) :: bulk_stencil
     type(stencil), allocatable :: left_stencils(2)
     type(stencil), allocatable :: right_stencils(2)
     class(tridiagsolv), allocatable :: tomsolv
  end type diffengine

  type(tri

contains

  function diffengine_constructor(bulk_sten, n, order, left_sten, right_sten) &
       & result(diffengine)
    character(*), intent(in) :: bulk_sten
    character(*), optional :: left_sten, right_sten
    integer, intent(in) :: n, order

    type(stencil) :: bulk_stencil
    type(stencil) :: :: left_stencils(2), right_stencils(2)


    diffengine%bulk_stencil = get_stencil(1, order)

    if (present(bulk_sten) .and. present(right_sten)) then
       is_periodic = .false.
       diffengine%left_stencils = get_boundary_stencils(1, order)
       diffengine%right_stencils = get_boundary_stencils(1, order, right=.true.)
    else if ((.not. present(left_stencil)) .and. &
         &   (.not. present(right_stencil)) then
       is_periodic = .true.
    else
       error stop "Both left and right boundary types must be specified."
    end if

    if (is_periodic) then
       diffengine%tomsolv = periodic_tridiagsolv( &
            & bulk_stencil%lower, bulk_stencil%upper, n)
    else
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
