module m_diffengine
  use m_tridiagsolv, only: tridiagsolv, periodic_tridiagsolv
  use m_stencil, only: stencil
  use m_stencil_definitions, only: get_boundary_stencils, get_stencil

  implicit none

  type :: diffengine_t
     type(stencil) :: bulk_stencil
     type(stencil) :: left_stencils(2)
     type(stencil) :: right_stencils(2)
     class(tridiagsolv), allocatable :: tomsolv
   contains
     procedure, public :: diff
  end type diffengine_t

  interface diffengine_t
     module procedure diffengine_constructor
  end interface diffengine_t

contains

  function diffengine_constructor(bulk_key, length, order, left_key, right_key) &
       & result(diffengine)
    character(*), intent(in) :: bulk_key
    character(*), optional :: left_key, right_key
    integer, intent(in) :: length, order
    type(diffengine_t) :: diffengine
    integer :: i, n
    logical :: is_periodic
    real :: lower_diag(length), upper_diag(length)


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

    if (is_periodic) then
       upper_diag = [(diffengine%bulk_stencil%get_upper(), i = 1, length)]
       lower_diag = [(diffengine%bulk_stencil%get_lower(), i = 1, length)]
       diffengine%tomsolv = periodic_tridiagsolv(lower_diag, upper_diag)
    else
       upper_diag = [ &
            & diffengine%left_stencils%get_upper(), &
            & (diffengine%bulk_stencil%get_upper(), i= 3, length - 2), &
            & reverse(diffengine%right_stencils%get_upper()) &
            ]
       lower_diag = [ &
            & diffengine%left_stencils%get_upper(), &
            & (diffengine%bulk_stencil%get_lower(), i = 3, n), &
            & reverse(diffengine%right_stencils%get_upper()) &
            ]
       diffengine%tomsolv = tridiagsolv(lower_diag, upper_diag)
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

  pure function reverse(x)
    real, intent(in) :: x(:)
    real, allocatable :: reverse(:)
    reverse = x(size(x):1:-1)
  end function reverse

end module m_diffengine
