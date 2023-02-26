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
    integer :: i
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
       associate(left => diffengine%left_stencils, &
            right => diffengine%right_stencils, &
            bulk => diffengine%bulk_stencil)
         upper_diag = [ &
              & left%get_upper(), &
              & (bulk%get_upper(), i=size(left), length - size(right)), &
              & reverse(right%get_upper()) &
              ]
         lower_diag = [ &
              & left%get_upper(), &
              & (bulk%get_lower(), i = size(left), length - size(right)), &
              & reverse(right%get_upper()) &
              ]
       end associate
       if (size(upper_diag) /= length) then
          error stop "upper-diag has the wrong size"
       end if
       if (size(lower_diag) /= length) then
          error stop "lower-diag has the wrong size"
       end if
       diffengine%tomsolv = tridiagsolv(lower_diag, upper_diag)
    end if

  end function diffengine_constructor


  subroutine diff(self, f, df, dx)
    class(diffengine_t), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(in) :: dx
    real, intent(out) :: df(:, :)
    integer :: i ! Loop counter
    integer :: pos ! current position along differentiation axis
    integer :: n ! Size of the differentiation pencil
    integer, allocatable :: lnodes_1(:), lnodes_2(:)
    integer, allocatable :: rnodes_1(:), rnodes_2(:)

    n = size(f, 2)

    pos = 1
    lnodes_1 = modulo(self%left_stencils(1)%nodes + (pos - 1), n) + 1
    pos = 2
    lnodes_2 = modulo(self%left_stencils(2)%nodes + (pos - 1), n) + 2
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
    rnodes_1 = modulo(self%right_stencils(1)%nodes + (pos - 1), n) + 1
    pos = n
    rnodes_2 = modulo(self%right_stencils(2)%nodes + (pos - 1), n) + 2
    !$omp simd
    do i = 1, size(f, 1)
       df(i, n - 1) = dot_product(self%right_stencils(1)%coeffs, f(i, rnodes_1))
       df(i, n) = dot_product(self%right_stencils(2)%coeffs, f(i, rnodes_2))
    end do
    !$omp end simd

    df = df / dx

    call self%tomsolv%solve(df, df)
  end subroutine diff

  pure function reverse(x)
    real, intent(in) :: x(:)
    real, allocatable :: reverse(:)
    reverse = x(size(x):1:-1)
  end function reverse

end module m_diffengine
