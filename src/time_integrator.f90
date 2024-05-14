module m_time_integrator
  use m_allocator, only: allocator_t, field_t, flist_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X

  implicit none

  private adams_bashforth_1st, adams_bashforth_2nd
  private adams_bashforth_3rd, adams_bashforth_4th

  type :: time_intg_t
    integer :: istep, nsteps, nsubsteps, order, nvars, nolds
    type(flist_t), allocatable :: olds(:, :)
    type(flist_t), allocatable :: curr(:)
    type(flist_t), allocatable :: deriv(:)
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
  contains
    procedure :: step
    procedure :: adams_bashforth_1st
    procedure :: adams_bashforth_2nd
    procedure :: adams_bashforth_3rd
    procedure :: adams_bashforth_4th
  end type time_intg_t

  interface time_intg_t
    module procedure constructor
  end interface time_intg_t

contains

  function constructor(backend, allocator, nvars)
    implicit none

    type(time_intg_t) :: constructor
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    integer, intent(in), optional :: nvars

    integer :: i, j

    constructor%backend => backend
    constructor%allocator => allocator

    if (present(nvars)) then
      constructor%nvars = nvars
    else
      constructor%nvars = 3
    end if

    constructor%istep = 0
    constructor%order = 1
    constructor%nolds = constructor%order - 1

    allocate (constructor%olds(constructor%nvars, constructor%nolds))
    allocate (constructor%curr(constructor%nvars))
    allocate (constructor%deriv(constructor%nvars))

    ! Request all the storage for old timesteps
    do i = 1, constructor%nvars
      do j = 1, constructor%nolds
        constructor%olds(i, j)%ptr => allocator%get_block(DIR_X)
      end do
    end do

  end function constructor

  subroutine step(self, u, v, w, du, dv, dw, dt)
    implicit none

    class(time_intg_t), intent(inout) :: self
    class(field_t), target, intent(inout) :: u, v, w
    class(field_t), target, intent(in) :: du, dv, dw

    real(dp), intent(in) :: dt

    integer :: order

    ! assign pointer to variables
    self%curr(1)%ptr => u
    self%curr(2)%ptr => v
    self%curr(3)%ptr => w

    ! assign pointer to variables
    self%deriv(1)%ptr => du
    self%deriv(2)%ptr => dv
    self%deriv(3)%ptr => dw

    order = min(self%istep + 1, self%order)
    select case (order)
    case (1)
      call self%adams_bashforth_1st(dt)
    case (2)
      call self%adams_bashforth_2nd(dt)
    case (3)
      call self%adams_bashforth_3rd(dt)
    case (4)
      call self%adams_bashforth_4th(dt)
    end select

    ! increment step counter
    self%istep = self%istep + 1
  end subroutine step

  subroutine adams_bashforth_1st(self, dt)
    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i

    do i = 1, self%nvars
      call self%backend%vecadd(dt, self%deriv(i)%ptr, 1._dp, self%curr(i)%ptr)

      ! for startup
      if (self%istep == 0 .and. self%order > 1) then
        ! update olds(1) with new derivative
        call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                                 self%olds(i, 1)%ptr)
      end if
    end do

  end subroutine adams_bashforth_1st

  subroutine adams_bashforth_2nd(self, dt)
    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i
    class(field_t), pointer :: ptr

    do i = 1, self%nvars
      call self%backend%vecadd(1.5_dp*dt, self%deriv(i)%ptr, 1._dp, &
                               self%curr(i)%ptr)
      call self%backend%vecadd(-0.5_dp*dt, self%olds(i, 1)%ptr, 1._dp, &
                               self%curr(i)%ptr)

      ! for startup
      if (self%istep == 1 .and. self%order > 2) then
        ! rotate pointers
        call rotate(self%olds(i, :), 2)
      end if

      ! update olds(1) with new derivative
      call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                               self%olds(i, 1)%ptr)
    end do

  end subroutine adams_bashforth_2nd

  subroutine adams_bashforth_3rd(self, dt)
    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i
    class(field_t), pointer :: ptr

    do i = 1, self%nvars
      ! update solution
      call self%backend%vecadd(23._dp/12._dp*dt, self%deriv(i)%ptr, &
                               1._dp, self%curr(i)%ptr)
      call self%backend%vecadd(-4._dp/3._dp*dt, self%olds(i, 1)%ptr, &
                               1._dp, self%curr(i)%ptr)
      call self%backend%vecadd(5._dp/12._dp*dt, self%olds(i, 2)%ptr, &
                               1._dp, self%curr(i)%ptr)

      ! for startup
      if (self%istep == 2 .and. self%order > 3) then
        ! rotate pointers
        call rotate(self%olds(i, :), 3)
        ! after startup
      else
        ! rotate pointers
        call rotate(self%olds(i, :), 2)
      end if

      ! update olds(1) with new derivative
      call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                               self%olds(i, 1)%ptr)
    end do

  end subroutine adams_bashforth_3rd

  subroutine adams_bashforth_4th(self, dt)
    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i
    class(field_t), pointer :: ptr

    do i = 1, self%nvars
      ! update solution
      call self%backend%vecadd(55._dp/24._dp*dt, self%deriv(i)%ptr, &
                               1._dp, self%curr(i)%ptr)
      call self%backend%vecadd(-59._dp/24._dp*dt, self%olds(i, 1)%ptr, &
                               1._dp, self%curr(i)%ptr)
      call self%backend%vecadd(37._dp/24._dp*dt, self%olds(i, 2)%ptr, &
                               1._dp, self%curr(i)%ptr)
      call self%backend%vecadd(3._dp/8._dp*dt, self%olds(i, 3)%ptr, &
                               1._dp, self%curr(i)%ptr)

      ! for startup
      if (self%istep == 3 .and. self%order > 4) then
        ! rotate pointers
        call rotate(self%olds(i, :), 4)
        ! after startup
      else
        ! rotate pointers
        call rotate(self%olds(i, :), 3)
      end if

      ! update olds(1) with new derivative
      call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                               self%olds(i, 1)%ptr)
    end do

  end subroutine adams_bashforth_4th

  subroutine rotate(sol, n)
    type(flist_t), intent(inout) :: sol(:)
    integer, intent(in) :: n

    integer :: i
    class(field_t), pointer :: ptr

    ! rotate pointer
    ptr => sol(n)%ptr
    do i = n, 2, -1
      sol(n)%ptr => sol(n - 1)%ptr
    end do
    sol(1)%ptr => ptr

  end subroutine rotate
end module m_time_integrator
