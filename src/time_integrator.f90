module m_time_integrator
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X
  use m_field, only: field_t, flist_t

  implicit none

  private adams_bashforth, runge_kutta

  type :: time_intg_t
    integer :: method, istep, istage, order, nstep, nstage, nvars, nolds
    real(dp) :: coeffs(4, 4)
    real(dp) :: rk_b(4, 4)
    real(dp) :: rk_a(3, 3, 4)
    character(len=3) :: sname
    type(flist_t), allocatable :: olds(:, :)
    type(flist_t), allocatable :: curr(:)
    type(flist_t), allocatable :: deriv(:)
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    procedure(stepper_func), pointer :: stepper => null()
  contains
    procedure :: finalize
    procedure :: step
    procedure :: runge_kutta
    procedure :: adams_bashforth
  end type time_intg_t

  interface time_intg_t
    module procedure init
  end interface time_intg_t

  abstract interface
    subroutine stepper_func(self, dt)
      import :: time_intg_t
      import :: dp
      implicit none

      class(time_intg_t), intent(inout) :: self
      real(dp), intent(in) :: dt
    end subroutine stepper_func
  end interface

contains

  subroutine finalize(self)
    implicit none

    !type(time_intg_t), intent(inout) :: self
    class(time_intg_t), intent(inout) :: self

    integer :: i, j

    ! Release all the storage for old timesteps
    do i = 1, self%nvars
      do j = 1, self%nolds
        call self%allocator%release_block(self%olds(i, j)%ptr)
      end do
    end do

    ! deallocate memory
    deallocate (self%olds)
    deallocate (self%curr)
    deallocate (self%deriv)

    print *, self%sname, ' time integrator deallocated'

  end subroutine finalize

  function init(backend, allocator, method, nvars)
    implicit none

    type(time_intg_t) :: init
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    character(3), intent(in) :: method
    integer, intent(in), optional :: nvars

    integer :: i, j, stat

    ! initialize Runge-Kutta coefficients
    ! rk1
    init%rk_a(:, 1, 1) = [0.0_dp, 0._dp, 0._dp]
    init%rk_a(:, 2, 1) = [0.0_dp, 0._dp, 0._dp]
    init%rk_a(:, 3, 1) = [0.0_dp, 0._dp, 0._dp]
    init%rk_b(:, 1) = [1._dp, 0._dp, 0._dp, 0._dp]

    ! rk2
    init%rk_a(:, 1, 2) = [0.5_dp, 0._dp, 0._dp]
    init%rk_a(:, 2, 2) = [0.0_dp, 0._dp, 0._dp]
    init%rk_a(:, 3, 2) = [0.0_dp, 0._dp, 0._dp]
    init%rk_b(:, 2) = [0._dp, 1._dp, 0._dp, 0._dp]

    ! rk3
    init%rk_a(:, 1, 3) = [0.5_dp, 0._dp, 0._dp]
    init%rk_a(:, 2, 3) = [0.0_dp, 3._dp/4._dp, 0._dp]
    init%rk_a(:, 3, 3) = [0.0_dp, 0._dp, 0._dp]
    init%rk_b(:, 3) = &
      [2._dp/9.0_dp, 1._dp/3._dp, 4._dp/9._dp, 0._dp]

    ! rk4
    init%rk_a(:, 1, 4) = [0.5_dp, 0._dp, 0._dp]
    init%rk_a(:, 2, 4) = [0._dp, 0.5_dp, 0._dp]
    init%rk_a(:, 3, 4) = [0._dp, 0._dp, 1._dp]
    init%rk_b(:, 4) = &
      [1._dp/6._dp, 1._dp/3._dp, 1._dp/3._dp, 1._dp/6._dp]

    ! initialize Adams-Bashforth coefficients
    ! ab1
    init%coeffs(:, 1) = [1._dp, 0._dp, 0._dp, 0._dp]
    ! ab2
    init%coeffs(:, 2) = [1.5_dp, -0.5_dp, 0._dp, 0._dp]
    ! ab3
    init%coeffs(:, 3) = &
      [23._dp/12._dp, -4._dp/3._dp, 5._dp/12._dp, 0._dp]
    ! ab4
    init%coeffs(:, 4) = &
      [55._dp/24._dp, -59._dp/24._dp, 37._dp/24._dp, -3._dp/8._dp]

    ! set variables
    init%backend => backend
    init%allocator => allocator
    init%sname = method

    if (init%sname(1:2) == 'AB') then
      read (init%sname(3:3), *, iostat=stat) init%order
      if (stat /= 0) error stop 'Error reading AB integration order'
      if (init%order >= 5) error stop 'Integration order >4 is not supported'
      init%nstep = init%order
      init%nstage = 1
      init%nolds = init%nstep - 1
      init%stepper => adams_bashforth
    else if (init%sname(1:2) == 'RK') then
      read (init%sname(3:3), *, iostat=stat) init%order
      if (stat /= 0) error stop 'Error reading RK integration order'
      if (init%order >= 5) error stop 'Integration order >4 is not supported'
      init%nstep = 1
      init%nstage = init%order
      init%nolds = init%nstage
      init%stepper => runge_kutta
    else
      print *, 'Integration method '//init%sname//' is not defined'
      error stop
    end if

    if (present(nvars)) then
      init%nvars = nvars
    else
      init%nvars = 3
    end if

    init%istep = 1
    init%istage = 1

    ! allocate memory
    allocate (init%olds(init%nvars, init%nolds))
    allocate (init%curr(init%nvars))
    allocate (init%deriv(init%nvars))

    ! Request all the storage for old timesteps
    do i = 1, init%nvars
      do j = 1, init%nolds
        init%olds(i, j)%ptr => allocator%get_block(DIR_X)
      end do
    end do

  end function init

  subroutine step(self, u, v, w, du, dv, dw, dt)
    implicit none

    class(time_intg_t), intent(inout) :: self
    class(field_t), target, intent(inout) :: u, v, w
    class(field_t), target, intent(in) :: du, dv, dw

    real(dp), intent(in) :: dt

    ! assign pointer to variables
    self%curr(1)%ptr => u
    self%curr(2)%ptr => v
    self%curr(3)%ptr => w

    ! assign pointer to variables
    self%deriv(1)%ptr => du
    self%deriv(2)%ptr => dv
    self%deriv(3)%ptr => dw

    call self%stepper(dt)

  end subroutine step

  subroutine runge_kutta(self, dt)
    implicit none

    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i, j

    ! update solution
    if (self%istage == self%nstage) then
      do i = 1, self%nvars
        ! update step solution from stage derivative
        if (self%nstage > 1) then
          call self%backend%vecadd(1.0_dp, self%olds(i, 1)%ptr, 0._dp, &
                                   self%curr(i)%ptr)
        end if

        do j = 1, self%nstage - 1
          call self%backend%vecadd(self%rk_b(j, self%nstage)*dt, &
                                   self%olds(i, j + 1)%ptr, &
                                   1._dp, self%curr(i)%ptr)
        end do
        call self%backend%vecadd(self%rk_b(self%nstage, self%nstage)*dt, &
                                 self%deriv(i)%ptr, &
                                 1._dp, self%curr(i)%ptr)

      end do

      ! reset stage counter
      self%istage = 1
    else
      do i = 1, self%nvars
        ! save step initial condition
        if (self%istage == 1) then
          call self%backend%vecadd(1.0_dp, self%curr(i)%ptr, 0._dp, &
                                   self%olds(i, 1)%ptr)
        end if

        ! save stage derivative
        call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                                 self%olds(i, self%istage + 1)%ptr)

        ! update stage solution
        if (self%istage > 1) then
          call self%backend%vecadd(1.0_dp, self%olds(i, 1)%ptr, 0._dp, &
                                   self%curr(i)%ptr)
        end if
        do j = 1, self%istage
          call self%backend%vecadd(self%rk_a(j, self%istage, self%nstage)*dt, &
                                   self%olds(i, j + 1)%ptr, &
                                   1._dp, self%curr(i)%ptr)
        end do
      end do

      ! increment stage counter
      self%istage = self%istage + 1
    end if

  end subroutine runge_kutta

  subroutine adams_bashforth(self, dt)
    implicit none

    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i, j
    integer :: nstep

    nstep = min(self%istep, self%nstep)
    do i = 1, self%nvars
      ! update solution
      call self%backend%vecadd(self%coeffs(1, nstep)*dt, &
                               self%deriv(i)%ptr, &
                               1._dp, self%curr(i)%ptr)
      do j = 2, nstep
        call self%backend%vecadd(self%coeffs(j, nstep)*dt, &
                                 self%olds(i, j - 1)%ptr, &
                                 1._dp, self%curr(i)%ptr)
      end do

      ! rotate pointers
      if (nstep < self%nstep) then
        ! for startup
        if (self%istep > 1) then
          call rotate(self%olds(i, :), nstep)
        end if
      else
        ! after startup
        if (self%nstep > 2) then
          call rotate(self%olds(i, :), nstep - 1)
        end if
      end if

      ! update olds(1) with new derivative
      if (self%nstep > 1) then
        call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                                 self%olds(i, 1)%ptr)
      end if
    end do

    ! increment step counter
    self%istep = self%istep + 1

  end subroutine adams_bashforth

  subroutine rotate(sol, n)
    implicit none

    type(flist_t), intent(inout) :: sol(:)
    integer, intent(in) :: n

    integer :: i
    class(field_t), pointer :: ptr

    ! rotate pointer
    ptr => sol(n)%ptr
    do i = n, 2, -1
      sol(i)%ptr => sol(i - 1)%ptr
    end do
    sol(1)%ptr => ptr

  end subroutine rotate

end module m_time_integrator
