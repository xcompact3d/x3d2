module m_time_integrator
  !! Time integration schemes for temporal advancement.
  !!
  !! This module provides explicit time integration methods for advancing
  !! solutions in time. It supports two families of schemes:
  !!
  !! **1. Runge-Kutta (RK) Methods**
  !! Multi-stage schemes that achieve high-order accuracy within a single
  !! timestep. Supported orders: RK1 (Euler), RK2, RK3, RK4. Each stage
  !! requires an evaluation of the right-hand side (derivative).
  !!
  !! **2. Adams-Bashforth (AB) Methods**
  !! Multi-step schemes that use derivative information from previous
  !! timesteps to achieve high-order accuracy. Supported orders: AB1, AB2,
  !! AB3, AB4. These methods are more memory-efficient than RK schemes
  !! for the same order of accuracy.
  !!
  !! The `time_intg_t` type encapsulates all integration state and provides
  !! a unified interface through the step procedure pointer, which routes
  !! to either runge_kutta() or adams_bashforth() based on the selected method.
  !!
  !! Old timestep/stage data is stored in the `olds` array and managed
  !! automatically through rotation mechanisms for AB methods.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X
  use m_field, only: field_t, flist_t

  implicit none

  private adams_bashforth, runge_kutta

  type :: time_intg_t
    !! Time integrator for explicit multi-step and multi-stage methods.
    !!
    !! This type encapsulates all data and methods needed for time integration
    !! of ordinary differential equations (ODEs) arising from spatial discretization
    !! of the Navier-Stokes equations:
    !!
    !! \[
    !! \frac{d\mathbf{u}}{dt} = \mathbf{F}(\mathbf{u}, t)
    !! \]
    !!
    !! where \(\mathbf{F}\) represents the spatial operators (advection, diffusion,
    !! pressure gradient, etc.).
    !!
    !! **Supported Methods:**
    !!
    !! - **Adams-Bashforth (AB1-AB4)**: Explicit multi-step methods using
    !!   previous timestep derivatives. Efficient (single evaluation per step)
    !!   but requires startup procedure for higher orders.
    !! - **Runge-Kutta (RK1-RK4)**: Explicit multi-stage methods using
    !!   intermediate stages within a timestep. Self-starting but requires
    !!   multiple evaluations per step.
    !!
    !! **Method Selection:**
    !!
    !! The `step` procedure pointer is bound at initialization to either
    !! `runge_kutta()` or `adams_bashforth()` based on the method name
    !! (e.g., "AB3" or "RK4"), enabling polymorphic time stepping.
    !!
    !! **Data Management:**
    !!
    !! - **AB methods**: Store previous timestep derivatives in `olds` array,
    !!   rotated each timestep to maintain history
    !! - **RK methods**: Store intermediate stage solutions in `olds` array,
    !!   overwritten within each timestep
    !!
    !! **Startup Procedure (AB only):**
    !!
    !! Higher-order AB methods (AB2-AB4) ramp up from first-order during initial
    !! timesteps until sufficient derivative history is available.
    integer :: method       !! Integration method identifier (unused, kept for compatibility)
    integer :: istep        !! Current timestep number (for AB startup ramping)
    integer :: istage       !! Current stage number within timestep (RK only)
    integer :: order        !! Order of accuracy of the scheme (1-4)
    integer :: nstep        !! Number of timesteps needed (AB: order, RK: 1)
    integer :: nstage       !! Number of stages per timestep (AB: 1, RK: order)
    integer :: nvars        !! Number of variables being integrated
    integer :: nolds        !! Number of old derivatives/solutions to store
    real(dp) :: coeffs(4, 4)  !! Adams-Bashforth coefficients [stage, order]
    real(dp) :: rk_b(4, 4)    !! Runge-Kutta final weights [stage, order]
    real(dp) :: rk_a(3, 3, 4) !! Runge-Kutta stage weights [from_stage, to_stage, order]
    character(len=3) :: sname !! Scheme name (e.g., 'AB3', 'RK4')
    type(flist_t), allocatable :: olds(:, :) !! Old derivatives/solutions [nvars, nolds]
    class(base_backend_t), pointer :: backend    !! Computational backend for operations
    class(allocator_t), pointer :: allocator     !! Memory allocator for field storage
    procedure(stepper_func), pointer :: step => null() !! Function pointer to integration method
  contains
    procedure :: finalize       !! Clean up and release allocated memory
    procedure :: runge_kutta    !! Runge-Kutta time integration implementation
    procedure :: adams_bashforth !! Adams-Bashforth time integration implementation
  end type time_intg_t

  interface time_intg_t
    module procedure init
  end interface time_intg_t

  abstract interface
    subroutine stepper_func(self, curr, deriv, dt)
      !! Abstract interface for time stepping functions.
      !!
      !! Defines the signature for integration methods (RK or AB).
      !! Each method takes the current solution, its derivative, and
      !! the timestep size, and updates the solution accordingly.
      import :: time_intg_t
      import :: dp
      import :: flist_t
      implicit none

      class(time_intg_t), intent(inout) :: self  !! Time integrator state
      type(flist_t), intent(inout) :: curr(:)    !! Current solution variables [nvars]
      type(flist_t), intent(in) :: deriv(:)      !! Time derivatives of variables [nvars]
      real(dp), intent(in) :: dt                 !! Timestep size
    end subroutine stepper_func
  end interface

contains

  subroutine finalize(self)
    !! Finalise time integrator and release allocated resources.
    !!
    !! Releases all field storage blocks used for storing old derivatives
    !! or stage solutions, and deallocates the olds array.
    implicit none

    !type(time_intg_t), intent(inout) :: self
    class(time_intg_t), intent(inout) :: self !! Time integrator to finalise

    integer :: i, j

    ! Release all the storage for old timesteps
    do i = 1, self%nvars
      do j = 1, self%nolds
        call self%allocator%release_block(self%olds(i, j)%ptr)
      end do
    end do

    ! deallocate memory
    deallocate (self%olds)

    print *, self%sname, ' time integrator deallocated'

  end subroutine finalize

  function init(backend, allocator, method, nvars)
    !! Initialise time integrator with specified method and coefficients.
    !!
    !! This constructor configures the time integration scheme based on the
    !! method string (e.g., 'AB3' or 'RK4'). It initialises all Runge-Kutta
    !! and Adams-Bashforth coefficients for orders 1-4, then selects the
    !! appropriate method and allocates storage for old derivatives or stages.
    !!
    !! **Supported Methods:**
    !! - AB1, AB2, AB3, AB4: Adams-Bashforth (explicit multi-step)
    !! - RK1, RK2, RK3, RK4: Runge-Kutta (explicit multi-stage)
    !!
    !! **RK Coefficients (Butcher tableau):**
    !! - RK1: Forward Euler
    !! - RK2: Midpoint method
    !! - RK3: Strong Stability Preserving RK3 (SSP-RK3)
    !! - RK4: Classical fourth-order Runge-Kutta
    !!
    !! **AB Coefficients:**
    !! Derived from polynomial extrapolation of previous derivatives.
    implicit none

    type(time_intg_t) :: init                         !! Initialised time integrator
    class(base_backend_t), pointer :: backend         !! Computational backend
    class(allocator_t), pointer :: allocator          !! Memory allocator
    character(3), intent(in) :: method                !! Integration method ('AB3', 'RK4', etc.)
    integer, intent(in) :: nvars                      !! Number of variables to integrate

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
      init%step => adams_bashforth
    else if (init%sname(1:2) == 'RK') then
      read (init%sname(3:3), *, iostat=stat) init%order
      if (stat /= 0) error stop 'Error reading RK integration order'
      if (init%order >= 5) error stop 'Integration order >4 is not supported'
      init%nstep = 1
      init%nstage = init%order
      init%nolds = init%nstage
      init%step => runge_kutta
    else
      print *, 'Integration method '//init%sname//' is not defined'
      error stop
    end if

    init%nvars = nvars

    init%istep = 1
    init%istage = 1

    ! allocate memory
    allocate (init%olds(init%nvars, init%nolds))

    ! Request all the storage for old timesteps
    do i = 1, init%nvars
      do j = 1, init%nolds
        init%olds(i, j)%ptr => allocator%get_block(DIR_X)
      end do
    end do

  end function init

  subroutine runge_kutta(self, curr, deriv, dt)
    !! Advance solution using Runge-Kutta method.
    !!
    !! Implements explicit Runge-Kutta schemes of orders 1-4. The general
    !! form for an s-stage RK method is:
    !!
    !! \[ k_i = f(t_n + c_i \Delta t, u_n + \Delta t \sum_{j=1}^{i-1} a_{ij} k_j) \]
    !! \[ u_{n+1} = u_n + \Delta t \sum_{i=1}^{s} b_i k_i \]
    !!
    !! Where \( k_i \) are stage derivatives, \( a_{ij} \) are stage weights,
    !! and \( b_i \) are final combination weights. This implementation stores
    !! stage derivatives in `olds(:, 2:nstage+1)` and the initial solution in
    !! `olds(:, 1)`.
    !!
    !! The subroutine is called once per stage. When `istage == nstage`, it
    !! computes the final solution and resets the stage counter.
    implicit none

    class(time_intg_t), intent(inout) :: self  !! Time integrator state
    type(flist_t), intent(inout) :: curr(:)    !! Current solution (updated)
    type(flist_t), intent(in) :: deriv(:)      !! Stage derivative
    real(dp), intent(in) :: dt                 !! Timestep size

    integer :: i, j

    ! update solution
    if (self%istage == self%nstage) then
      do i = 1, self%nvars
        ! update step solution from stage derivative
        if (self%nstage > 1) then
          call self%backend%veccopy(curr(i)%ptr, self%olds(i, 1)%ptr)
        end if

        do j = 1, self%nstage - 1
          call self%backend%vecadd(self%rk_b(j, self%nstage)*dt, &
                                   self%olds(i, j + 1)%ptr, &
                                   1._dp, curr(i)%ptr)
        end do
        call self%backend%vecadd(self%rk_b(self%nstage, self%nstage)*dt, &
                                 deriv(i)%ptr, &
                                 1._dp, curr(i)%ptr)

      end do

      ! reset stage counter
      self%istage = 1
    else
      do i = 1, self%nvars
        ! save step initial condition
        if (self%istage == 1) then
          call self%backend%veccopy(self%olds(i, 1)%ptr, curr(i)%ptr)
        end if

        ! save stage derivative
        call self%backend%veccopy(self%olds(i, self%istage + 1)%ptr, &
                                  deriv(i)%ptr)

        ! update stage solution
        if (self%istage > 1) then
          call self%backend%veccopy(curr(i)%ptr, self%olds(i, 1)%ptr)
        end if
        do j = 1, self%istage
          call self%backend%vecadd(self%rk_a(j, self%istage, self%nstage)*dt, &
                                   self%olds(i, j + 1)%ptr, &
                                   1._dp, curr(i)%ptr)
        end do
      end do

      ! increment stage counter
      self%istage = self%istage + 1
    end if

  end subroutine runge_kutta

  subroutine adams_bashforth(self, curr, deriv, dt)
    !! Advance solution using Adams-Bashforth method.
    !!
    !! Implements explicit Adams-Bashforth schemes of orders 1-4. These
    !! multi-step methods use derivatives from previous timesteps:
    !!
    !! \[ u_{n+1} = u_n + \Delta t \sum_{i=0}^{s-1} b_i f_{n-i} \]
    !!
    !! Where \( f_{n-i} \) are stored derivatives from previous steps and
    !! \( b_i \) are the Adams-Bashforth coefficients. The method has an
    !! automatic startup phase: for the first `order` steps, it uses a
    !! lower-order scheme (e.g., AB2 uses AB1 on step 1, then AB2 on step 2+).
    !!
    !! Old derivatives are stored in `olds(:, 1:nstep-1)` and rotated after
    !! each step. The current derivative is used directly and then stored
    !! in `olds(:, 1)` for the next timestep.
    implicit none

    class(time_intg_t), intent(inout) :: self  !! Time integrator state
    type(flist_t), intent(inout) :: curr(:)    !! Current solution (updated)
    type(flist_t), intent(in) :: deriv(:)      !! Current time derivative
    real(dp), intent(in) :: dt                 !! Timestep size

    integer :: i, j
    integer :: nstep

    nstep = min(self%istep, self%nstep)
    do i = 1, self%nvars
      ! update solution
      call self%backend%vecadd(self%coeffs(1, nstep)*dt, &
                               deriv(i)%ptr, &
                               1._dp, curr(i)%ptr)
      do j = 2, nstep
        call self%backend%vecadd(self%coeffs(j, nstep)*dt, &
                                 self%olds(i, j - 1)%ptr, &
                                 1._dp, curr(i)%ptr)
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
        call self%backend%veccopy(self%olds(i, 1)%ptr, deriv(i)%ptr)
      end if
    end do

    ! increment step counter
    self%istep = self%istep + 1

  end subroutine adams_bashforth

  subroutine rotate(sol, n)
    !! Rotate pointer array for Adams-Bashforth old derivatives.
    !!
    !! Shifts pointers in the array to make room for a new derivative:
    !! sol(i) <- sol(i-1) for i from n down to 2, and sol(1) gets the
    !! old sol(n). This implements a circular buffer for old derivatives
    !! without copying data - only pointers are reassigned.
    !!
    !! Example for n=3: [new, old1, old2] becomes [?, new, old1]
    !! (where ? will be filled with the newest derivative)
    implicit none

    type(flist_t), intent(inout) :: sol(:)  !! Array of field list pointers to rotate
    integer, intent(in) :: n                !! Number of elements to rotate

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
