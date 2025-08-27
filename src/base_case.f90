module m_base_case
  !! Provides the base case for running a simulation. New cases are
  !! implemented by extending this to specify the initial and boundary
  !! conditions, forcing terms and case-specific postprocessing and analysis.

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Z, DIR_C, VERT
  use m_field, only: field_t, flist_t
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, init
  use m_checkpoint_manager, only: checkpoint_manager_t
  use mpi, only: MPI_COMM_WORLD

  implicit none

  type, abstract :: base_case_t
    class(solver_t), allocatable :: solver
    type(checkpoint_manager_t) :: checkpoint_mgr
  contains
    procedure(boundary_conditions), deferred :: boundary_conditions
    procedure(initial_conditions), deferred :: initial_conditions
    procedure(forcings), deferred :: forcings
    procedure(pre_correction), deferred :: pre_correction
    procedure(postprocess), deferred :: postprocess
    procedure :: case_init
    procedure :: case_finalise
    procedure :: set_init
    procedure :: run
    procedure :: print_enstrophy
    procedure :: print_div_max_mean
  end type base_case_t

  abstract interface
    subroutine boundary_conditions(self)
      !! Applies case-specific boundary coinditions
      import :: base_case_t
      implicit none

      class(base_case_t) :: self
    end subroutine boundary_conditions

    subroutine initial_conditions(self)
      !! Sets case-specific initial conditions
      import :: base_case_t
      implicit none

      class(base_case_t) :: self
    end subroutine initial_conditions

    subroutine forcings(self, du, dv, dw, iter)
      !! Applies case-specific or model realated forcings after transeq
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      integer, intent(in) :: iter
    end subroutine forcings

    subroutine pre_correction(self, u, v, w)
      !! Applies case-specific pre-correction to the velocity fields before
      !! pressure correction
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self
      class(field_t), intent(inout) :: u, v, w
    end subroutine pre_correction

    subroutine postprocess(self, iter, t)
      !! Triggers case-specific postprocessings at user specified intervals
      import :: base_case_t
      import :: dp
      implicit none

      class(base_case_t) :: self
      integer, intent(in) :: iter
      real(dp), intent(in) :: t
    end subroutine postprocess
  end interface

contains

  subroutine case_init(self, backend, mesh, host_allocator)
    implicit none

    class(base_case_t) :: self
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator

    self%solver = init(backend, mesh, host_allocator)

    call self%checkpoint_mgr%init(MPI_COMM_WORLD)
    if (self%checkpoint_mgr%is_restart()) then
      call self%checkpoint_mgr%handle_restart(self%solver, MPI_COMM_WORLD)
    else
      call self%initial_conditions()
    end if

  end subroutine case_init

  subroutine case_finalise(self)
    class(base_case_t) :: self

    if (self%solver%mesh%par%is_root()) print *, 'run end'

    call self%checkpoint_mgr%finalise()
  end subroutine case_finalise

  subroutine set_init(self, field, field_func)
    implicit none

    class(base_case_t) :: self
    class(field_t), intent(inout) :: field

    interface
      pure function field_func(coords) result(r)
        import dp
        implicit none

        real(dp), intent(in) :: coords(3)
        real(dp) :: r
      end function field_func
    end interface

    class(field_t), pointer :: field_init

    integer :: i, j, k, dims(3)
    real(dp) :: coords(3)

    dims = self%solver%mesh%get_dims(VERT)
    field_init => self%solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%solver%mesh%get_coordinates(i, j, k)
          field_init%data(i, j, k) = field_func(coords)
        end do
      end do
    end do

    call self%solver%backend%set_field_data(field, field_init%data)

    call self%solver%host_allocator%release_block(field_init)

  end subroutine set_init

  subroutine print_enstrophy(self, u, v, w)
    !! Reports the enstrophy
    implicit none

    class(base_case_t), intent(in) :: self
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: du, dv, dw
    real(dp) :: enstrophy

    du => self%solver%backend%allocator%get_block(DIR_X, VERT)
    dv => self%solver%backend%allocator%get_block(DIR_X, VERT)
    dw => self%solver%backend%allocator%get_block(DIR_X, VERT)

    call self%solver%curl(du, dv, dw, u, v, w)

    enstrophy = 0.5_dp*(self%solver%backend%scalar_product(du, du) &
                        + self%solver%backend%scalar_product(dv, dv) &
                        + self%solver%backend%scalar_product(dw, dw)) &
                /self%solver%ngrid

    if (self%solver%mesh%par%is_root()) print *, 'enstrophy:', enstrophy

    call self%solver%backend%allocator%release_block(du)
    call self%solver%backend%allocator%release_block(dv)
    call self%solver%backend%allocator%release_block(dw)

  end subroutine print_enstrophy

  subroutine print_div_max_mean(self, u, v, w)
    !! Reports the div(u) at cell centres
    implicit none

    class(base_case_t), intent(in) :: self
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: div_u
    real(dp) :: div_u_max, div_u_mean

    div_u => self%solver%backend%allocator%get_block(DIR_Z)

    call self%solver%divergence_v2p(div_u, u, v, w)

    call self%solver%backend%field_max_mean(div_u_max, div_u_mean, div_u)
    if (self%solver%mesh%par%is_root()) &
      print *, 'div u max mean:', div_u_max, div_u_mean

    call self%solver%backend%allocator%release_block(div_u)

  end subroutine print_div_max_mean

  subroutine run(self)
    !! Runs the solver forwards in time from t=t_0 to t=T, performing
    !! postprocessing/IO and reporting diagnostics.
    implicit none

    class(base_case_t), intent(inout) :: self

    type(flist_t), allocatable :: curr(:)
    type(flist_t), allocatable :: deriv(:)

    real(dp) :: t
    integer :: i, iter, sub_iter, start_iter

    if (self%checkpoint_mgr%is_restart()) then
      t = self%solver%current_iter*self%solver%dt
      if (self%solver%mesh%par%is_root()) &
        ! for restarts current_iter is read from the checkpoint file
        print *, 'Continuing from iteration:', &
        self%solver%current_iter, 'at time ', t
    else
      self%solver%current_iter = 0
      if (self%solver%mesh%par%is_root()) print *, 'initial conditions'
      t = 0._dp
    end if

    call self%postprocess(self%solver%current_iter, t)
    start_iter = self%solver%current_iter + 1

    if (self%solver%mesh%par%is_root()) print *, 'start run'

    allocate (curr(self%solver%time_integrator%nvars))
    allocate (deriv(self%solver%time_integrator%nvars))

    curr(1)%ptr => self%solver%u
    curr(2)%ptr => self%solver%v
    curr(3)%ptr => self%solver%w
    do i = 1, self%solver%nspecies
      curr(3 + i)%ptr => self%solver%species(i)%ptr
    end do

    do iter = start_iter, self%solver%n_iters
      do sub_iter = 1, self%solver%time_integrator%nstage
        ! first apply case-specific BCs
        call self%boundary_conditions()

        do i = 1, self%solver%nvars
          deriv(i)%ptr => self%solver%backend%allocator%get_block(DIR_X)
        end do

        call self%solver%transeq(deriv, curr)

        ! models that introduce source terms handled here
        call self%forcings(deriv(1)%ptr, deriv(2)%ptr, deriv(3)%ptr, iter)

        ! time integration
        call self%solver%time_integrator%step(curr, deriv, self%solver%dt)

        do i = 1, self%solver%nvars
          call self%solver%backend%allocator%release_block(deriv(i)%ptr)
        end do

        call self%pre_correction(self%solver%u, self%solver%v, self%solver%w)
        if (self%solver%ibm_on) then
          call self%solver%ibm%body(self%solver%u, self%solver%v, &
                                    self%solver%w)
        end if

        call self%solver%pressure_correction(self%solver%u, self%solver%v, &
                                             self%solver%w)
      end do

      self%solver%current_iter = iter

      if (mod(iter, self%solver%n_output) == 0) then
        t = iter*self%solver%dt

        call self%postprocess(iter, t)
      end if

      call self%checkpoint_mgr%handle_io_step(self%solver, &
                                              iter, MPI_COMM_WORLD)
    end do

    call self%case_finalise

    ! deallocate memory
    deallocate (curr)
    deallocate (deriv)

  end subroutine run

end module m_base_case
