module m_base_case
  !! Provides the base case for running a simulation. New cases are
  !! implemented by extending this to specify the initial and boundary
  !! conditions, forcing terms and case-specific postprocessing and analysis.

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_C, VERT, MPI_X3D2_DP
  use m_monitoring, only: monitoring_t
  use m_field, only: field_t, flist_t
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, init
  use m_postprocess, only: compute_derived_fields, compute_pressure_vert
  use m_config, only: has_output_field
  use m_io_manager, only: io_manager_t
  use mpi, only: MPI_COMM_WORLD, MPI_Wtime, MPI_Reduce, MPI_MAX

  implicit none

  type, abstract :: base_case_t
    class(solver_t), allocatable :: solver
    type(io_manager_t) :: io_mgr
    type(monitoring_t) :: monitoring
    ! Persistent device BC fields (DIR_X, VERT), shared across all cases.
    ! Allocated by the derived case on first use, refilled per substep,
    ! released only at program end. A given case populates only the
    ! face/end pairs it needs; the rest stay null().
    class(field_t), pointer :: bc_start_u_x => null()
    class(field_t), pointer :: bc_start_v_x => null()
    class(field_t), pointer :: bc_start_w_x => null()
    class(field_t), pointer :: bc_end_u_x => null()
    class(field_t), pointer :: bc_end_v_x => null()
    class(field_t), pointer :: bc_end_w_x => null()
    class(field_t), pointer :: bc_start_u_y => null()
    class(field_t), pointer :: bc_start_v_y => null()
    class(field_t), pointer :: bc_start_w_y => null()
    class(field_t), pointer :: bc_end_u_y => null()
    class(field_t), pointer :: bc_end_v_y => null()
    class(field_t), pointer :: bc_end_w_y => null()

  contains
    procedure(define_BC), deferred :: define_BC
    procedure(initial_conditions), deferred :: initial_conditions
    procedure(forcings), deferred :: forcings
    procedure(apply_BC), deferred :: apply_BC
    procedure(postprocess), deferred :: postprocess
    procedure :: case_init
    procedure :: case_finalise
    procedure :: set_init
    procedure :: run
  end type base_case_t

  abstract interface
    subroutine define_BC(self)
      !! Applies case-specific boundary coinditions
      import :: base_case_t
      implicit none

      class(base_case_t) :: self
    end subroutine define_BC

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

    subroutine apply_BC(self, u, v, w)
      !! Applies case-specific pre-correction to the velocity fields before
      !! pressure correction
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self
      class(field_t), intent(inout) :: u, v, w
    end subroutine apply_BC

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

    call self%io_mgr%init(self%solver, MPI_COMM_WORLD)

    ! Tell the solver to persist pressure if output is enabled
    self%solver%keep_pressure = &
      has_output_field(self%io_mgr%snapshot_mgr%config, 'pressure') &
      .and. self%io_mgr%snapshot_mgr%config%snapshot_freq > 0

    if (self%io_mgr%is_restart()) then
      call self%io_mgr%handle_restart(self%solver, MPI_COMM_WORLD)
    else
      call self%initial_conditions()
    end if

    call self%monitoring%init(self%solver)

  end subroutine case_init

  subroutine case_finalise(self)
    class(base_case_t) :: self

    if (self%solver%mesh%par%is_root()) print *, 'run end'

    call self%monitoring%finalise()
    call self%io_mgr%finalise()
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

  subroutine run(self)
    !! Runs the solver forwards in time from t=t_0 to t=T, performing
    !! postprocessing/IO and reporting diagnostics.
    implicit none

    class(base_case_t), intent(inout) :: self

    type(flist_t), allocatable :: curr(:)
    type(flist_t), allocatable :: deriv(:)

    real(dp) :: t
    integer :: i, iter, sub_iter, start_iter
    logical :: output_vorticity, output_qcriterion
    ! t_start/t_end/t_step0/t_step1 hold raw MPI_Wtime readings, which
    ! are double precision; subtract in double, then convert
    ! it to the working precision.
    double precision :: t_start, t_end, t_step0, t_step1
    real(dp) :: t_total_local, t_total, t_step_local, t_step
    logical :: do_timing_report, is_root
    integer :: ierr

    ! Error: ‘t_step0’ may be used uninitialized in this function [-Werror=maybe-uninitialized]
    t_step0 = 0.0_dp

    output_vorticity = &
      has_output_field(self%io_mgr%snapshot_mgr%config, 'vorticity') &
      .and. self%io_mgr%snapshot_mgr%config%snapshot_freq > 0
    output_qcriterion = &
      has_output_field(self%io_mgr%snapshot_mgr%config, 'qcriterion') &
      .and. self%io_mgr%snapshot_mgr%config%snapshot_freq > 0

    if (self%io_mgr%is_restart()) then
      t = self%solver%current_iter*self%solver%dt
      if (self%solver%mesh%par%is_root()) &
        ! for restarts current_iter is read from the checkpoint file
        print *, 'Continuing from iteration:', &
        self%solver%current_iter, 'at time ', t
      if (self%solver%n_iters <= self%solver%current_iter) then
        error stop &
          'Restart requires n_iters greater than the restart iteration.'
      end if
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

    is_root = self%solver%mesh%par%is_root()
    t_start = MPI_Wtime()

    do iter = start_iter, self%solver%n_iters
      ! Report per-step timing on the final step and on every n_output-th
      ! step. This is computed identically on every rank so that the
      ! reduction below stays collective; the nested guard keeps mod()
      ! from being evaluated when n_output is 0.
      do_timing_report = (iter == self%solver%n_iters)
      if (self%solver%n_output > 0) then
        if (mod(iter, self%solver%n_output) == 0) do_timing_report = .true.
      end if

      ! Skip the first step: it carries one-off setup costs that would
      ! make the per-step figure unrepresentative.
      if (do_timing_report .and. iter > start_iter) then
        t_step0 = MPI_Wtime()
      end if
      do sub_iter = 1, self%solver%time_integrator%nstage
        ! first apply case-specific BCs
        call self%define_BC()

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

        call self%apply_BC(self%solver%u, self%solver%v, self%solver%w)
        if (self%solver%ibm_on) then
          call self%solver%ibm%body(self%solver%u, self%solver%v, &
                                    self%solver%w)
        end if

        call self%solver%pressure_correction(self%solver%u, self%solver%v, &
                                             self%solver%w)
      end do

      ! All ranks measure and reduce; only root prints. Gating the
      ! measurement on is_root would break the collective MPI_Reduce.
      if (do_timing_report .and. iter > start_iter) then
        t_step1 = MPI_Wtime()
        t_step_local = real(t_step1 - t_step0, 4)

        call MPI_Reduce(t_step_local, t_step, 1, MPI_X3D2_DP, MPI_MAX, &
                        0, MPI_COMM_WORLD, ierr)

        if (is_root) then
          print *, 'Time for this time step (s):', real(t_step, 4)
        end if
      end if

      self%solver%current_iter = iter

      ! Compute pressure on VERT grid if output_pressure is enabled
      if (self%solver%keep_pressure) then
        call compute_pressure_vert(self%solver)
      end if

      ! Compute postprocess fields (vorticity magnitude, Q-criterion)
      if (output_vorticity .or. output_qcriterion) then
        call compute_derived_fields(self%solver, &
                                    output_vorticity, &
                                    output_qcriterion)
      end if

      call self%io_mgr%update_stats(self%solver, iter)

      if (self%solver%n_output > 0) then
        if (mod(iter, self%solver%n_output) == 0) then
          t = iter*self%solver%dt

          call self%postprocess(iter, t)
        end if
      end if
      call self%io_mgr%handle_io_step(self%solver, &
                                      iter, MPI_COMM_WORLD)
    end do

    ! total + averaged-per-step over the whole run
    t_end = MPI_Wtime()
    t_total_local = real(t_end - t_start, 4)

    call MPI_Reduce(t_total_local, t_total, 1, MPI_X3D2_DP, MPI_MAX, &
                    0, MPI_COMM_WORLD, ierr)

    if (is_root) then
      print *, '==========================================================='
      print *, 'Averaged time per step (s):', &
        real(t_total/(self%solver%n_iters - (start_iter - 1)), 4)
      print *, 'Total wallclock (s):', real(t_total, 4)
      print *, 'Total wallclock (m):', real(t_total/60.0_dp, 4)
      print *, 'Total wallclock (h):', real(t_total/3600.0_dp, 4)
    end if
    call self%case_finalise

    ! deallocate memory
    deallocate (curr)
    deallocate (deriv)

  end subroutine run

end module m_base_case
