module m_base_case
  !! Base class for flow simulation cases.
  !!
  !! This abstract base class provides the framework for implementing specific
  !! flow cases (channel, TGV, generic, etc.). New cases extend this class and
  !! override deferred procedures to specify:
  !! - **Initial conditions**: Set velocity and other field initial states
  !! - **Boundary conditions**: Apply physical boundary conditions each timestep
  !! - **Forcing terms**: Add body forces or model-specific source terms
  !! - **Pre-correction**: Modify velocity before pressure correction (e.g., IBM)
  !! - **Postprocessing**: Compute statistics, output diagnostics, etc.
  !!
  !! **Simulation Workflow:**
  !! The `run()` method orchestrates the time integration loop:
  !! 1. Apply boundary conditions
  !! 2. Advance solution one timestep via solver%step()
  !! 3. Write checkpoints/snapshots (via checkpoint_mgr)
  !! 4. Perform case-specific postprocessing
  !! 5. Repeat until final time reached
  !!
  !! **Time Integration:**
  !! Each timestep involves multiple stages (for RK) or steps (for AB):
  !! - Transport equation (transeq) computes velocity derivatives
  !! - Forcing terms applied after transeq
  !! - Pre-correction modifies velocity (e.g., for immersed boundaries)
  !! - Pressure correction enforces incompressibility
  !!
  !! **Restart Capability:**
  !! The checkpoint manager handles restart from saved states automatically
  !! if a restart file is detected.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Z, DIR_C, VERT
  use m_field, only: field_t, flist_t
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, init
  use m_io_manager, only: io_manager_t
  use mpi, only: MPI_COMM_WORLD

  implicit none

  type, abstract :: base_case_t
    !! Abstract base type for flow cases.
    !!
    !! Derived types must implement all deferred procedures to define
    !! case-specific behaviour.
    class(solver_t), allocatable :: solver       !! Incompressible Navier-Stokes solver
    type(io_manager_t) :: checkpoint_mgr         !! Checkpoint and snapshot manager
  contains
    procedure(boundary_conditions), deferred :: boundary_conditions !! Apply BCs (deferred)
    procedure(initial_conditions), deferred :: initial_conditions   !! Set ICs (deferred)
    procedure(forcings), deferred :: forcings                       !! Add forcing terms (deferred)
    procedure(pre_correction), deferred :: pre_correction           !! Pre-pressure correction (deferred)
    procedure(postprocess), deferred :: postprocess                 !! Case-specific analysis (deferred)
    procedure :: case_init          !! Initialise case and solver
    procedure :: case_finalise      !! Clean up and finalise
    procedure :: set_init           !! Set initial conditions and prepare for run
    procedure :: run                !! Main time integration loop
    procedure :: print_enstrophy    !! Print enstrophy diagnostic
    procedure :: print_div_max_mean !! Print divergence diagnostics
  end type base_case_t

  abstract interface
    subroutine boundary_conditions(self)
      !! Abstract interface for applying boundary conditions.
      !!
      !! Called each timestep before computing derivatives. Implementations
      !! should set velocity and scalar values at domain boundaries according
      !! to the physical boundary conditions (Dirichlet, Neumann, periodic, etc.).
      import :: base_case_t
      implicit none

      class(base_case_t) :: self !! Case instance
    end subroutine boundary_conditions

    subroutine initial_conditions(self)
      !! Abstract interface for setting initial conditions.
      !!
      !! Called once during initialisation to set the initial state of velocity
      !! and scalar fields. Implementations should populate u, v, w (and species
      !! if present) with case-appropriate initial values.
      import :: base_case_t
      implicit none

      class(base_case_t) :: self !! Case instance
    end subroutine initial_conditions

    subroutine forcings(self, du, dv, dw, iter)
      !! Abstract interface for applying forcing terms.
      !!
      !! Called after transport equation (transeq) but before pressure correction.
      !! Add body forces, source terms, or model-specific forcings (e.g., mean
      !! pressure gradient for channel flow, immersed boundary forces, etc.).
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self                  !! Case instance
      class(field_t), intent(inout) :: du, dv, dw !! Velocity derivatives to modify
      integer, intent(in) :: iter                 !! Current iteration number
    end subroutine forcings

    subroutine pre_correction(self, u, v, w)
      !! Abstract interface for pre-pressure correction modifications.
      !!
      !! Called after forcings but before pressure correction. Used for operations
      !! that need to modify the velocity field before enforcing incompressibility,
      !! such as immersed boundary method (IBM) velocity corrections.
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self                 !! Case instance
      class(field_t), intent(inout) :: u, v, w   !! Velocity fields to modify
    end subroutine pre_correction

    subroutine postprocess(self, iter, t)
      !! Abstract interface for case-specific postprocessing.
      !!
      !! Called at user-specified intervals during time integration. Implement
      !! this to compute statistics, output diagnostics, write custom data files,
      !! or perform any case-specific analysis.
      import :: base_case_t
      import :: dp
      implicit none

      class(base_case_t) :: self    !! Case instance
      integer, intent(in) :: iter   !! Current iteration number
      real(dp), intent(in) :: t     !! Current simulation time
    end subroutine postprocess
  end interface

contains

  subroutine case_init(self, backend, mesh, host_allocator)
    !! Initialise case with solver and checkpoint manager.
    !!
    !! Creates the solver instance and initialises the checkpoint/snapshot
    !! manager. If a restart file is detected, loads the saved state.
    implicit none

    class(base_case_t) :: self                              !! Case instance
    class(base_backend_t), target, intent(inout) :: backend !! Computational backend
    type(mesh_t), target, intent(inout) :: mesh             !! Mesh with decomposition
    type(allocator_t), target, intent(inout) :: host_allocator !! Host memory allocator

    self%solver = init(backend, mesh, host_allocator)

    call self%checkpoint_mgr%init(MPI_COMM_WORLD)
    if (self%checkpoint_mgr%is_restart()) then
      call self%checkpoint_mgr%handle_restart(self%solver, MPI_COMM_WORLD)
    else
      call self%initial_conditions()
    end if

  end subroutine case_init

  subroutine case_finalise(self)
    !! Finalise the case and clean up resources.
    !!
    !! Performs cleanup operations at the end of a simulation run:
    !! - Finalises the checkpoint manager (closes files, flushes buffers)
    !! - Prints completion message on root process
    !!
    !! This should be called after the main time integration loop completes.
    class(base_case_t) :: self !! Case instance to finalise

    if (self%solver%mesh%par%is_root()) print *, 'run end'

    call self%checkpoint_mgr%finalise()
  end subroutine case_finalise

  subroutine set_init(self, field, field_func)
    !! Initialise a field using an analytical function.
    !!
    !! This utility subroutine sets a field's values by evaluating a
    !! user-provided pure function at each grid point. The function
    !! is evaluated on the host, then transferred to the backend device
    !! (if using GPU backend).
    !!
    !! **Usage Example:**
    !! ```fortran
    !! call self%set_init(self%solver%u, u_initial)
    !! ```
    !! where `u_initial` is a pure function taking coordinates [x,y,z]
    !! and returning the initial velocity value.
    !!
    !! This is commonly used in `initial_conditions()` implementations
    !! to set velocity or scalar fields from analytical expressions.
    implicit none

    class(base_case_t) :: self                  !! Case instance
    class(field_t), intent(inout) :: field      !! Field to initialise

    interface
      pure function field_func(coords) result(r)
        !! Pure function defining field values at each point.
        import dp
        implicit none

        real(dp), intent(in) :: coords(3)  !! Spatial coordinates [x, y, z]
        real(dp) :: r                      !! Field value at this location
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
    !! Compute and print the volume-averaged enstrophy.
    !!
    !! Enstrophy is a measure of the rotational kinetic energy density:
    !! \[ E = \frac{1}{2V} \int_V |\omega|^2 \, dV = \frac{1}{2V} \int_V |\nabla \times \mathbf{u}|^2 \, dV \]
    !!
    !! where \( \omega = \nabla \times \mathbf{u} \) is the vorticity.
    !!
    !! This diagnostic is useful for monitoring:
    !! - Flow transition to turbulence (enstrophy increases)
    !! - Energy cascade to small scales
    !! - Numerical stability (sudden spikes indicate problems)
    !! - Comparison with theoretical predictions (e.g., TGV decay)
    !!
    !! Only the root MPI rank prints the result.
    implicit none

    class(base_case_t), intent(in) :: self  !! Case instance
    class(field_t), intent(in) :: u, v, w   !! Velocity components

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
    !! Compute and print maximum and mean divergence.
    !!
    !! For incompressible flow, the velocity divergence should be zero:
    !! \[ \nabla \cdot \mathbf{u} = 0 \]
    !!
    !! This diagnostic reports:
    !! - **Maximum divergence**: Largest local violation of incompressibility
    !! - **Mean divergence**: Volume-averaged divergence (should be near machine zero)
    !!
    !! **Purpose:**
    !! - Monitor quality of pressure correction (divergence should be ~ 1e-10 or smaller)
    !! - Detect numerical issues (large divergence indicates solver problems)
    !! - Verify proper boundary condition implementation
    !! - Check convergence of iterative Poisson solvers
    !!
    !! Divergence is computed at cell centres from vertex velocities using
    !! staggered derivatives and interpolation.
    !!
    !! Only the root MPI rank prints the result.
    implicit none

    class(base_case_t), intent(in) :: self  !! Case instance
    class(field_t), intent(in) :: u, v, w   !! Velocity components

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
    !! Main time integration loop for the simulation.
    !!
    !! Advances the solution from initial time t=t_0 to final time t=T,
    !! orchestrating all aspects of the simulation:
    !!
    !! **Each Timestep:**
    !! 1. Apply boundary conditions
    !! 2. Compute derivatives and advance via time_integrator%step()
    !! 3. Handle checkpointing and snapshot output (via checkpoint_mgr)
    !! 4. Perform case-specific postprocessing
    !! 5. Print diagnostics (divergence, enstrophy)
    !!
    !! **Time Integration Stages:**
    !! For multi-stage methods (RK), each timestep involves multiple stages.
    !! The solver%step() method handles the stage-by-stage advancement,
    !! calling transeq, forcings, pre_correction, and pressure_correction
    !! at appropriate points.
    !!
    !! **Restart Support:**
    !! If a restart file is detected, continues from the saved iteration
    !! and time rather than starting from t=0.
    implicit none

    class(base_case_t), intent(inout) :: self !! Case instance

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
