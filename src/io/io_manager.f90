module m_io_manager
  !! High-level manager orchestrating checkpoint and snapshot operations.
  !!
  !! This module acts as a facade to the I/O subsystem, simplifying the main
  !! simulation loop by providing a single point of contact for all I/O-related
  !! actions. The main program only needs to interact with `io_manager_t`, which
  !! delegates tasks to specialised checkpoint and snapshot managers.
  !!
  !! **Responsibilities:**
  !!
  !! - Initialise checkpoint and snapshot managers
  !! - Coordinate restart from checkpoints
  !! - Orchestrate periodic checkpoint and snapshot writes
  !! - Finalise I/O operations and clean up resources
  !!
  !! **Usage Pattern:**
  !!
  !! ```fortran
  !! type(io_manager_t) :: io_mgr
  !! call io_mgr%init(comm)
  !! if (io_mgr%is_restart()) call io_mgr%handle_restart(solver, comm)
  !! do timestep = 1, n_steps
  !!   call io_mgr%handle_io_step(solver, timestep, comm)
  !! end do
  !! call io_mgr%finalise()
  !! ```
  use m_checkpoint_manager, only: checkpoint_manager_t
  use m_snapshot_manager, only: snapshot_manager_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: io_manager_t

  type :: io_manager_t
    !! Unified manager for checkpoint and snapshot operations.
    !!
    !! Contains both checkpoint and snapshot managers and provides
    !! a simplified interface for the main simulation loop.
    type(checkpoint_manager_t) :: checkpoint_mgr !! Manages restart and checkpoint files
    type(snapshot_manager_t) :: snapshot_mgr     !! Manages visualisation output files
  contains
    procedure :: init => io_init                   !! Initialise I/O managers
    procedure :: handle_restart => io_handle_restart !! Load restart data if needed
    procedure :: handle_io_step => io_handle_step  !! Process checkpoints/snapshots for timestep
    procedure :: finalise => io_finalise           !! Finalise and clean up
    procedure :: is_restart => io_is_restart       !! Check if simulation is restarting
  end type io_manager_t

contains

  subroutine io_init(self, comm)
    !! Initialise checkpoint and snapshot managers.
    !!
    !! Sets up both managers by passing the MPI communicator. Each manager
    !! reads its configuration and prepares for I/O operations.
    implicit none

    class(io_manager_t), intent(inout) :: self !! I/O manager instance
    integer, intent(in) :: comm                 !! MPI communicator

    call self%checkpoint_mgr%init(comm)
    call self%snapshot_mgr%init(comm)
  end subroutine io_init

  subroutine io_handle_restart(self, solver, comm)
    !! Handle restart by loading checkpoint data.
    !!
    !! Delegates to the checkpoint manager to load solver state from
    !! the most recent checkpoint file. Should only be called if
    !! `is_restart()` returns true.
    implicit none

    class(io_manager_t), intent(inout) :: self    !! I/O manager instance
    class(solver_t), intent(inout) :: solver      !! Solver to load state into
    integer, intent(in), optional :: comm         !! MPI communicator (optional)

    call self%checkpoint_mgr%handle_restart(solver, comm)
  end subroutine io_handle_restart

  subroutine io_handle_step(self, solver, timestep, comm)
    !! Handle I/O operations for current timestep.
    !!
    !! Checks if checkpoint or snapshot output is required at this timestep
    !! and writes data accordingly. Typically called at the end of each
    !! timestep in the main simulation loop.
    implicit none

    class(io_manager_t), intent(inout) :: self !! I/O manager instance
    class(solver_t), intent(in) :: solver      !! Solver containing current state
    integer, intent(in) :: timestep             !! Current timestep number
    integer, intent(in), optional :: comm       !! MPI communicator (optional)

    call self%checkpoint_mgr%handle_checkpoint_step(solver, timestep, comm)
    call self%snapshot_mgr%handle_snapshot_step(solver, timestep, comm)
  end subroutine io_handle_step

  function io_is_restart(self) result(is_restart)
    !! Check if simulation is restarting from checkpoint.
    !!
    !! Queries the checkpoint manager to determine if a restart file
    !! exists and should be loaded.
    implicit none

    class(io_manager_t), intent(in) :: self !! I/O manager instance
    logical :: is_restart                    !! True if restarting from checkpoint

    is_restart = self%checkpoint_mgr%is_restart()
  end function io_is_restart

  subroutine io_finalise(self)
    !! Finalise I/O operations and clean up resources.
    !!
    !! Closes any open files and releases resources held by both
    !! checkpoint and snapshot managers. Should be called at the end
    !! of the simulation.
    implicit none

    class(io_manager_t), intent(inout) :: self !! I/O manager instance

    call self%checkpoint_mgr%finalise()
    call self%snapshot_mgr%finalise()
  end subroutine io_finalise

end module m_io_manager
