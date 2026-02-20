module m_io_manager
!! @brief Provides a high-level manager that orchestrates all checkpoint and
!! snapshot operations.
!!
!! @details This module acts as a facade to the I/O subsystem.
!! Its purpose is to simplify the main simulation loop by providing
!! a single point of contact for all I/O-related actions. The mainprogram only
!! needs to interact with the `io_manager_t` type, which then delegates tasks
!! to the specialised checkpoint and snapshot managers.
  use m_checkpoint_manager, only: checkpoint_manager_t
  use m_snapshot_manager, only: snapshot_manager_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: io_manager_t

  type :: io_manager_t
    type(checkpoint_manager_t) :: checkpoint_mgr
    type(snapshot_manager_t) :: snapshot_mgr
  contains
    procedure :: init => io_init
    procedure :: handle_restart => io_handle_restart
    procedure :: handle_io_step => io_handle_step
    procedure :: finalise => io_finalise
    procedure :: is_restart => io_is_restart
  end type io_manager_t

contains

  subroutine io_init(self, comm)
    class(io_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    call self%checkpoint_mgr%init(comm)
    call self%snapshot_mgr%init(comm)
  end subroutine io_init

  subroutine io_handle_restart(self, solver, comm)
    class(io_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm

    call self%checkpoint_mgr%handle_restart(solver, comm)
  end subroutine io_handle_restart

  subroutine io_handle_step(self, solver, timestep, comm)
    class(io_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    call self%checkpoint_mgr%handle_checkpoint_step(solver, timestep, comm)
    call self%snapshot_mgr%handle_snapshot_step(solver, timestep, comm)
  end subroutine io_handle_step

  function io_is_restart(self) result(is_restart)
    class(io_manager_t), intent(in) :: self
    logical :: is_restart

    is_restart = self%checkpoint_mgr%is_restart()
  end function io_is_restart

  subroutine io_finalise(self)
    class(io_manager_t), intent(inout) :: self

    call self%checkpoint_mgr%finalise()
    call self%snapshot_mgr%finalise()
  end subroutine io_finalise

end module m_io_manager
