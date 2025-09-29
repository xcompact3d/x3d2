module m_io_manager
!! Simple facade using existing checkpoint_io module
  use m_checkpoint_manager, only: checkpoint_manager_t 
  use m_solver, only: solver_t

  implicit none

  private
  public :: io_manager_t

  type :: io_manager_t
    type(checkpoint_manager_t) :: impl
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
    call self%impl%init(comm)
  end subroutine io_init

  subroutine io_handle_restart(self, solver, comm)
    class(io_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm
    call self%impl%handle_restart(solver, comm)
  end subroutine io_handle_restart

  subroutine io_handle_step(self, solver, timestep, comm)
    class(io_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
    call self%impl%handle_io_step(solver, timestep, comm)
  end subroutine io_handle_step

  function io_is_restart(self) result(is_restart)
    class(io_manager_t), intent(in) :: self
    logical :: is_restart
    is_restart = self%impl%is_restart()
  end function io_is_restart

  subroutine io_finalise(self)
    class(io_manager_t), intent(inout) :: self
    call self%impl%finalise()
  end subroutine io_finalise

end module m_io_manager
