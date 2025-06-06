module m_checkpoint_manager_dummy
!! Null implementation that does nothing if ADIOS2 is not enabled
  use m_solver, only: solver_t
  use m_config, only: checkpoint_config_t

  implicit none

  private
  public :: checkpoint_manager_dummy_t

  type, abstract :: checkpoint_manager_dummy_t
    type(checkpoint_config_t) :: checkpoint_cfg
    logical :: is_restart = .false.
  contains
    procedure(init_interface), deferred :: init
    procedure(handle_restart_interface), deferred :: handle_restart
    procedure(handle_io_step_interface), deferred :: handle_io_step
    procedure(finalise_interface), deferred :: finalise
  end type checkpoint_manager_dummy_t

  abstract interface
    subroutine init_interface(self, comm)
      import :: checkpoint_manager_dummy_t
      class(checkpoint_manager_dummy_t), intent(inout) :: self
      integer, intent(in) :: comm
    end subroutine init_interface

    subroutine handle_restart_interface(self, solver, comm)
      import :: checkpoint_manager_dummy_t, solver_t
      class(checkpoint_manager_dummy_t), intent(inout) :: self
      class(solver_t), intent(inout) :: solver
      integer, intent(in), optional :: comm
    end subroutine handle_restart_interface

    subroutine handle_io_step_interface(self, solver, timestep, comm)
      import :: checkpoint_manager_dummy_t, solver_t
      class(checkpoint_manager_dummy_t), intent(inout) :: self
      class(solver_t), intent(in) :: solver
      integer, intent(in) :: timestep
      integer, intent(in), optional :: comm
    end subroutine handle_io_step_interface

    subroutine finalise_interface(self)
      import :: checkpoint_manager_dummy_t
      class(checkpoint_manager_dummy_t), intent(inout) :: self
    end subroutine finalise_interface
  end interface

end module m_checkpoint_manager_dummy

module m_checkpoint_manager_impl
  use mpi, only: MPI_Comm_rank, MPI_Barrier
  use m_common, only: get_argument
  use m_config, only: checkpoint_config_t
  use m_solver, only: solver_t
  use m_checkpoint_manager_dummy, only: checkpoint_manager_dummy_t
  use iso_fortran_env, only: output_unit

  implicit none

  private
  public :: checkpoint_manager_impl_t

  type, extends(checkpoint_manager_dummy_t) :: checkpoint_manager_impl_t
  contains
    procedure :: init
    procedure :: handle_restart
    procedure :: handle_io_step
    procedure :: finalise
  end type checkpoint_manager_impl_t

contains

  subroutine print_error(comm, message)
    integer, intent(in) :: comm
    character(len=*), intent(in) :: message
    integer :: rank, ierr

    call MPI_Comm_rank(comm, rank, ierr)
    if (rank == 0) then
      print *, "ERROR: ", trim(message)
      flush (output_unit)
    end if
    call MPI_Barrier(comm, ierr)
  end subroutine print_error

  subroutine init(self, comm)
    class(checkpoint_manager_impl_t), intent(inout) :: self
    integer, intent(in) :: comm
    integer :: ierr

    self%checkpoint_cfg = checkpoint_config_t()
    call self%checkpoint_cfg%read(nml_file=get_argument(1))

    ! print error if checkpointing is requested
    if (self%checkpoint_cfg%checkpoint_freq > 0 .or. &
        self%checkpoint_cfg%snapshot_freq > 0 .or. &
        self%checkpoint_cfg%restart_from_checkpoint) then
      call print_error( &
        comm, "Checkpoint functionality requested but ADIOS2 is not enabled. &
        & Recompile with -DENABLE_ADIOS2=ON to enable checkpointing")
      call MPI_Barrier(comm, ierr)
      stop 1
    end if
  end subroutine init

  subroutine handle_restart(self, solver, comm)
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm
  end subroutine handle_restart

  subroutine handle_io_step(self, solver, timestep, comm)
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
  end subroutine handle_io_step

  subroutine finalise(self)
    class(checkpoint_manager_impl_t), intent(inout) :: self
  end subroutine finalise

end module m_checkpoint_manager_impl

module m_checkpoint_manager
  use m_checkpoint_manager_impl, only: checkpoint_manager_impl_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: checkpoint_manager_t, create_checkpoint_manager

  type, public :: checkpoint_manager_t
    type(checkpoint_manager_impl_t) :: impl
  contains
    procedure :: init => cm_init
    procedure :: handle_restart => cm_handle_restart
    procedure :: handle_io_step => cm_handle_io_step
    procedure :: finalise => cm_finalise
    procedure :: is_restart => cm_is_restart
  end type checkpoint_manager_t

contains
  function create_checkpoint_manager(comm) result(mgr)
    integer, intent(in) :: comm
    type(checkpoint_manager_t) :: mgr

    call mgr%init(comm)
  end function create_checkpoint_manager

  subroutine cm_init(self, comm)
    class(checkpoint_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    call self%impl%init(comm)
  end subroutine cm_init

  subroutine cm_handle_restart(self, solver, comm)
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm

    call self%impl%handle_restart(solver, comm)
  end subroutine cm_handle_restart

  subroutine cm_handle_io_step(self, solver, timestep, comm)
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    call self%impl%handle_io_step(solver, timestep, comm)
  end subroutine cm_handle_io_step

  subroutine cm_finalise(self)
    class(checkpoint_manager_t), intent(inout) :: self

    call self%impl%finalise()
  end subroutine cm_finalise

  function cm_is_restart(self) result(is_restart)
    class(checkpoint_manager_t), intent(in) :: self
    logical :: is_restart

    is_restart = self%impl%checkpoint_cfg%restart_from_checkpoint
  end function cm_is_restart

end module m_checkpoint_manager
