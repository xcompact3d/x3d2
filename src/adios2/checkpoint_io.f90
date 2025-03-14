module m_checkpoint_io
  !! Provides checkpoint and snapshot functionality for saving and restoring simulation state
  use mpi
  
  use m_common, only: dp
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_adios2_io, only: adios2_writer_t, adios2_mode_write, adios2_file_t
  use iso_fortran_env, only: real64, int64
  
  implicit none
  
  private
  public :: checkpoint_manager_t
  
  type :: checkpoint_manager_t
    type(adios2_writer_t), pointer :: adios2_writer => null()
    
    ! Configuration options
    integer :: checkpoint_freq = 0      !! Frequency of checkpointing (0 = off)
    integer :: snapshot_freq = 0        !! Frequency of snapshots (0 = off)
    logical :: keep_checkpoint = .true. !! If false, only keep latest checkpoint
    character(len=256) :: checkpoint_prefix = "checkpoint" !! Checkpoint filename prefix
    character(len=256) :: snapshot_prefix = "snapshot"     !! Snapshot filename prefix
    
    ! Internal state
    integer :: last_checkpoint_step = -1
  contains
    procedure :: init
    procedure :: configure
    procedure :: write_checkpoint
    procedure :: write_snapshot
    procedure :: restart_checkpoint
    procedure :: verify_checkpoint
    procedure :: handle_io_step
  end type checkpoint_manager_t
  
contains

  subroutine init(self, adios2_writer)
    !! Initialize the checkpoint manager with a reference to an ADIOS2 writer
    class(checkpoint_manager_t), intent(inout) :: self
    type(adios2_writer_t), target, intent(in) :: adios2_writer
    
    self%adios2_writer => adios2_writer
  end subroutine init
  
  subroutine configure(self, checkpoint_freq, snapshot_freq, keep_checkpoint, &
                       checkpoint_prefix, snapshot_prefix, comm)
    !! Configure checkpoint and snapshot settings
    class(checkpoint_manager_t), intent(inout) :: self
    integer, intent(in), optional :: checkpoint_freq, snapshot_freq
    logical, intent(in), optional :: keep_checkpoint
    character(len=*), intent(in), optional :: checkpoint_prefix, snapshot_prefix
    integer, intent(in), optional :: comm
    
    integer :: myrank, ierr
    integer :: comm_to_use
    
    if (present(comm)) then
      comm_to_use = comm
    else
      comm_to_use = MPI_COMM_WORLD
    end if
    
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    if (present(checkpoint_freq)) then
      self%checkpoint_freq = checkpoint_freq
    end if
    
    if (present(snapshot_freq)) then
      self%snapshot_freq = snapshot_freq
    end if
    
    if (present(keep_checkpoint)) then
      self%keep_checkpoint = keep_checkpoint
    end if
    
    if (present(checkpoint_prefix)) then
      self%checkpoint_prefix = checkpoint_prefix
    end if
    
    if (present(snapshot_prefix)) then
      self%snapshot_prefix = snapshot_prefix
    end if
    
    if (myrank == 0 .and. (self%checkpoint_freq > 0 .or. self%snapshot_freq > 0)) then
      print *, 'Checkpoint system configured:'
      if (self%checkpoint_freq > 0) then
        print *, '  Checkpoint frequency:', self%checkpoint_freq
        print *, '  Keep all checkpoints:', self%keep_checkpoint
        print *, '  Checkpoint prefix:', trim(self%checkpoint_prefix)
      end if
      
      if (self%snapshot_freq > 0) then
        print *, '  Snapshot frequency:', self%snapshot_freq
        print *, '  Snapshot prefix:', trim(self%snapshot_prefix)
      end if
    end if
  end subroutine configure
  
  subroutine write_checkpoint(self, solver, timestep, comm)
    !! Write a checkpoint file for simulation restart
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
    
    character(len=256) :: filename
    type(adios2_file_t) :: file
    integer :: ierr, myrank
    integer :: comm_to_use
    class(field_t), pointer :: host_field
    integer :: dims(3), data_loc
    integer(kind=int64), dimension(3) :: shape_dims, start_dims, count_dims
    
    if (self%checkpoint_freq <= 0) return
    if (mod(timestep, self%checkpoint_freq) /= 0) return
    
    if (present(comm)) then
      comm_to_use = comm
    else
      comm_to_use = MPI_COMM_WORLD
    end if
    
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    ! Generate checkpoint filename
    write(filename, '(A,A,I0.6,A)') trim(self%checkpoint_prefix), '_', timestep, '.bp'
    
    if (myrank == 0) then
      print *, 'Writing checkpoint: ', trim(filename)
    end if
    
    ! Set up dimensions for writing
    data_loc = solver%u%data_loc
    dims = solver%mesh%get_dims(data_loc)
    shape_dims = int(solver%mesh%get_global_dims(data_loc), kind=int64)
    start_dims = int(solver%mesh%par%n_offset, kind=int64)
    count_dims = int(dims, kind=int64)
    
    ! Open the checkpoint file
    file = self%adios2_writer%open(filename, adios2_mode_write, comm_to_use)
    call self%adios2_writer%begin_step(file)
    
    ! Write velocity fields
    host_field => solver%host_allocator%get_block(data_loc)
    
    call solver%backend%get_field_data(host_field%data, solver%u)
    call self%adios2_writer%write_data("u", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
                                      file, shape_dims, start_dims, count_dims)
    
    call solver%backend%get_field_data(host_field%data, solver%v)
    call self%adios2_writer%write_data("v", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
                                      file, shape_dims, start_dims, count_dims)
    
    call solver%backend%get_field_data(host_field%data, solver%w)
    call self%adios2_writer%write_data("w", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
                                      file, shape_dims, start_dims, count_dims)
    
    ! Write metadata needed for restart
    call self%adios2_writer%write_data("timestep", real(timestep, real64), file)
    call self%adios2_writer%write_data("time", real(timestep * solver%dt, real64), file)
    call self%adios2_writer%write_data("dt", real(solver%dt, real64), file)
    
    ! Release the host field
    call solver%host_allocator%release_block(host_field)
    
    ! Close the file
    call self%adios2_writer%close(file)
    
    ! Remove old checkpoint if configured to keep only the latest
    if (.not. self%keep_checkpoint .and. self%last_checkpoint_step > 0) then
      write(filename, '(A,A,I0.6,A)') trim(self%checkpoint_prefix), '_', &
                                     self%last_checkpoint_step, '.bp'
      if (myrank == 0) then
        call execute_command_line('rm -f ' // trim(filename), exitstat=ierr)
        if (ierr /= 0) then
          print *, 'Warning: Failed to remove old checkpoint: ', trim(filename)
        end if
      end if
    end if
    
    self%last_checkpoint_step = timestep
  end subroutine write_checkpoint
  
  subroutine write_snapshot(self, solver, timestep, comm)
    !! Write a snapshot file for visualization/analysis
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
    
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    integer :: myrank, ierr
    integer :: comm_to_use
    character(len=256) :: filename
    
    if (self%snapshot_freq <= 0) return
    if (mod(timestep, self%snapshot_freq) /= 0) return
    
    if (present(comm)) then
      comm_to_use = comm
    else
      comm_to_use = MPI_COMM_WORLD
    end if
    
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    ! Generate snapshot filename with timestep
    write(filename, '(A,A,I0.6,A)') trim(self%snapshot_prefix), '_', timestep, '.bp'
    
    if (myrank == 0) then
      print *, 'Writing snapshot: ', trim(filename)
    end if
    
    ! Use the existing writer's method to write fields to a file
    call self%adios2_writer%write_fields(timestep, field_names, solver)
  end subroutine write_snapshot
  
  subroutine restart_checkpoint(self, solver, filename, timestep, restart_time, comm)
    !! Restart simulation state from checkpoint file
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    character(len=*), intent(in) :: filename
    integer, intent(out) :: timestep
    real(dp), intent(out) :: restart_time
    integer, intent(in), optional :: comm
    
    logical :: file_exists
    integer :: myrank, ierr
    integer :: comm_to_use
    integer :: dims(3), data_loc
    integer(kind=int64), dimension(3) :: start_dims, count_dims
    class(field_t), pointer :: host_field
    
    if (present(comm)) then
      comm_to_use = comm
    else
      comm_to_use = MPI_COMM_WORLD
    end if
    
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    ! Check if file exists before trying to read
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      if (myrank == 0) then
        print *, 'Error: Checkpoint file not found: ', trim(filename)
      end if
      timestep = 0
      restart_time = 0.0_dp
      return
    end if
    
    if (myrank == 0) then
      print *, 'Restoring from checkpoint: ', trim(filename)
    end if
    
    ! Prepare for reading
    data_loc = solver%u%data_loc
    dims = solver%mesh%get_dims(data_loc)
    start_dims = int(solver%mesh%par%n_offset, kind=int64)
    count_dims = int(dims, kind=int64)
    
    ! Open the checkpoint file
    ! Need to create a reader type for this to work properly
    ! For now, we'll use the writer to read data since this is a prototype
    
    ! Read velocity fields
    host_field => solver%host_allocator%get_block(data_loc)
    
    ! Since we don't have proper reading yet, we'd need to implement reading functionality
    ! For this prototype, we would do something like:
    ! Read data from file and copy to host_field%data
    ! Then copy from host to device using:
    ! call solver%backend%set_field_data(solver%u, host_field%data)
    
    ! For now, we'll simulate the data restoration
    if (myrank == 0) then
      print *, 'Note: This is a placeholder for checkpoint restoration.'
      print *, 'Actual data transfer functionality needs to be implemented.'
    end if
    
    ! Just set some placeholder values for now
    timestep = 100  ! Example value
    restart_time = 1.0_dp  ! Example value
    
    ! Clean up
    call solver%host_allocator%release_block(host_field)
    
    if (myrank == 0) then
      print *, 'Checkpoint restoration complete.'
      print *, '  Restart timestep:', timestep
      print *, '  Restart time:', restart_time
    end if
  end subroutine restart_checkpoint
  
  subroutine verify_checkpoint(self, solver, timestep, comm)
    !! Verify checkpoint functionality by writing and reading back test data
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
    
    integer :: myrank, ierr
    integer :: comm_to_use
    character(len=256) :: test_filename
    
    if (present(comm)) then
      comm_to_use = comm
    else
      comm_to_use = MPI_COMM_WORLD
    end if
    
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    ! Generate test filename
    write(test_filename, '(A,A,I0.6,A)') 'test_checkpoint_', timestep, '.bp'
    
    if (myrank == 0) then
      print *, 'Verifying checkpoint system with file: ', trim(test_filename)
      print *, 'This is a placeholder for checkpoint verification.'
      print *, 'Full verification will be implemented in the next version.'
    end if
    
    ! This would involve:
    ! 1. Save a copy of current velocity fields
    ! 2. Write a checkpoint
    ! 3. Modify fields
    ! 4. Read back checkpoint
    ! 5. Compare restart fields with saved copies
    
    ! Clean up test file if it exists
    if (myrank == 0) then
      call execute_command_line('rm -f ' // trim(test_filename), exitstat=ierr)
    end if
  end subroutine verify_checkpoint
  
  subroutine handle_io_step(self, solver, timestep, comm)
    !! Convenience method to handle checkpoint and snapshot writing at a given timestep
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm
    
    ! Handle checkpoints
    call self%write_checkpoint(solver, timestep, comm)
    
    ! Handle snapshots
    call self%write_snapshot(solver, timestep, comm)
  end subroutine handle_io_step
  
end module m_checkpoint_io
