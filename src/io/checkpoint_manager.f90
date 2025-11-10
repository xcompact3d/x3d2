module m_checkpoint_manager
! @brief Manages the creation and restoration of simulation checkpoints
!! for restart capabilities.
!!
!! @details This module is responsible for periodically saving the full, unstrided
!! simulation state to a file. This allows a simulation to be stopped and resumed
!! from the exact state it was in.
!!
!! Key features include:
!! - Reading all checkpoint settings from a configuration file
!! - Periodically writing the full-resolution simulation state
!! - Handling the full logic for restarting a simulation from
!! a specified checkpoint file.
!! - A safe-write strategy that writes to a temporary file first,
!!   then atomically renames it to the final filename to
!! prevent corrupted checkpoints.
!! - Optional cleanup of old checkpoint files to conserve disk space.
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Abort
  use m_common, only: dp, i8, DIR_X, get_argument
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: reader_session_t, writer_session_t
  use m_config, only: checkpoint_config_t
  use m_io_field_utils, only: field_buffer_map_t, field_ptr_t, &
                              setup_field_arrays, cleanup_field_arrays, &
                              stride_data_to_buffer, get_output_dimensions, &
                              prepare_field_buffers, cleanup_field_buffers, &
                              write_single_field_to_buffer

  implicit none

  private
  public :: checkpoint_manager_t

  type :: checkpoint_manager_t
    type(checkpoint_config_t) :: config
    integer :: last_checkpoint_step = -1
    integer, dimension(3) :: full_resolution = [1, 1, 1]
    type(field_buffer_map_t), allocatable :: field_buffers(:)
    integer(i8), dimension(3) :: last_shape_dims = 0
    integer, dimension(3) :: last_stride_factors = 0
    integer(i8), dimension(3) :: last_output_shape = 0
  contains
    procedure :: init
    procedure :: handle_restart
    procedure :: handle_checkpoint_step
    procedure :: is_restart
    procedure :: finalise
    procedure, private :: write_checkpoint
    procedure, private :: restart_checkpoint
    procedure, private :: write_fields
    procedure, private :: cleanup_output_buffers
  end type checkpoint_manager_t

contains

  subroutine init(self, comm)
    !! Initialise checkpoint manager
    class(checkpoint_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    self%config = checkpoint_config_t()
    call self%config%read(nml_file=get_argument(1))

    if (self%config%checkpoint_freq > 0) then
      call configure_output(self, comm)
    end if
  end subroutine init

  subroutine configure_output(self, comm)
    !! Configure checkpoint output settings
    use m_io_backend, only: get_default_backend, IO_BACKEND_DUMMY
    class(checkpoint_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    integer :: myrank, ierr

    call MPI_Comm_rank(comm, myrank, ierr)

    if (myrank == 0 .and. get_default_backend() /= IO_BACKEND_DUMMY) then
      print *, 'Checkpoint frequency: ', self%config%checkpoint_freq
      print *, 'Keep all checkpoints: ', self%config%keep_checkpoint
      print *, 'Checkpoint prefix: ', trim(self%config%checkpoint_prefix)
    end if
  end subroutine configure_output

  function is_restart(self) result(restart)
    !! Check if this is a restart run
    class(checkpoint_manager_t), intent(in) :: self
    logical :: restart

    restart = self%config%restart_from_checkpoint
  end function is_restart

  subroutine handle_restart(self, solver, comm)
    !! Handle restart from checkpoint
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm

    character(len=256) :: restart_file
    integer :: restart_timestep
    real(dp) :: restart_time

    restart_file = trim(self%config%restart_file)
    if (solver%mesh%par%is_root()) then
      print *, 'Restarting from checkpoint: ', restart_file
    end if

    call self%restart_checkpoint(solver, restart_file, restart_timestep, &
                                 restart_time, comm)

    solver%current_iter = restart_timestep

    if (solver%mesh%par%is_root()) then
      print *, 'Successfully restarted from checkpoint at iteration ', &
        restart_timestep, ' with time ', restart_time
    end if
  end subroutine handle_restart

  subroutine handle_checkpoint_step(self, solver, timestep, comm)
    !! Handle checkpoint writing at a given timestep
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    integer :: comm_to_use

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm

    call self%write_checkpoint(solver, timestep, comm_to_use)
  end subroutine handle_checkpoint_step

  subroutine write_checkpoint(self, solver, timestep, comm)
    !! Write a checkpoint file for simulation restart
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in) :: comm

    character(len=256) :: filename, temp_filename, old_filename
    integer :: ierr, myrank
    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    real(dp) :: simulation_time
    logical :: file_exists
    type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
    integer :: data_loc, is_ab, i, j, n_total_vars
    type(writer_session_t) :: writer_session
    integer :: nolds_total, idx
    character(len=16), allocatable :: old_field_names(:)
    type(field_ptr_t), allocatable :: old_field_ptrs(:)

    if (self%config%checkpoint_freq <= 0) return
    if (mod(timestep, self%config%checkpoint_freq) /= 0) return

    call MPI_Comm_rank(comm, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%config%checkpoint_prefix), '_', timestep, '.bp'
    write (temp_filename, '(A,A)') &
      trim(self%config%checkpoint_prefix), '_temp.bp'

    call writer_session%open(temp_filename, comm)
    if (writer_session%is_session_functional() .and. myrank == 0) then
      print *, 'Writing checkpoint: ', trim(filename)
    end if

    simulation_time = timestep*solver%dt
    data_loc = solver%u%data_loc
    call writer_session%write_data("timestep", timestep)
    call writer_session%write_data("time", real(simulation_time, dp))
    call writer_session%write_data("dt", real(solver%dt, dp))
    call writer_session%write_data("data_loc", data_loc)

    n_total_vars = size(field_names)

    call setup_field_arrays(solver, field_names, field_ptrs, host_fields)

    call self%write_fields( &
      field_names, host_fields, &
      solver, writer_session, data_loc &
      )

    ! serialise time integrator metadata
    if (solver%time_integrator%sname(1:2) == 'AB') then
      is_ab = 1
    else
      is_ab = 0
    end if
    if (myrank == 0 .and. is_ab == 1) then
      print *, 'DEBUG CHECKPOINT: istep=', solver%time_integrator%istep, ' at iteration=', timestep
    end if
    call writer_session%write_data('ti_is_ab', is_ab)
    call writer_session%write_data('ti_order', solver%time_integrator%order)
    call writer_session%write_data('ti_istep', solver%time_integrator%istep)
    call writer_session%write_data('ti_nstep', solver%time_integrator%nstep)

    ! for AB methods with order >1, keep derivative history olds(i,j)
    if (is_ab == 1 .and. solver%time_integrator%order > 1) then
      ! Debug output before entering block
      if (myrank == 0) then
        print *, 'DEBUG: Writing', solver%time_integrator%nolds, 'olds for', n_total_vars, 'variables'
      end if
      
      nolds_total = solver%time_integrator%nolds*n_total_vars
      allocate (old_field_names(nolds_total))
      allocate (old_field_ptrs(nolds_total))
      idx = 0
      do i = 1, n_total_vars
        do j = 1, solver%time_integrator%nolds
          idx = idx + 1
          write (old_field_names(idx), '(A,"_rhs_old",I0)') &
            trim(field_names(i)), j
          old_field_ptrs(idx)%ptr => solver%time_integrator%olds(i, j)%ptr
        end do
      end do

      ! Olds fields are in DIR_X with blocked SIMD layout
      ! Use backend%get_field_data to extract to natural [nx,ny,nz] layout
      block
        integer :: phys_dims(3), rank_id, ierr_local
        integer(i8), dimension(3) :: phys_start, phys_count, phys_shape
        real(dp), allocatable :: temp_data_array(:,:,:,:)  ! One buffer per olds field
        
        call MPI_Comm_rank(comm, rank_id, ierr_local)
        
        ! Use unpadded cell dimensions for physical data
        phys_dims = solver%mesh%grid%cell_dims
        phys_start = int(solver%mesh%par%n_offset, i8)
        phys_count = int(phys_dims, i8)
        phys_shape = int(solver%mesh%grid%global_cell_dims, i8)
        
        ! Allocate array to hold ALL olds fields (cannot reuse buffer - ADIOS2 deferred writes!)
        allocate(temp_data_array(phys_dims(1), phys_dims(2), phys_dims(3), nolds_total))
        
        ! Extract all olds fields to natural layout
        do idx = 1, nolds_total
          call solver%backend%get_field_data(temp_data_array(:,:,:,idx), old_field_ptrs(idx)%ptr, DIR_X)
        end do
        
        ! Debug output for first field
        if (rank_id == 0) then
          print *, 'DEBUG WRITE: Saving', nolds_total, 'olds fields'
          do idx = 1, min(6, nolds_total)  ! Print first 6 fields
            print *, 'DEBUG WRITE:', trim(old_field_names(idx)), '(1,1,1)=', temp_data_array(1,1,1,idx), &
              ' min/max=', minval(temp_data_array(:,:,:,idx)), maxval(temp_data_array(:,:,:,idx))
          end do
        end if
        
        ! Now write all fields (data is stable in temp_data_array)
        do idx = 1, nolds_total
          call writer_session%write_data( &
            trim(old_field_names(idx)), &
            temp_data_array(:,:,:,idx), &
            start_dims=phys_start, count_dims=phys_count, &
            shape_dims=phys_shape &
            )
        end do
        
        ! DO NOT deallocate temp_data_array here - ADIOS2 needs it until close()
        call writer_session%close()
        deallocate(temp_data_array)
      end block
    else
      call writer_session%close()
    end if

    ! clean up buffers after session close (ADIOS2 deferred writes need them until end_step)
    call self%cleanup_output_buffers()
    if (allocated(old_field_ptrs)) then
      deallocate (old_field_ptrs)
    end if
    if (allocated(old_field_names)) then
      deallocate (old_field_names)
    end if

    call cleanup_field_arrays(solver, field_ptrs, host_fields)

    if (myrank == 0) then
      inquire (file=trim(temp_filename), exist=file_exists)
      if (file_exists) then
        ! Move temporary file to final checkpoint filename
        call execute_command_line('mv '//trim(temp_filename)//' '// &
                                  trim(filename))

        inquire (file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
          print *, 'ERROR: Checkpoint file not created: ', trim(filename)
        end if
      else
        ! temp file doesn't exist - skip file operations silently
      end if

      ! Remove old checkpoint if configured to keep only the latest
      if (.not. self%config%keep_checkpoint &
          .and. self%last_checkpoint_step > 0) then
        write (old_filename, '(A,A,I0.6,A)') &
          trim(self%config%checkpoint_prefix), '_', &
          self%last_checkpoint_step, '.bp'
        inquire (file=trim(old_filename), exist=file_exists)
        if (file_exists) then
          call execute_command_line('rm -rf '//trim(old_filename), &
                                    exitstat=ierr)
          if (ierr /= 0) then
            print *, 'WARNING: failed to remove old checkpoint: ', &
              trim(old_filename)
          end if
        end if
      end if
    end if

    self%last_checkpoint_step = timestep
  end subroutine write_checkpoint

  subroutine restart_checkpoint( &
    self, solver, filename, timestep, restart_time, comm &
    )
    !! Restart simulation state from checkpoint file
    class(checkpoint_manager_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    character(len=*), intent(in) :: filename
    integer, intent(out) :: timestep
    real(dp), intent(out) :: restart_time
    integer, intent(in) :: comm

    type(reader_session_t) :: reader_session
    integer :: ierr, myrank, data_loc
    integer :: dims(3)
    integer(i8), dimension(3) :: start_dims, count_dims
    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    logical :: file_exists
    integer :: ti_is_ab, ti_order, ti_istep, ti_nstep
    logical :: have_ti_meta
    character(len=16), allocatable :: var_names(:)
    integer :: n_total_vars
    character(len=64) :: old_name
    integer :: i, j

    call MPI_Comm_rank(comm, myrank, ierr)

    inquire (file=filename, exist=file_exists)
    if (.not. file_exists) then
      if (solver%mesh%par%is_root()) then
        print *, 'ERROR: Checkpoint file not found: ', trim(filename)
      end if
      call MPI_Abort(comm, 1, ierr)
      return
    end if

    call reader_session%open(filename, comm)
    
    if (myrank == 0) then
      print *, 'DEBUG: Checkpoint file opened successfully'
    end if
    
    call reader_session%read_data("timestep", timestep)
    call reader_session%read_data("time", restart_time)
    call reader_session%read_data("data_loc", data_loc)
    ! restore dt if present
    block
      real(dp) :: chk_dt
      call reader_session%read_data("dt", chk_dt)
      solver%dt = chk_dt
    end block
    ! attempt to read time integrator metadata
    have_ti_meta = .true.
    block
      call reader_session%read_data('ti_is_ab', ti_is_ab)
      call reader_session%read_data('ti_order', ti_order)
      call reader_session%read_data('ti_istep', ti_istep)
      call reader_session%read_data('ti_nstep', ti_nstep)
    end block

    dims = solver%mesh%get_dims(data_loc)
    start_dims = int(solver%mesh%par%n_offset, i8)
    count_dims = int(dims, i8)

    call solver%u%set_data_loc(data_loc)
    call solver%v%set_data_loc(data_loc)
    call solver%w%set_data_loc(data_loc)

    block
      real(dp), allocatable, target :: field_data_u(:, :, :)
      real(dp), allocatable, target :: field_data_v(:, :, :)
      real(dp), allocatable, target :: field_data_w(:, :, :)

      allocate (field_data_u(count_dims(1), count_dims(2), count_dims(3)))
      allocate (field_data_v(count_dims(1), count_dims(2), count_dims(3)))
      allocate (field_data_w(count_dims(1), count_dims(2), count_dims(3)))
      call reader_session%read_data("u", field_data_u, &
        start_dims=start_dims, count_dims=count_dims)
      call reader_session%read_data("v", field_data_v, &
        start_dims=start_dims, count_dims=count_dims)
      call reader_session%read_data("w", field_data_w, &
        start_dims=start_dims, count_dims=count_dims)
      call solver%backend%set_field_data(solver%u, field_data_u)
      call solver%backend%set_field_data(solver%v, field_data_v)
      call solver%backend%set_field_data(solver%w, field_data_w)
    end block

    ! restore AB derivative history if metadata indicates AB and order>1
    if (have_ti_meta) then
      if (ti_is_ab == 1) then
        if (solver%mesh%par%is_root()) then
          print *, 'DEBUG RESTORE: istep=', ti_istep, ' at iteration=', timestep
        end if
        solver%time_integrator%istep = ti_istep
        solver%time_integrator%nstep = ti_nstep
        if (solver%time_integrator%order /= ti_order) then
          if (solver%mesh%par%is_root()) then
            print *, 'WARNING: checkpoint AB order differs from current &
                  & solver config; using checkpoint order'
          end if
          solver%time_integrator%order = ti_order
        end if
        if (ti_order > 1) then
          ! Olds fields use blocked SIMD layout - use backend to restore
          block
            integer :: phys_dims(3)
            integer(i8), dimension(3) :: phys_start, phys_count
            real(dp), allocatable, target :: old_field(:, :, :)
            
            ! Use unpadded cell dimensions
            phys_dims = solver%mesh%grid%cell_dims
            phys_start = int(solver%mesh%par%n_offset, i8)
            phys_count = int(phys_dims, i8)
            
            allocate (old_field(phys_count(1), phys_count(2), phys_count(3)))
            
            n_total_vars = size(field_names)
            allocate (var_names(n_total_vars))
            var_names = field_names
            
            ! Zero out olds fields before restoring to ensure padding is initialized
            if (solver%mesh%par%is_root()) then
              print *, 'DEBUG: Zeroing olds fields before restore'
            end if
            do i = 1, n_total_vars
              do j = 1, solver%time_integrator%nolds
                call solver%time_integrator%olds(i, j)%ptr%fill(0.0_dp)
              end do
            end do
            
            ! Read old derivative data and restore using backend
            do i = 1, n_total_vars
              do j = 1, solver%time_integrator%nolds
                write (old_name, '(A,"_rhs_old",I0)') trim(var_names(i)), j
                
                call reader_session%read_data(trim(old_name), old_field, &
                  start_dims=phys_start, count_dims=phys_count)
                
                if (solver%mesh%par%is_root()) then
                  print *, 'DEBUG READ:', trim(old_name), '(1,1,1)=', old_field(1,1,1)
                end if
                
                ! Use backend to convert natural layout back to blocked layout
                call solver%backend%set_field_data( &
                  solver%time_integrator%olds(i, j)%ptr, old_field, DIR_X)
                
                ! Verify restoration with min/max
                if (solver%mesh%par%is_root()) then
                  print *, 'DEBUG: ', trim(old_name), ' restored min/max=', &
                    minval(solver%time_integrator%olds(i,j)%ptr%data), &
                    maxval(solver%time_integrator%olds(i,j)%ptr%data)
                end if
              end do
            end do
            deallocate (var_names)
            deallocate (old_field)
          end block
        end if
      end if
    end if

    call reader_session%close()
  end subroutine restart_checkpoint

  subroutine write_fields( &
    self, field_names, host_fields, solver, writer_session, data_loc &
    )
    !! Write field data for checkpoints (no striding)
    class(checkpoint_manager_t), intent(inout) :: self
    character(len=*), dimension(:), intent(in) :: field_names
    class(field_ptr_t), dimension(:), target, intent(in) :: host_fields
    class(solver_t), intent(in) :: solver
    type(writer_session_t), intent(inout) :: writer_session
    integer, intent(in) :: data_loc

    integer :: i_field
    integer(i8), dimension(3) :: output_start, output_count
    integer, dimension(3) :: output_dims_local

    ! Prepare buffers for full resolution (no striding for checkpoints)
    call prepare_field_buffers( &
      solver, self%full_resolution, field_names, data_loc, &
      self%field_buffers, self%last_shape_dims, self%last_stride_factors, &
      self%last_output_shape &
      )

    ! Calculate output dimensions for writing
    call get_output_dimensions( &
      int(solver%mesh%get_global_dims(data_loc), i8), &
      int(solver%mesh%par%n_offset, i8), &
      int(solver%mesh%get_dims(data_loc), i8), &
      self%full_resolution, &
      self%last_output_shape, output_start, output_count, &
      output_dims_local, &
      self%last_shape_dims, self%last_stride_factors, &
      self%last_output_shape &
      )

    do i_field = 1, size(field_names)
      call write_single_field_to_buffer( &
        trim(field_names(i_field)), host_fields(i_field)%ptr, &
        solver, self%full_resolution, data_loc, &
        self%field_buffers, self%last_shape_dims, self%last_stride_factors, &
        self%last_output_shape &
        )

      call writer_session%write_data( &
        trim(field_names(i_field)), &
        self%field_buffers(i_field)%buffer, &
        self%last_output_shape, &
        output_start, output_count &
        )
    end do
  end subroutine write_fields

  subroutine cleanup_output_buffers(self)
    !! Clean up dynamic field buffers
    class(checkpoint_manager_t), intent(inout) :: self

    call cleanup_field_buffers(self%field_buffers)
  end subroutine cleanup_output_buffers

  subroutine finalise(self)
    !! Clean up checkpoint manager
    class(checkpoint_manager_t), intent(inout) :: self

    call self%cleanup_output_buffers()
  end subroutine finalise

end module m_checkpoint_manager

