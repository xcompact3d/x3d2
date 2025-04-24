module m_checkpoint_manager_base
  use m_solver, only: solver_t
  use m_config, only: checkpoint_config_t
  
  implicit none
  
  private
  public :: checkpoint_manager_base_t
  
  type, abstract :: checkpoint_manager_base_t
    type(checkpoint_config_t) :: checkpoint_cfg
    logical :: is_restart = .false.
  contains
    procedure(init_interface), deferred :: init
    procedure(handle_restart_interface), deferred :: handle_restart
    procedure(handle_io_step_interface), deferred :: handle_io_step
    procedure(finalise_interface), deferred :: finalise
  end type checkpoint_manager_base_t

  abstract interface
    subroutine init_interface(self, comm)
      import :: checkpoint_manager_base_t
      class(checkpoint_manager_base_t), intent(inout) :: self
      integer, intent(in) :: comm
    end subroutine init_interface

    subroutine handle_restart_interface(self, solver, comm)
      import :: checkpoint_manager_base_t, solver_t
      class(checkpoint_manager_base_t), intent(inout) :: self
      class(solver_t), intent(inout) :: solver
      integer, intent(in), optional :: comm
    end subroutine handle_restart_interface

    subroutine handle_io_step_interface(self, solver, timestep, comm)
      import :: checkpoint_manager_base_t, solver_t
      class(checkpoint_manager_base_t), intent(inout) :: self
      class(solver_t), intent(in) :: solver
      integer, intent(in) :: timestep
      integer, intent(in), optional :: comm
    end subroutine handle_io_step_interface

    subroutine finalise_interface(self)
      import :: checkpoint_manager_base_t
      class(checkpoint_manager_base_t), intent(inout) :: self
    end subroutine finalise_interface
  end interface
end module m_checkpoint_manager_base

module m_checkpoint_manager_impl
!! Implementation of checkpoint manager when ADIOS2 is enabled
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Allreduce, &
                 MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MIN, MPI_Abort, &
                 MPI_INTEGER, MPI_SUM, MPI_IN_PLACE
  use m_common, only: dp, i8, DIR_C, VERT, get_argument
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_adios2_io, only: adios2_writer_t, adios2_reader_t, adios2_file_t, &
                         adios2_mode_write, adios2_mode_read
  use m_config, only: checkpoint_config_t
  use m_checkpoint_manager_base, only: checkpoint_manager_base_t
  
  implicit none
  
  private
  public :: checkpoint_manager_adios2_t
  
  type, extends(checkpoint_manager_base_t) :: checkpoint_manager_adios2_t
    type(adios2_writer_t) :: adios2_writer
    integer :: last_checkpoint_step = -1
    integer, dimension(3) :: output_stride = [2, 2, 2]  !! Spatial stride for snapshot output
    integer :: output_precision = dp                    !! Output precision for snapshot
  contains
    procedure :: init
    procedure :: handle_restart
    procedure :: handle_io_step
    procedure :: finalise
    procedure, private :: write_checkpoint
    procedure, private :: write_snapshot
    procedure, private :: restart_checkpoint
    procedure, private :: configure
    procedure, private :: write_fields
    procedure, private :: stride_data
  end type checkpoint_manager_adios2_t
  
contains

  subroutine init(self, comm)
    !! Initialise checkpoint manager with a reference to an ADIOS2 writer
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    integer, intent(in) :: comm
    
    call self%adios2_writer%init(comm, "checkpoint_writer")
    
    self%checkpoint_cfg = checkpoint_config_t()
    call self%checkpoint_cfg%read(nml_file=get_argument(1))
    
    call self%configure( &
      self%checkpoint_cfg%checkpoint_freq, &
      self%checkpoint_cfg%snapshot_freq, &
      self%checkpoint_cfg%keep_checkpoint, &
      self%checkpoint_cfg%checkpoint_prefix, &
      self%checkpoint_cfg%snapshot_prefix, &
      comm)
  end subroutine init

  subroutine handle_restart(self, solver, comm)
    !! Check if a restart is needed and handle it
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm

    character(len=256) :: restart_file
    integer :: restart_timestep
    real(dp) :: restart_time

    if (.not. self%checkpoint_cfg%restart_from_checkpoint) return
    self%is_restart = .true.

    restart_file = trim(self%checkpoint_cfg%restart_file)
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
  
  subroutine handle_io_step(self, solver, timestep, comm)
    !! Method to handle checkpoint and snapshot writing at a given timestep
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    integer :: comm_to_use

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm

    call self%write_checkpoint(solver, timestep, comm)
    call self%write_snapshot(solver, timestep, comm)
  end subroutine handle_io_step
  
  subroutine configure(self, checkpoint_freq, snapshot_freq, keep_checkpoint, &
                     checkpoint_prefix, snapshot_prefix, comm)
    !! Configure checkpoint and snapshot settings
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    integer, intent(in), optional :: checkpoint_freq, snapshot_freq
    logical, intent(in), optional :: keep_checkpoint
    character(len=*), intent(in), optional :: checkpoint_prefix, &
                                              snapshot_prefix
    integer, intent(in), optional :: comm
    
    integer :: comm_to_use, myrank, ierr
    
    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)
    
    if (present(checkpoint_freq)) &
        self%checkpoint_cfg%checkpoint_freq = checkpoint_freq
    if (present(snapshot_freq)) &
        self%checkpoint_cfg%snapshot_freq = snapshot_freq
    if (present(keep_checkpoint)) &
        self%checkpoint_cfg%keep_checkpoint = keep_checkpoint
    if (present(checkpoint_prefix)) &
        self%checkpoint_cfg%checkpoint_prefix = checkpoint_prefix
    if (present(snapshot_prefix)) &
        self%checkpoint_cfg%snapshot_prefix = snapshot_prefix
    
    if (myrank == 0) then
      if (self%checkpoint_cfg%checkpoint_freq > 0) then
        print *, 'Checkpoint frequency: ', self%checkpoint_cfg%checkpoint_freq
        print *, 'Keep all checkpoints: ', self%checkpoint_cfg%keep_checkpoint
        print *, 'Checkpoint prefix: ', &
                  trim(self%checkpoint_cfg%checkpoint_prefix)
      end if

      if (self%checkpoint_cfg%snapshot_freq > 0) then
        print *, 'Snapshot frequency: ', self%checkpoint_cfg%snapshot_freq
        print *, 'Snapshot prefix: ', trim(self%checkpoint_cfg%snapshot_prefix)
      end if
    end if
  end subroutine configure
  
  subroutine write_checkpoint(self, solver, timestep, comm)
    !! Write a checkpoint file for simulation restart
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=256) :: filename
    type(adios2_file_t) :: file
    integer :: ierr, myrank
    integer :: comm_to_use
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    real(dp) :: simulation_time

    if (self%checkpoint_cfg%checkpoint_freq <= 0) return
    if (mod(timestep, self%checkpoint_cfg%checkpoint_freq) /= 0) return

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)

    write(filename, '(A,A,I0.6,A)') &
      trim(self%checkpoint_cfg%checkpoint_prefix), '_', timestep, '.bp'
    if (myrank == 0) print *, 'Writing checkpoint: ', trim(filename)

    file = self%adios2_writer%open(filename, adios2_mode_write, comm_to_use)
    call self%adios2_writer%begin_step(file)

    simulation_time = timestep * solver%dt
    call self%adios2_writer%write_data("timestep", timestep, file)
    call self%adios2_writer%write_data("time", real(simulation_time, dp), file)
    call self%adios2_writer%write_data("dt", real(solver%dt, dp), file)

    call self%write_fields(field_names, solver, file, use_stride=.false.)
    call self%adios2_writer%close(file)

    ! Remove old checkpoint if configured to keep only the latest
    if (.not. self%checkpoint_cfg%keep_checkpoint &
      .and. self%last_checkpoint_step > 0) then
      write(filename, '(A,A,I0.6,A)') &
        trim(self%checkpoint_cfg%checkpoint_prefix), '_', &
        self%last_checkpoint_step, '.bp'
      if (myrank == 0) then
        call execute_command_line('rm -rf ' // trim(filename), exitstat=ierr)
        if (ierr /= 0) then
          print *, 'Warning: Failed to remove old checkpoint: ', trim(filename)
        end if
      end if
    end if

    self%last_checkpoint_step = timestep
  end subroutine write_checkpoint
  
  subroutine write_snapshot(self, solver, timestep, comm)
   !! Write a snapshot file for visualisation
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    integer :: myrank, ierr, i, j, k
    integer :: comm_to_use
    character(len=256) :: filename
    type(adios2_file_t) :: file
    integer :: nx, ny, nz
    integer :: dims(3), data_loc
    real(dp), dimension(3) :: min_coords, max_coords, local_min, local_max
    real(dp) :: dx, dy, dz
    real(dp), dimension(3) :: coords
    integer :: global_dims(3)
    integer :: global_nx, global_ny, global_nz

    if (self%checkpoint_cfg%snapshot_freq <= 0) return
    if (mod(timestep, self%checkpoint_cfg%snapshot_freq) /= 0) return

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)

    write(filename, '(A,A,I0.6,A)') &
      trim(self%checkpoint_cfg%snapshot_prefix), '_', timestep, '.bp'
    if (myrank == 0) print *, 'Writing snapshot: ', trim(filename)
    file = self%adios2_writer%open(filename, adios2_mode_write, comm_to_use)
    call self%adios2_writer%begin_step(file)

    data_loc = solver%u%data_loc
    dims = solver%mesh%get_dims(data_loc)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    global_dims = solver%mesh%get_global_dims(data_loc)
    global_nx = global_dims(1)
    global_ny = global_dims(2)
    global_nz = global_dims(3)

    local_min = [huge(0.0_dp), huge(0.0_dp), huge(0.0_dp)]
    local_max = [-huge(0.0_dp), -huge(0.0_dp), -huge(0.0_dp)]

    ! sample corners of the local domain to find bounds
    do k = 1, nz, max(1, nz-1)  ! just check first and last
      do j = 1, ny, max(1, ny-1)
        do i = 1, nx, max(1, nx-1)
          coords = solver%mesh%get_coordinates(i, j, k)

          local_min(1) = min(local_min(1), coords(1))
          local_min(2) = min(local_min(2), coords(2))
          local_min(3) = min(local_min(3), coords(3))

          local_max(1) = max(local_max(1), coords(1))
          local_max(2) = max(local_max(2), coords(2))
          local_max(3) = max(local_max(3), coords(3))
        end do
      end do
    end do

    call MPI_Allreduce(local_min, min_coords, 3, &
         MPI_DOUBLE_PRECISION, MPI_MIN, comm_to_use, ierr)
    call MPI_Allreduce(local_max, max_coords, 3, &
         MPI_DOUBLE_PRECISION, MPI_MAX, comm_to_use, ierr)

    ! calculate approximate spacing
    dx = (max_coords(1) - min_coords(1)) / max(1, global_nx - 1)
    dy = (max_coords(2) - min_coords(2)) / max(1, global_ny - 1)
    dz = (max_coords(3) - min_coords(3)) / max(1, global_nz - 1)

    ! paraview-specific mesh metadata
    if (myrank == 0) then
      call self%adios2_writer%write_attribute("mesh/time_varying", &
                                             "false", file)

      call self%adios2_writer%write_attribute("mesh/format", "vtk", file)
      call self%adios2_writer%write_attribute("mesh/type", "structured", file)
      call self%adios2_writer%write_attribute( &
        "mesh/origin/x", trim(real_to_str(min_coords(1))), file)
      call self%adios2_writer%write_attribute( &
        "mesh/origin/y", trim(real_to_str(min_coords(2))), file)
      call self%adios2_writer%write_attribute( &
        "mesh/origin/z", trim(real_to_str(min_coords(3))), file)
      call self%adios2_writer%write_attribute( &
        "mesh/spacing/x", trim(real_to_str(dx)), file)
      call self%adios2_writer%write_attribute( &
        "mesh/spacing/y", trim(real_to_str(dy)), file)
      call self%adios2_writer%write_attribute( &
        "mesh/spacing/z", trim(real_to_str(dz)), file)

      call self%adios2_writer%write_attribute("ParaView/MeshDimensions", &
                                              "u", file)

      call self%adios2_writer%write_attribute("mesh/dimensionality", "3", file)
      call self%adios2_writer%write_attribute("u/mesh", "mesh", file)
      call self%adios2_writer%write_attribute("v/mesh", "mesh", file)
      call self%adios2_writer%write_attribute("w/mesh", "mesh", file)
    end if

    call self%write_fields(field_names, solver, file, use_stride=.true.)
    call self%adios2_writer%close(file)
  end subroutine write_snapshot

  function stride_data(self, input_data, dims, stride) result(strided_data)
    !! Stride the input data based on the specified stride
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    real(dp), dimension(:,:,:), intent(in) :: input_data
    integer, dimension(3), intent(in) :: dims
    integer, dimension(3), intent(in) :: stride

    real(dp), dimension(:,:,:), allocatable :: strided_data
    integer :: i, j, k, is, js, ks
    integer :: strided_dims(3)

    strided_dims(1) = (dims(1) + stride(1) - 1) / stride(1)
    strided_dims(2) = (dims(2) + stride(2) - 1) / stride(2)
    strided_dims(3) = (dims(3) + stride(3) - 1) / stride(3)

    allocate(strided_data(strided_dims(1), strided_dims(2), &
                          strided_dims(3)))

    is = 1
    do i = 1, dims(1), stride(1)
      js = 1
      do j = 1, dims(2), stride(2)
        ks = 1
        do k = 1, dims(3), stride(3)
          strided_data(is, js, ks) = input_data(i, j, k)
          ks = ks + 1
        end do
        js = js + 1
      end do
      is = is + 1
    end do

  end function stride_data
  
  subroutine restart_checkpoint( &
    self, solver, filename, timestep, restart_time, comm &
    )
    !! Restart simulation state from checkpoint file
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    character(len=*), intent(in) :: filename
    integer, intent(out) :: timestep
    real(dp), intent(out) :: restart_time
    integer, intent(in) :: comm

    type(adios2_reader_t) :: reader
    type(adios2_file_t) :: file
    real(dp), allocatable :: field_data(:,:,:)
    integer :: i, ierr
    integer :: dims(3)
    class(field_t), pointer :: host_field
    integer(i8), dimension(3) :: start_dims, count_dims
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    logical :: file_exists

    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      if (solver%mesh%par%is_root()) then
        print *, 'ERROR: Checkpoint file not found: ', trim(filename)
      end if
      call MPI_Abort(comm, 1, ierr)
      return
    end if

    call reader%init(comm, "checkpoint_reader")

    dims = solver%mesh%get_dims(VERT)
    start_dims = int(solver%mesh%par%n_offset, i8)
    count_dims = int(dims, i8)

    file = reader%open(filename, adios2_mode_read, comm)
    call reader%begin_step(file)
    call reader%read_data("timestep", timestep, file)
    call reader%read_data("time", restart_time, file)

    do i = 1, size(field_names)
      if (allocated(field_data)) deallocate(field_data)

      call reader%read_data(field_names(i), &
           field_data, file, start_dims, count_dims)

      if (.not. allocated(field_data)) then
        if (solver%mesh%par%is_root()) then
          print *, 'ERROR: Failed to read field: ', trim(field_names(i))
        end if
        call MPI_Abort(comm, 1, ierr)
        return
      end if

      select case (trim(field_names(i)))
        case ("u")
          host_field => solver%host_allocator%get_block(DIR_C)
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = field_data
          call solver%backend%set_field_data(solver%u, host_field%data)
          call solver%host_allocator%release_block(host_field)
        case ("v")
          host_field => solver%host_allocator%get_block(DIR_C)
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = field_data
          call solver%backend%set_field_data(solver%v, host_field%data)
          call solver%host_allocator%release_block(host_field)
        case ("w")
          host_field => solver%host_allocator%get_block(DIR_C)
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = field_data
          call solver%backend%set_field_data(solver%w, host_field%data)
          call solver%host_allocator%release_block(host_field)
        case default
          call self%adios2_writer%handle_error( & 
            1, "Invalid field name"//trim(field_names(i)))
      end select
    end do

    call solver%u%set_data_loc(VERT)
    call solver%v%set_data_loc(VERT)
    call solver%w%set_data_loc(VERT)

    call reader%close(file)
    call reader%finalise()
    if (allocated(field_data)) deallocate(field_data)
  end subroutine restart_checkpoint
  
  subroutine write_fields(self, field_names, solver, file, use_stride)
    !! Write field data, optionally with striding
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    character(len=*), dimension(:), intent(in) :: field_names
    class(solver_t), intent(in) :: solver
    type(adios2_file_t), intent(inout) :: file
    logical, intent(in), optional :: use_stride

    class(field_t), pointer :: host_field, field_ptr
    integer :: i_field, dims(3), data_loc
    integer(i8), dimension(3)  :: shape_dims, start_dims, count_dims
    integer(i8), dimension(3) :: output_shape, output_start, output_count
    real(dp), dimension(:,:,:), allocatable :: data_to_write
    logical :: apply_stride

    apply_stride = .false.
    if (present(use_stride)) apply_stride = use_stride


    do i_field = 1, size(field_names)
      select case (trim(field_names(i_field)))
        case ("u")
          field_ptr => solver%u
          data_loc = solver%u%data_loc
        case ("v")
          field_ptr => solver%v
          data_loc = solver%v%data_loc
        case ("w")
          field_ptr => solver%w
          data_loc = solver%w%data_loc
        case default
          call self%adios2_writer%handle_error( & 
            1, "Invalid field name"//trim(field_names(i_field)))
          cycle
      end select

      dims = solver%mesh%get_dims(data_loc)
      shape_dims = int(solver%mesh%get_global_dims(data_loc), i8)
      start_dims = int(solver%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)

      if (apply_stride) then
        output_shape = (shape_dims + self%output_stride - 1) / &
                        self%output_stride
        output_start = start_dims / self%output_stride
        output_count = (count_dims + self%output_stride - 1) / &
                       self%output_stride
      else
        output_shape = shape_dims
        output_start = start_dims
        output_count = count_dims
      end if

      host_field => solver%host_allocator%get_block(DIR_C, data_loc)

      if (allocated(data_to_write)) deallocate(data_to_write)

      call solver%backend%get_field_data(host_field%data, field_ptr)

      if (apply_stride) then
        data_to_write = self%stride_data(&
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
          dims, self%output_stride &
          )
      else
        allocate(data_to_write(dims(1), dims(2), dims(3)))
        data_to_write = host_field%data(1:dims(1), 1:dims(2), 1:dims(3))
      end if

      select case (trim(field_names(i_field)))
        case ("u")
          call self%adios2_writer%write_data( &
            "u", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
            file, shape_dims, start_dims, count_dims &
            )
        case ("v")
          call self%adios2_writer%write_data( &
            "v", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
            file, shape_dims, start_dims, count_dims &
            )
        case ("w")
          call self%adios2_writer%write_data( &
            "w", host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
            file, shape_dims, start_dims, count_dims &
            )
      end select

      if (allocated(data_to_write)) deallocate(data_to_write)
    end do
    call solver%host_allocator%release_block(host_field)
  end subroutine write_fields
  
  subroutine finalise(self)
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    
    call self%adios2_writer%finalise()
  end subroutine finalise

  function real_to_str(val) result(str)
    real(dp), intent(in) :: val
    character(len=30) :: str
    
    write(str, '(ES15.7E3)') val
    str = adjustl(str)
  end function real_to_str
end module m_checkpoint_manager_impl

module m_checkpoint_manager
  !! Public facade and factory function for checkpoint manager
  use m_checkpoint_manager_impl, only: checkpoint_manager_adios2_t
  use m_solver, only: solver_t
  
  implicit none
  
  private
  public :: checkpoint_manager_t, create_checkpoint_manager
  
  type, public :: checkpoint_manager_t
    type(checkpoint_manager_adios2_t) :: impl
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
    
    is_restart = self%impl%is_restart
  end function cm_is_restart

end module m_checkpoint_manager
