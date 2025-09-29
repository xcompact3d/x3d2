module m_checkpoint_manager_base
!! Base module defining abstract interface for checkpoint functionality
  use m_solver, only: solver_t
  use m_config, only: checkpoint_config_t

  implicit none

  private
  public :: checkpoint_manager_base_t

  type, abstract :: checkpoint_manager_base_t
    type(checkpoint_config_t) :: checkpoint_cfg
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
!! Implementation of checkpoint manager using generic I/O abstraction
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Abort
  use m_common, only: dp, i8, DIR_C, VERT, get_argument
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: allocate_io_writer, io_session_t
  use m_io_base, only: io_writer_t, io_file_t, io_mode_write
  use m_config, only: checkpoint_config_t
  use m_checkpoint_manager_base, only: checkpoint_manager_base_t
  use m_io_field_utils, only: field_buffer_map_t, field_ptr_t, stride_data, stride_data_to_buffer, &
                              get_output_dimensions, generate_coordinates, &
                              setup_field_arrays, cleanup_field_arrays

  implicit none

  private
  public :: checkpoint_manager_impl_t

  type, extends(checkpoint_manager_base_t) :: checkpoint_manager_impl_t
    class(io_writer_t), pointer :: checkpoint_writer => null()               !! Writer for checkpoints
    class(io_writer_t), pointer :: snapshot_writer => null()                 !! Writer for snapshots
    integer :: last_checkpoint_step = -1
    integer, dimension(3) :: output_stride = [1, 1, 1]              !! Stride factors for snapshots (default: full resolution)
    integer, dimension(3) :: full_resolution = [1, 1, 1]           !! Full resolution factors for checkpoints (no downsampling)
    type(field_buffer_map_t), allocatable :: field_buffers(:)       !! Dynamic field buffer mapping for race-free async I/O
    real(dp), dimension(:, :, :), allocatable :: coords_x, coords_y, coords_z
    integer(i8), dimension(3) :: last_shape_dims = 0
    integer, dimension(3) :: last_stride_factors = 0
    integer(i8), dimension(3) :: last_output_shape = 0
    character(len=4096) :: vtk_xml = ""                             !! VTK XML string for ParaView compatibility
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
    procedure, private :: cleanup_output_buffers
    procedure, private :: generate_vtk_xml
  end type checkpoint_manager_impl_t

contains

  subroutine init(self, comm)
    !! Initialise checkpoint manager with a reference to a generic I/O writer
    class(checkpoint_manager_impl_t), intent(inout) :: self
    integer, intent(in) :: comm

    call allocate_io_writer(self%checkpoint_writer)
    call self%checkpoint_writer%init(comm, "checkpoint_writer")
    call allocate_io_writer(self%snapshot_writer)
    call self%snapshot_writer%init(comm, "snapshot_writer")

    self%checkpoint_cfg = checkpoint_config_t()
    call self%checkpoint_cfg%read(nml_file=get_argument(1))

    call self%configure( &
      self%checkpoint_cfg%checkpoint_freq, &
      self%checkpoint_cfg%snapshot_freq, &
      self%checkpoint_cfg%keep_checkpoint, &
      self%checkpoint_cfg%checkpoint_prefix, &
      self%checkpoint_cfg%snapshot_prefix, &
      self%checkpoint_cfg%output_stride, &
      comm &
      )
  end subroutine init

  subroutine handle_restart(self, solver, comm)
    !! Check if a restart is needed and handle it
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    integer, intent(in), optional :: comm

    character(len=256) :: restart_file
    integer :: restart_timestep
    real(dp) :: restart_time

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
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    integer :: comm_to_use

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm

    call self%write_checkpoint(solver, timestep, comm_to_use)
    call self%write_snapshot(solver, timestep, comm_to_use)
  end subroutine handle_io_step

  subroutine configure( &
    self, checkpoint_freq, snapshot_freq, keep_checkpoint, &
    checkpoint_prefix, snapshot_prefix, output_stride, comm &
    )
    !! Configure checkpoint and snapshot settings
    class(checkpoint_manager_impl_t), intent(inout) :: self
    integer, intent(in), optional :: checkpoint_freq, snapshot_freq
    logical, intent(in), optional :: keep_checkpoint
    character(len=*), intent(in), optional :: checkpoint_prefix, &
                                              snapshot_prefix
    integer, dimension(3), intent(in), optional :: output_stride
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
    if (present(output_stride)) &
      self%output_stride = self%checkpoint_cfg%output_stride

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
        print *, 'Output stride: ', self%output_stride
      end if
    end if
  end subroutine configure

  subroutine write_checkpoint(self, solver, timestep, comm)
    !! Write a checkpoint file for simulation restart
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=256) :: filename, temp_filename, old_filename
    integer :: ierr, myrank
    integer :: comm_to_use, i
    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    real(dp) :: simulation_time
    logical :: file_exists
    type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
    integer, parameter :: num_fields = size(field_names)
    integer :: data_loc
    type(io_session_t) :: io_session

    if (self%checkpoint_cfg%checkpoint_freq <= 0) return
    if (mod(timestep, self%checkpoint_cfg%checkpoint_freq) /= 0) return

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%checkpoint_cfg%checkpoint_prefix), '_', timestep, '.bp'
    write (temp_filename, '(A,A)') &
      trim(self%checkpoint_cfg%checkpoint_prefix), '_temp.bp'
    if (myrank == 0) print *, 'Writing checkpoint: ', trim(filename)

    call io_session%open(temp_filename, comm_to_use, io_mode_write)

    simulation_time = timestep*solver%dt
    data_loc = solver%u%data_loc
    call io_session%write_data("timestep", timestep)
    call io_session%write_data("time", real(simulation_time, dp))
    call io_session%write_data("dt", real(solver%dt, dp))
    call io_session%write_data("data_loc", data_loc)

    call setup_field_arrays(solver, field_names, field_ptrs, host_fields)

    call self%write_fields( &
      field_names, host_fields, &
      solver, io_session, data_loc, use_stride=.false. &
      )

    call io_session%close()

    call cleanup_field_arrays(solver, field_ptrs, host_fields)

    if (myrank == 0) then
      ! Move temporary file to final checkpoint filename
      call execute_command_line('mv '//trim(temp_filename)//' '// &
                                trim(filename))

      inquire (file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
        print *, 'ERROR: Checkpoint file not created: ', trim(filename)
      end if

      ! Remove old checkpoint if configured to keep only the latest
      if (.not. self%checkpoint_cfg%keep_checkpoint &
          .and. self%last_checkpoint_step > 0) then
        write (old_filename, '(A,A,I0.6,A)') &
          trim(self%checkpoint_cfg%checkpoint_prefix), '_', &
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

  subroutine write_snapshot(self, solver, timestep, comm)
   !! Write a snapshot file for visualisation
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    integer :: myrank, ierr, comm_to_use, i
    character(len=256) :: filename
    integer(i8), dimension(3) :: output_shape_dims
    integer, dimension(3) :: global_dims, output_dims
    type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
    integer, parameter :: num_fields = size(field_names)
    real(dp), dimension(3) :: origin, original_spacing, output_spacing
    real(dp) :: simulation_time
    logical :: snapshot_uses_stride = .true. ! snapshots always use striding
    type(io_session_t) :: io_session

    if (self%checkpoint_cfg%snapshot_freq <= 0) return
    if (mod(timestep, self%checkpoint_cfg%snapshot_freq) /= 0) return

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%checkpoint_cfg%snapshot_prefix), '_', timestep, '.bp'

    global_dims = solver%mesh%get_global_dims(VERT)
    origin = solver%mesh%get_coordinates(1, 1, 1)
    original_spacing = solver%mesh%geo%d

    if (snapshot_uses_stride) then
      output_spacing = original_spacing*real(self%output_stride, dp)
      do i = 1, size(global_dims)
        output_dims(i) = (global_dims(i) + self%output_stride(i) - 1)/ &
                         self%output_stride(i)
      end do
    else
      output_spacing = original_spacing
      output_dims = global_dims
    end if
    output_shape_dims = int(output_dims, i8)

    call self%generate_vtk_xml( &
      output_shape_dims, field_names, origin, output_spacing &
      )

    call io_session%open(filename, comm_to_use, io_mode_write)
    if (myrank == 0) print *, 'Creating snapshot file: ', trim(filename)

    ! Write VTK XML attributes for ParaView compatibility
    if (myrank == 0) then
      call io_session%write_attribute("vtk.xml", self%vtk_xml)
    end if

    simulation_time = timestep*solver%dt
    if (myrank == 0) then
      print *, 'Writing snapshot for time =', simulation_time, &
        ' iteration =', timestep
    end if

    call io_session%write_data("time", real(simulation_time, dp))

    call setup_field_arrays(solver, field_names, field_ptrs, host_fields)

    call self%write_fields( &
      field_names, host_fields, &
      solver, io_session, solver%u%data_loc, &
      use_stride=snapshot_uses_stride &
      )

    call io_session%close()

    call cleanup_field_arrays(solver, field_ptrs, host_fields)
  end subroutine write_snapshot

  subroutine generate_vtk_xml(self, dims, fields, origin, spacing)
    !! Generate VTK XML string for ImageData format for ParaView's ADIOS2VTXReader
    class(checkpoint_manager_impl_t), intent(inout) :: self
    integer(i8), dimension(3), intent(in) :: dims
    character(len=*), dimension(:), intent(in) :: fields
    real(dp), dimension(3), intent(in) :: origin, spacing

    character(len=4096) :: xml
    character(len=96) :: extent_str, origin_str, spacing_str
    integer :: i

    ! VTK uses (x,y,z) order, extent defines grid size from 0 to N-1
    write (extent_str, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      '0 ', dims(3) - 1, ' 0 ', dims(2) - 1, ' 0 ', dims(1) - 1

    write (origin_str, '(G0, 1X, G0, 1X, G0)') origin
    write (spacing_str, '(G0, 1X, G0, 1X, G0)') spacing

    xml = '<?xml version="1.0"?>'//new_line('a')// &
          '<VTKFile type="ImageData" version="0.1">'//new_line('a')// &
          '  <ImageData WholeExtent=" '//trim(adjustl(extent_str))//'" ' &
          //'Origin="'//trim(adjustl(origin_str))//'" '// &
          'Spacing="'//trim(adjustl(spacing_str))//'">'//new_line('a')// &
          '    <Piece Extent="'//trim(adjustl(extent_str))//'">' &
          //new_line('a')//'      <PointData>'//new_line('a')

    do i = 1, size(fields)
      xml = trim(xml)//'      <DataArray Name="'//trim(fields(i))//'">'// &
            trim(fields(i))//'</DataArray>'//new_line('a')
    end do

    xml = trim(xml)//'        <DataArray Name="TIME">time</DataArray>' &
          //new_line('a')

    xml = trim(xml)//'      </PointData>'//new_line('a')// &
          '    </Piece>'//new_line('a')// &
          '  </ImageData>'//new_line('a')// &
          '</VTKFile>'

    self%vtk_xml = xml
  end subroutine generate_vtk_xml

  subroutine restart_checkpoint( &
    self, solver, filename, timestep, restart_time, comm &
    )
    !! Restart simulation state from checkpoint file
    class(checkpoint_manager_impl_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    character(len=*), intent(in) :: filename
    integer, intent(out) :: timestep
    real(dp), intent(out) :: restart_time
    integer, intent(in) :: comm

    type(io_session_t) :: io_session
    integer :: i, ierr
    integer :: dims(3)
    integer(i8), dimension(3) :: start_dims, count_dims
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    logical :: file_exists
    integer :: data_loc

    inquire (file=filename, exist=file_exists)
    if (.not. file_exists) then
      if (solver%mesh%par%is_root()) then
        print *, 'ERROR: Checkpoint file not found: ', trim(filename)
      end if
      call MPI_Abort(comm, 1, ierr)
      return
    end if

    call io_session%open(filename, comm)
    call io_session%read_data("timestep", timestep)
    call io_session%read_data("time", restart_time)
    call io_session%read_data("data_loc", data_loc)

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

      allocate(field_data_u(count_dims(1), count_dims(2), count_dims(3)))
      allocate(field_data_v(count_dims(1), count_dims(2), count_dims(3)))
      allocate(field_data_w(count_dims(1), count_dims(2), count_dims(3)))
      
      call io_session%read_data("u", field_data_u)
      call io_session%read_data("v", field_data_v)
      call io_session%read_data("w", field_data_w)
      
      call solver%backend%set_field_data(solver%u, field_data_u)
      call solver%backend%set_field_data(solver%v, field_data_v)
      call solver%backend%set_field_data(solver%w, field_data_w)
    end block

    call io_session%close()
  end subroutine restart_checkpoint

  subroutine write_fields( &
    self, field_names, host_fields, solver, io_session, data_loc, &
    use_stride &
    )
    !! Write field data, optionally with striding, using separate buffers for true async I/O
    class(checkpoint_manager_impl_t), intent(inout) :: self
    character(len=*), dimension(:), intent(in) :: field_names
    class(field_ptr_t), dimension(:), target, intent(in) :: host_fields
    class(solver_t), intent(in) :: solver
    type(io_session_t), intent(inout) :: io_session
    integer, intent(in) :: data_loc
    logical, intent(in), optional :: use_stride

    integer :: i_field
    logical :: apply_stride

    apply_stride = .false.
    if (present(use_stride)) apply_stride = use_stride

    ! pre-allocate separate buffers for each field to enable race-free async I/O
    if (apply_stride) then
      call prepare_field_buffers(solver, self%output_stride, field_names, &
                                 data_loc)
    else
      call prepare_field_buffers(solver, self%full_resolution, field_names, &
                                 data_loc)
    end if

    do i_field = 1, size(field_names)
      call write_single_field( &
        trim(field_names(i_field)), host_fields(i_field)%ptr, &
        data_loc)
    end do

  contains

    subroutine prepare_field_buffers( &
      solver, stride_factors, field_names, data_loc &
      )
      class(solver_t), intent(in) :: solver
      integer, dimension(3), intent(in) :: stride_factors
      character(len=*), dimension(:), intent(in) :: field_names
      integer, intent(in) :: data_loc
      integer :: dims(3), output_dims_local(3), i
      integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
      integer(i8), dimension(3) :: output_shape, output_start, output_count

      dims = solver%mesh%get_dims(data_loc)
      shape_dims = int(solver%mesh%get_global_dims(data_loc), i8)
      start_dims = int(solver%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)

      call get_output_dimensions( &
        shape_dims, start_dims, count_dims, stride_factors, &
        output_shape, output_start, output_count, &
        output_dims_local, &
        self%last_shape_dims, self%last_stride_factors, self%last_output_shape &
        )

      if (allocated(self%field_buffers)) deallocate (self%field_buffers)
      allocate (self%field_buffers(size(field_names)))

      do i = 1, size(field_names)
        self%field_buffers(i)%field_name = trim(field_names(i))
        allocate ( &
          self%field_buffers(i)%buffer( &
          output_dims_local(1), &
          output_dims_local(2), &
          output_dims_local(3)))
      end do
    end subroutine prepare_field_buffers

    subroutine write_single_field(field_name, host_field, data_loc)
      character(len=*), intent(in) :: field_name
      class(field_t), pointer :: host_field
      integer, intent(in) :: data_loc

      integer, dimension(3) :: output_dims_local, stride_factors
      integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
      integer(i8), dimension(3) :: output_shape, output_start, output_count
      integer :: dims(3)
      integer :: buffer_idx
      logical :: buffer_found

      dims = solver%mesh%get_dims(data_loc)
      shape_dims = int(solver%mesh%get_global_dims(data_loc), i8)
      start_dims = int(solver%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)

      if (apply_stride) then
        stride_factors = self%output_stride
      else
        !no striding for checkpointing
        stride_factors = self%full_resolution
      end if

      call get_output_dimensions( &
        shape_dims, start_dims, count_dims, stride_factors, &
        output_shape, output_start, output_count, &
        output_dims_local, &
        self%last_shape_dims, self%last_stride_factors, self%last_output_shape &
        )

      ! find the matching buffer for this field
      buffer_found = .false.
      do buffer_idx = 1, size(self%field_buffers)
        if (trim(self%field_buffers(buffer_idx)%field_name) &
            == trim(field_name)) then
          buffer_found = .true.
          exit
        end if
      end do

      if (buffer_found) then
        ! use the dedicated buffer for this field for true async I/O
        call stride_data_to_buffer( &
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), dims, &
          stride_factors, self%field_buffers(buffer_idx)%buffer, &
          output_dims_local &
          )

        call io_session%write_data( &
          field_name, self%field_buffers(buffer_idx)%buffer, &
          start_dims=output_start, count_dims=output_count &
          )
      else
        ! ERROR: All fields must have pre-allocated buffers for race-free async I/O
        print *, 'INTERNAL ERROR: No dedicated buffer found for field: ', trim(field_name)
        print *, 'Available buffers:'
        do buffer_idx = 1, size(self%field_buffers)
          print *, '  - ', trim(self%field_buffers(buffer_idx)%field_name)
        end do
        print *, 'This violates the race-free async I/O design.'
        print *, 'All fields must be pre-allocated in prepare_field_buffers().'
        error stop 'Missing field buffer - race condition risk'
      end if

    end subroutine write_single_field

  end subroutine write_fields

  subroutine cleanup_output_buffers(self)
    !! Clean up dynamic field buffers
    class(checkpoint_manager_impl_t), intent(inout) :: self
    integer :: i

    if (allocated(self%field_buffers)) then
      do i = 1, size(self%field_buffers)
        if (allocated(self%field_buffers(i)%buffer)) then
          deallocate (self%field_buffers(i)%buffer)
        end if
      end do
      deallocate (self%field_buffers)
    end if
  end subroutine cleanup_output_buffers

  subroutine finalise(self)
    class(checkpoint_manager_impl_t), intent(inout) :: self

    call self%cleanup_output_buffers()
    call self%checkpoint_writer%finalise()
    call self%snapshot_writer%finalise()

    if (associated(self%checkpoint_writer)) deallocate(self%checkpoint_writer)
    if (associated(self%snapshot_writer)) deallocate(self%snapshot_writer)
  end subroutine finalise

end module m_checkpoint_manager_impl

module m_checkpoint_manager
  !! Public facade and factory function for checkpoint manager
  use m_checkpoint_manager_impl, only: checkpoint_manager_impl_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: checkpoint_manager_t

  type :: checkpoint_manager_t
    type(checkpoint_manager_impl_t) :: impl
  contains
    procedure :: init => cm_init
    procedure :: handle_restart => cm_handle_restart
    procedure :: handle_io_step => cm_handle_io_step
    procedure :: finalise => cm_finalise
    procedure :: is_restart => cm_is_restart
  end type checkpoint_manager_t

contains

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
