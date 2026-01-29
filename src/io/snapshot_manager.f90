module m_snapshot_manager
  !! Manages creation of simulation snapshots for post-processing and visualisation.
  !!
  !! This module periodically writes simulation data to files intended for
  !! analysis and visualisation. Unlike checkpoints (full-resolution for exact
  !! restarts), snapshots can be strided to reduce file size while retaining
  !! sufficient resolution for visualisation.
  !!
  !! **Key Differences from Checkpoints:**
  !! - **Purpose**: Visualisation/analysis vs exact restart
  !! - **Resolution**: Can be strided (e.g., every 2nd point) vs full resolution
  !! - **Frequency**: Typically more frequent than checkpoints
  !! - **File Management**: Single persistent file with multiple timesteps vs
  !!   separate files per checkpoint
  !!
  !! **Features:**
  !! - Configurable spatial striding to reduce output size
  !! - Persistent file handle (stays open across multiple writes)
  !! - Generates VTK-compatible XML for ParaView visualisation
  !! - Writes velocity fields at each snapshot interval
  !!
  !! **Configuration:**
  !! Controlled via `checkpoint_config_t` read from input namelist:
  !! - snapshot_freq: write interval (iterations)
  !! - snapshot_prefix: filename prefix
  !! - output_stride: spatial stride factors [sx, sy, sz]
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank
  use m_common, only: dp, i8, DIR_C, VERT, get_argument
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: writer_session_t
  use m_config, only: checkpoint_config_t
  use m_io_field_utils, only: field_buffer_map_t, field_ptr_t, &
                              setup_field_arrays, cleanup_field_arrays, &
                              stride_data_to_buffer, get_output_dimensions, &
                              prepare_field_buffers, cleanup_field_buffers, &
                              write_single_field_to_buffer

  implicit none

  private
  public :: snapshot_manager_t

  type :: snapshot_manager_t
    !! Manager for snapshot file operations (periodic visualisation output).
    !!
    !! Handles periodic writing of visualisation data with optional striding.
    !! Maintains a persistent file handle that stays open across multiple
    !! snapshot writes for efficient I/O.
    type(checkpoint_config_t) :: config                  !! Configuration settings
    integer, dimension(3) :: output_stride = [1, 1, 1]   !! Spatial stride factors [sx, sy, sz]
    type(field_buffer_map_t), allocatable :: field_buffers(:) !! Buffers for field data I/O
    integer(i8), dimension(3) :: last_shape_dims = 0     !! Shape dimensions from last write
    integer, dimension(3) :: last_stride_factors = 0     !! Stride factors from last write
    integer(i8), dimension(3) :: last_output_shape = 0   !! Output shape from last write
    character(len=4096) :: vtk_xml = ""                  !! VTK XML metadata for ParaView
    logical :: is_snapshot_file_open = .false.           !! File handle state flag
    type(writer_session_t) :: snapshot_writer            !! I/O session writer
    logical :: convert_to_sp = .false.                   !! Flag for single precision snapshots
  contains
    procedure :: init                          !! Initialise snapshot manager
    procedure :: handle_snapshot_step          !! Write snapshot if needed at timestep
    procedure :: finalise                      !! Clean up and finalise
    procedure, private :: write_snapshot       !! Write snapshot file (internal)
    procedure, private :: write_fields         !! Write field data to file (internal)
    procedure, private :: cleanup_output_buffers !! Free output buffers (internal)
    procedure, private :: generate_vtk_xml     !! Generate VTK XML metadata (internal)
    procedure, private :: open_snapshot_file   !! Open snapshot file (internal)
    procedure, private :: close_snapshot_file  !! Close snapshot file (internal)
  end type snapshot_manager_t

contains

  subroutine init(self, comm)
    !! Initialise snapshot manager from configuration.
    !!
    !! Reads snapshot settings from input namelist and configures
    !! output if snapshot frequency is positive. Prints snapshot
    !! settings including stride factors on root process.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance
    integer, intent(in) :: comm                       !! MPI communicator

    self%config = checkpoint_config_t()
    call self%config%read(nml_file=get_argument(1))

    if (self%config%snapshot_freq > 0) then
      call configure_output(self, comm)
    end if
  end subroutine init

  subroutine configure_output(self, comm)
    !! Configure and print snapshot output settings.
    !!
    !! Displays snapshot configuration on root process including
    !! frequency, file prefix, and output stride factors.
    use m_io_backend, only: get_default_backend, IO_BACKEND_DUMMY
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance
    integer, intent(in) :: comm                       !! MPI communicator

    integer :: myrank, ierr

    call MPI_Comm_rank(comm, myrank, ierr)

    self%output_stride = self%config%output_stride
    self%convert_to_sp = self%config%snapshot_sp

    if (myrank == 0 .and. get_default_backend() /= IO_BACKEND_DUMMY) then
      print *, 'Snapshot frequency: ', self%config%snapshot_freq
      print *, 'Snapshot prefix: ', trim(self%config%snapshot_prefix)
      print *, 'Output stride: ', self%output_stride
      print *, 'Snapshot precision: ', merge('Single', 'Double', &
                                             self%config%snapshot_sp)
    end if
  end subroutine configure_output

  subroutine handle_snapshot_step(self, solver, timestep, comm)
    !! Write snapshot if frequency condition is met.
    !!
    !! Checks if current timestep is a snapshot interval (divisible by
    !! snapshot_freq) and writes snapshot if so. Called each timestep
    !! from main simulation loop.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance
    class(solver_t), intent(in) :: solver            !! Solver containing current state
    integer, intent(in) :: timestep                   !! Current timestep number
    integer, intent(in), optional :: comm             !! MPI communicator (optional)

    integer :: comm_to_use

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm

    call self%write_snapshot(solver, timestep, comm_to_use)
  end subroutine handle_snapshot_step

  subroutine write_snapshot(self, solver, timestep, comm)
    !! Write a snapshot file for visualisation.
    !!
    !! Uses a persistent file that stays open across multiple snapshots.
    !! Each snapshot is written as a separate timestep within the file.
    !! Data can be strided according to output_stride configuration.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance
    class(solver_t), intent(in) :: solver            !! Solver containing field data
    integer, intent(in) :: timestep                   !! Current timestep number
    integer, intent(in) :: comm                       !! MPI communicator

    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    integer :: myrank, ierr
    character(len=256) :: filename
    integer(i8), dimension(3) :: output_shape_dims
    integer, dimension(3) :: global_dims, output_dims
    type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
    real(dp), dimension(3) :: origin, original_spacing, output_spacing
    real(dp) :: simulation_time
    logical :: snapshot_uses_stride = .true.
    integer :: i

    if (self%config%snapshot_freq <= 0) return
    if (mod(timestep, self%config%snapshot_freq) /= 0) return

    call MPI_Comm_rank(comm, myrank, ierr)

    write (filename, '(A,A)') trim(self%config%snapshot_prefix), '.bp'

    ! Open snapshot file on first call (check for existence)
    if (.not. self%is_snapshot_file_open) then
      call self%open_snapshot_file(filename, comm)
    else
      ! For subsequent snapshots, begin a new step
      call self%snapshot_writer%begin_step()
    end if

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

    simulation_time = timestep*solver%dt
    if (self%snapshot_writer%is_session_functional() .and. myrank == 0) then
      print *, 'Writing snapshot for time =', simulation_time, &
        ' iteration =', timestep
    end if

    ! Write VTK XML attributes for ParaView compatibility (only on first step)
    if (timestep == self%config%snapshot_freq .and. myrank == 0) then
      call self%snapshot_writer%write_attribute("vtk.xml", self%vtk_xml)
    end if

    call self%snapshot_writer%write_data("time", real(simulation_time, dp))

    call setup_field_arrays(solver, field_names, field_ptrs, host_fields)

    call self%write_fields( &
      field_names, host_fields, &
      solver, self%snapshot_writer, solver%u%data_loc &
      )

    call self%snapshot_writer%end_step()

    call cleanup_field_arrays(solver, field_ptrs, host_fields)
  end subroutine write_snapshot

  subroutine generate_vtk_xml(self, dims, fields, origin, spacing)
    !! Generate VTK XML metadata for ParaView visualization (internal).
    !!
    !! Creates VTK ImageData XML string that describes the structured grid
    !! for ParaView's ADIOS2VTXReader. This enables direct visualization of
    !! ADIOS2 files in ParaView without conversion.
    !!
    !! **VTK ImageData Format:**
    !! - Defines structured rectilinear grid with uniform spacing
    !! - Extent: grid dimensions from 0 to N-1 in (z,y,x) order
    !! - Origin: physical coordinates of first grid point
    !! - Spacing: grid resolution (dx, dy, dz)
    !! - Point data: velocity fields (u, v, w) stored at grid points
    !!
    !! **Note:** VTK uses (x,y,z) order while X3D2 uses (z,y,x) internally,
    !! requiring dimension reordering in the extent string.
    class(snapshot_manager_t), intent(inout) :: self        !! Snapshot manager instance
    integer(i8), dimension(3), intent(in) :: dims           !! Grid dimensions [nx, ny, nz]
    character(len=*), dimension(:), intent(in) :: fields    !! Field names ["u", "v", "w"]
    real(dp), dimension(3), intent(in) :: origin            !! Grid origin [x0, y0, z0]
    real(dp), dimension(3), intent(in) :: spacing           !! Grid spacing [dx, dy, dz]

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

  subroutine write_fields( &
    self, field_names, host_fields, solver, writer_session, data_loc &
    )
    !! Write field data with optional striding for snapshots (internal).
    !!
    !! Writes field data with spatial striding to reduce file size while
    !! maintaining sufficient resolution for visualization. The procedure:
    !! 1. Prepare field buffers with configured stride factors
    !! 2. Calculate strided output dimensions and hyperslab selection
    !! 3. For each field (u, v, w):
    !!    - Copy strided field data to output buffer
    !!    - Write buffer to file with proper hyperslab parameters
    !!
    !! **Striding:** Unlike checkpoints (full resolution), snapshots can
    !! subsample data. For example, stride [2,2,2] writes every 2nd point
    !! in each direction, reducing file size by factor of 8.
    !!
    !! **Parallel I/O:** Each MPI rank writes its strided local subdomain
    !! using hyperslab selection to assemble the strided global field.
    class(snapshot_manager_t), intent(inout) :: self       !! Snapshot manager instance
    character(len=*), dimension(:), intent(in) :: field_names !! Field names ["u", "v", "w"]
    class(field_ptr_t), dimension(:), target, intent(in) :: host_fields !! Field pointers
    class(solver_t), intent(in) :: solver                  !! Solver containing mesh info
    type(writer_session_t), intent(inout) :: writer_session !! I/O writer session
    integer, intent(in) :: data_loc                        !! Data location (VERT or CELL)

    integer :: i_field
    integer(i8), dimension(3) :: output_start, output_count
    integer, dimension(3) :: output_dims_local

    ! Prepare buffers with striding for snapshots
    call prepare_field_buffers( &
      solver, self%output_stride, field_names, data_loc, &
      self%field_buffers, self%last_shape_dims, self%last_stride_factors, &
      self%last_output_shape &
      )

    ! Calculate output dimensions for writing
    call get_output_dimensions( &
      int(solver%mesh%get_global_dims(data_loc), i8), &
      int(solver%mesh%par%n_offset, i8), &
      int(solver%mesh%get_dims(data_loc), i8), &
      self%output_stride, &
      self%last_output_shape, output_start, output_count, &
      output_dims_local, &
      self%last_shape_dims, self%last_stride_factors, &
      self%last_output_shape &
      )

    do i_field = 1, size(field_names)
      call write_single_field_to_buffer( &
        trim(field_names(i_field)), host_fields(i_field)%ptr, &
        solver, self%output_stride, data_loc, &
        self%field_buffers, self%last_shape_dims, self%last_stride_factors, &
        self%last_output_shape &
        )

      call writer_session%write_data( &
        trim(field_names(i_field)), &
        self%field_buffers(i_field)%buffer, &
        self%last_output_shape, &
        output_start, output_count, &
        self%convert_to_sp)
    end do
  end subroutine write_fields

  subroutine cleanup_output_buffers(self)
    !! Clean up dynamically allocated field buffers (internal).
    !!
    !! Frees memory allocated for field I/O buffers. Called during
    !! finalisation to prevent memory leaks.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance

    call cleanup_field_buffers(self%field_buffers)
  end subroutine cleanup_output_buffers

  subroutine finalise(self)
    !! Finalise snapshot manager and free resources.
    !!
    !! Cleans up all dynamically allocated buffers and closes the
    !! persistent snapshot file. Should be called at the end of
    !! simulation or when snapshot manager is no longer needed.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance

    call self%cleanup_output_buffers()
    call self%close_snapshot_file()
  end subroutine finalise

  subroutine open_snapshot_file(self, filename, comm)
    !! Open persistent snapshot file for appending timesteps (internal).
    !!
    !! Opens or creates a snapshot file that remains open across multiple
    !! snapshot writes. Each snapshot is written as a new timestep within
    !! the same file, enabling efficient time-series visualization.
    !!
    !! **Persistent File Strategy:**
    !! - File opened once at first snapshot
    !! - Remains open for subsequent snapshots (append mode)
    !! - Each write adds a new timestep to the file
    !! - Closed only during finalisation
    !!
    !! **Benefits:** Reduces file open/close overhead and keeps all snapshots
    !! in a single file for easy ParaView animation.
    !!
    !! **ADIOS2 Behaviour:** Automatically handles both creating new files
    !! and appending to existing ones based on file existence.
    class(snapshot_manager_t), intent(inout) :: self !! Snapshot manager instance
    character(len=*), intent(in) :: filename         !! Snapshot file path
    integer, intent(in) :: comm                      !! MPI communicator

    logical :: file_exists
    integer :: myrank, ierr

    call MPI_Comm_rank(comm, myrank, ierr)

    if (myrank == 0) then
      inquire (file=trim(filename), exist=file_exists)
      if (file_exists) then
        print *, 'Appending to existing snapshot file: ', trim(filename)
      else
        print *, 'Creating new snapshot file: ', trim(filename)
      end if
    end if

    call self%snapshot_writer%open(filename, comm)

    self%is_snapshot_file_open = .true.
  end subroutine open_snapshot_file

  subroutine close_snapshot_file(self)
    !! Close the persistent snapshot file
    class(snapshot_manager_t), intent(inout) :: self

    if (self%is_snapshot_file_open) then
      call self%snapshot_writer%close()
      self%is_snapshot_file_open = .false.
    end if
  end subroutine close_snapshot_file

end module m_snapshot_manager
