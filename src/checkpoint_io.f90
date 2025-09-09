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
!! Implementation of checkpoint manager when ADIOS2 is enabled
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Abort
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
    real(dp), dimension(:, :, :), allocatable :: strided_buffer
    real(dp), dimension(:, :, :), allocatable :: coords_x, coords_y, coords_z
    integer(i8), dimension(3) :: last_shape_dims = 0
    integer(i8), dimension(3) :: last_strided_shape = 0
    character(len=4096) :: vtk_xml = ""                 !! VTK XML string for ParaView compatibility
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
    procedure, private :: stride_data_to_buffer
    procedure, private :: cleanup_strided_buffers
    procedure, private :: get_strided_dimensions
    procedure, private :: generate_coordinates
    procedure, private :: generate_vtk_xml
  end type checkpoint_manager_adios2_t

  type :: field_ptr_t
    class(field_t), pointer :: ptr => null()
  end type field_ptr_t

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
      self%checkpoint_cfg%output_stride, &
      comm &
      )
  end subroutine init

  subroutine handle_restart(self, solver, comm)
    !! Check if a restart is needed and handle it
    class(checkpoint_manager_adios2_t), intent(inout) :: self
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

  subroutine configure( &
    self, checkpoint_freq, snapshot_freq, keep_checkpoint, &
    checkpoint_prefix, snapshot_prefix, output_stride, comm &
    )
    !! Configure checkpoint and snapshot settings
    class(checkpoint_manager_adios2_t), intent(inout) :: self
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
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=256) :: filename, temp_filename, old_filename
    type(adios2_file_t) :: file
    integer :: ierr, myrank
    integer :: comm_to_use
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    real(dp) :: simulation_time
    logical :: file_exists
    type(field_ptr_t), allocatable :: field_ptrs(:)

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

    file = self%adios2_writer%open(temp_filename, adios2_mode_write, &
                                   comm_to_use)
    call self%adios2_writer%begin_step(file)

    simulation_time = timestep*solver%dt
    call self%adios2_writer%write_data("timestep", timestep, file)
    call self%adios2_writer%write_data("time", real(simulation_time, dp), file)
    call self%adios2_writer%write_data("dt", real(solver%dt, dp), file)

    allocate (field_ptrs(3))
    field_ptrs(1)%ptr => solver%u
    field_ptrs(2)%ptr => solver%v
    field_ptrs(3)%ptr => solver%w

    call self%write_fields( &
      field_names, field_ptrs, &
      solver, file, use_stride=.false. &
      )
    call self%adios2_writer%close(file)

    deallocate (field_ptrs)

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
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    integer :: myrank, ierr, comm_to_use, i
    character(len=256) :: filename
    type(adios2_file_t) :: file
    integer(i8), dimension(3) :: strided_shape_dims
    integer, dimension(3) :: global_dims, strided_dims
    type(field_ptr_t), allocatable :: field_ptrs(:)
    real(dp), dimension(3) :: origin, original_spacing, strided_spacing, &
                              coords_max
    real(dp) :: simulation_time

    if (self%checkpoint_cfg%snapshot_freq <= 0) return
    if (mod(timestep, self%checkpoint_cfg%snapshot_freq) /= 0) return

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm
    call MPI_Comm_rank(comm_to_use, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%checkpoint_cfg%snapshot_prefix), '_', timestep, '.bp'
    if (myrank == 0) print *, 'Writing snapshot: ', trim(filename)

    global_dims = solver%mesh%get_global_dims(VERT)
    origin = solver%mesh%get_coordinates(1, 1, 1)
    coords_max = solver%mesh%geo%L
    original_spacing = solver%mesh%geo%d
    strided_spacing = original_spacing*real(self%output_stride, dp)

    do i = 1, 3
      strided_dims(i) = (global_dims(i) + self%output_stride(i) - 1)/ &
                        self%output_stride(i)
    end do
    strided_shape_dims = int(strided_dims, i8)

    file = self%adios2_writer%open(filename, adios2_mode_write, comm_to_use)
    call self%adios2_writer%begin_step(file)

    simulation_time = timestep*solver%dt
    call self%adios2_writer%write_data("time", real(simulation_time, dp), file)

    call self%generate_vtk_xml( &
      strided_shape_dims, field_names, origin, strided_spacing &
      )

    if (myrank == 0) then
      call self%adios2_writer%write_attribute("vtk.xml", self%vtk_xml, file)
    end if

    allocate (field_ptrs(3))
    field_ptrs(1)%ptr => solver%u
    field_ptrs(2)%ptr => solver%v
    field_ptrs(3)%ptr => solver%w

    call self%write_fields( &
      field_names, field_ptrs, &
      solver, file, use_stride=.true. &
      )
    deallocate (field_ptrs)

    call self%adios2_writer%close(file)
  end subroutine write_snapshot

  subroutine generate_vtk_xml(self, dims, fields, origin, spacing)
    !! Generate VTK XML string for ImageData format for ParaView's ADIOS2VTXReader
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    integer(i8), dimension(3), intent(in) :: dims
    character(len=*), dimension(:), intent(in) :: fields
    real(dp), dimension(3), intent(in) :: origin, spacing

    character(len=4096) :: xml
    character(len=96) :: extent_str, origin_str, spacing_str
    integer :: i

    ! for ImageData, the Extent defines the grid size (from 0 to N-1).
    write (extent_str, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
      '0 ', dims(1) - 1, ' 0 ', dims(2) - 1, ' 0 ', dims(3) - 1

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
      xml = trim(xml)//'        <DataArray Name="'//trim(fields(i))// &
            '" />'//new_line('a')
    end do

    xml = trim(xml)//'        <DataArray Name="TIME">time</DataArray>' &
          //new_line('a')

    xml = trim(xml)//'      </PointData>'//new_line('a')// &
          '    </Piece>'//new_line('a')// &
          '  </ImageData>'//new_line('a')// &
          '</VTKFile>'

    self%vtk_xml = xml
  end subroutine generate_vtk_xml

  subroutine generate_coordinates( &
    self, solver, file, shape_dims, start_dims, count_dims, data_loc &
    )
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    type(adios2_file_t), intent(inout) :: file
    integer(i8), dimension(3), intent(in) :: shape_dims, start_dims, count_dims
    integer, intent(in) :: data_loc

    integer :: i, nx, ny, nz
    real(dp), dimension(3) :: coords
    integer(i8), dimension(3) :: x_shape, y_shape, z_shape
    integer(i8), dimension(3) :: x_start, y_start, z_start
    integer(i8), dimension(3) :: x_count, y_count, z_count

    nx = int(count_dims(1))
    ny = int(count_dims(2))
    nz = int(count_dims(3))

    ! coordinates are structured as 3D arrays for ParaView ADIOS2 reader compatibility
    if (.not. allocated(self%coords_x) .or. size(self%coords_x) /= nx) then
      if (allocated(self%coords_x)) deallocate (self%coords_x)
      allocate (self%coords_x(nx, 1, 1))
    end if

    if (.not. allocated(self%coords_y) .or. size(self%coords_y) /= ny) then
      if (allocated(self%coords_y)) deallocate (self%coords_y)
      allocate (self%coords_y(1, ny, 1))
    end if

    if (.not. allocated(self%coords_z) .or. size(self%coords_z) /= nz) then
      if (allocated(self%coords_z)) deallocate (self%coords_z)
      allocate (self%coords_z(1, 1, nz))
    end if

    do i = 1, nx
      coords = solver%mesh%get_coordinates(i, 1, 1, data_loc)
      self%coords_x(i, 1, 1) = coords(1)
    end do

    do i = 1, ny
      coords = solver%mesh%get_coordinates(1, i, 1, data_loc)
      self%coords_y(1, i, 1) = coords(2)
    end do

    do i = 1, nz
      coords = solver%mesh%get_coordinates(1, 1, i, data_loc)
      self%coords_z(1, 1, i) = coords(3)
    end do

    x_shape = [1_i8, 1_i8, shape_dims(1)]
    x_start = [0_i8, 0_i8, start_dims(1)]
    x_count = [1_i8, 1_i8, count_dims(1)]

    y_shape = [1_i8, shape_dims(2), 1_i8]
    y_start = [0_i8, start_dims(2), 0_i8]
    y_count = [1_i8, count_dims(2), 1_i8]

    z_shape = [shape_dims(3), 1_i8, 1_i8]
    z_start = [start_dims(3), 0_i8, 0_i8]
    z_count = [count_dims(3), 1_i8, 1_i8]

    call self%adios2_writer%write_data( &
      "coordinates/x", self%coords_x, file, x_shape, x_start, x_count &
      )
    call self%adios2_writer%write_data( &
      "coordinates/y", self%coords_y, file, y_shape, y_start, y_count &
      )
    call self%adios2_writer%write_data( &
      "coordinates/z", self%coords_z, file, z_shape, z_start, z_count &
      )
  end subroutine generate_coordinates

  function stride_data( &
    self, input_data, dims, stride, strided_dims_out) &
    result(strided_data)
    !! stride the input data based on the specified stride
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(in) :: input_data
    integer, dimension(3), intent(in) :: dims
    integer, dimension(3), intent(in) :: stride
    integer, dimension(3), intent(out) :: strided_dims_out

    real(dp), dimension(:, :, :), allocatable :: strided_data
    integer :: i_stride, j_stride, k_stride
    integer :: i_max, j_max, k_max

    if (all(stride == 1)) then
      allocate (strided_data(dims(1), dims(2), dims(3)))
      strided_data = input_data
      strided_dims_out = dims
      return
    end if

    i_stride = stride(1); j_stride = stride(2); k_stride = stride(3)

    i_max = (dims(1) - 1)/i_stride + 1
    j_max = (dims(2) - 1)/j_stride + 1
    k_max = (dims(3) - 1)/k_stride + 1

    strided_dims_out = [i_max, j_max, k_max]
    allocate (strided_data(i_max, j_max, k_max))

    strided_data = input_data(1:dims(1):i_stride, &
                              1:dims(2):j_stride, 1:dims(3):k_stride)
  end function stride_data

  subroutine stride_data_to_buffer( &
    self, input_data, dims, stride, out_buffer, strided_dims_out &
    )
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(in) :: input_data
    integer, dimension(3), intent(in) :: dims
    integer, dimension(3), intent(in) :: stride
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: out_buffer
    integer, dimension(3), intent(out) :: strided_dims_out

    integer :: i_stride, j_stride, k_stride
    integer :: i_max, j_max, k_max

    if (all(stride == 1)) then
      if (allocated(out_buffer)) then
        if (size(out_buffer, 1) /= dims(1) &
            .or. size(out_buffer, 2) /= dims(2) .or. &
            size(out_buffer, 3) /= dims(3)) then
          deallocate (out_buffer)
          allocate (out_buffer(dims(1), dims(2), dims(3)))
        end if
      else
        allocate (out_buffer(dims(1), dims(2), dims(3)))
      end if
      out_buffer = input_data
      strided_dims_out = dims
      return
    end if

    i_stride = stride(1); j_stride = stride(2); k_stride = stride(3)

    i_max = (dims(1) + i_stride - 1)/i_stride
    j_max = (dims(2) + j_stride - 1)/j_stride
    k_max = (dims(3) + k_stride - 1)/k_stride

    strided_dims_out = [i_max, j_max, k_max]

    if (allocated(out_buffer)) then
      if (size(out_buffer, 1) /= i_max &
          .or. size(out_buffer, 2) /= j_max .or. &
          size(out_buffer, 3) /= k_max) then
        deallocate (out_buffer)
        allocate (out_buffer(i_max, j_max, k_max))
      end if
    else
      allocate (out_buffer(i_max, j_max, k_max))
    end if

    out_buffer = input_data(1:dims(1):i_stride, &
                            1:dims(2):j_stride, 1:dims(3):k_stride)
  end subroutine stride_data_to_buffer

  subroutine get_strided_dimensions( &
    self, shape_dims, start_dims, count_dims, stride_factors, &
    strided_shape, strided_start, strided_count, strided_dims_local)

    class(checkpoint_manager_adios2_t), intent(inout) :: self
    integer(i8), dimension(3), intent(in) :: shape_dims, start_dims, count_dims
    integer, dimension(3), intent(in) :: stride_factors
    integer(i8), dimension(3), intent(out) :: strided_shape, strided_start
    integer(i8), dimension(3), intent(out) :: strided_count
    integer, dimension(3), intent(out) :: strided_dims_local

    if (all(shape_dims == self%last_shape_dims) .and. &
        all(self%last_strided_shape > 0)) then
      strided_shape = self%last_strided_shape
    else
      strided_shape = [(shape_dims(3) + stride_factors(3) - 1_i8) &
                       /int(stride_factors(3), i8), &
                       (shape_dims(2) + stride_factors(2) - 1_i8) &
                       /int(stride_factors(2), i8), &
                       (shape_dims(1) + stride_factors(1) - 1_i8) &
                       /int(stride_factors(1), i8)]

      self%last_shape_dims = shape_dims
      self%last_strided_shape = strided_shape
    end if

    strided_start = [start_dims(3)/int(stride_factors(3), i8), &
                     start_dims(2)/int(stride_factors(2), i8), &
                     start_dims(1)/int(stride_factors(1), i8)]

    strided_dims_local = [(int(count_dims(1)) + stride_factors(1) - 1) &
                          /stride_factors(1), &
                          (int(count_dims(2)) + stride_factors(2) - 1) &
                          /stride_factors(2), &
                          (int(count_dims(3)) + stride_factors(3) - 1) &
                          /stride_factors(3)]

    strided_count = [int(strided_dims_local(3), i8), &
                     int(strided_dims_local(2), i8), &
                     int(strided_dims_local(1), i8)]
  end subroutine get_strided_dimensions

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
    real(dp), allocatable :: field_data(:, :, :)
    integer :: i, ierr
    integer :: dims(3)
    class(field_t), pointer :: host_field
    integer(i8), dimension(3) :: start_dims, count_dims
    character(len=*), parameter :: field_names(3) = ["u", "v", "w"]
    logical :: file_exists

    inquire (file=filename, exist=file_exists)
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

    allocate (field_data(dims(1), dims(2), dims(3)))

    do i = 1, size(field_names)
      call reader%read_data(field_names(i), &
                            field_data, file, start_dims, count_dims)

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
    if (allocated(field_data)) deallocate (field_data)
  end subroutine restart_checkpoint

  subroutine write_fields(self, field_names, fields, solver, file, use_stride)
    !! Write field data, optionally with striding
    class(checkpoint_manager_adios2_t), intent(inout) :: self
    character(len=*), dimension(:), intent(in) :: field_names
    class(field_ptr_t), dimension(:), target, intent(in) :: fields
    class(solver_t), intent(in) :: solver
    type(adios2_file_t), intent(inout) :: file
    logical, intent(in), optional :: use_stride

    class(field_t), pointer :: host_field
    integer :: i_field
    logical :: apply_stride

    call self%cleanup_strided_buffers()

    apply_stride = .false.
    if (present(use_stride)) apply_stride = use_stride

    do i_field = 1, size(field_names)
      host_field => solver%host_allocator%get_block( &
                    DIR_C, fields(i_field)%ptr%data_loc &
                    )
      call solver%backend%get_field_data(host_field%data, fields(i_field)%ptr)

      call write_single_field(trim(field_names(i_field)), host_field, &
                              fields(i_field)%ptr%data_loc)
      call solver%host_allocator%release_block(host_field)
    end do

  contains

    subroutine write_single_field(field_name, host_field, data_loc)
      character(len=*), intent(in) :: field_name
      class(field_t), pointer :: host_field
      integer, intent(in) :: data_loc
      integer, dimension(3) :: strided_dims_local, stride_factors

      integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
      integer(i8), dimension(3) :: strided_shape, strided_start, strided_count
      real(dp), dimension(:, :, :), pointer :: field_data
      integer :: dims(3)

      dims = solver%mesh%get_dims(data_loc)
      shape_dims = int(solver%mesh%get_global_dims(data_loc), i8)
      start_dims = int(solver%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)

      field_data => host_field%data(1:dims(1), 1:dims(2), 1:dims(3))

      if (use_stride) then
        stride_factors = self%output_stride

        call self%get_strided_dimensions( &
          shape_dims, start_dims, count_dims, stride_factors, &
          strided_shape, strided_start, strided_count, &
          strided_dims_local &
          )

        call self%stride_data_to_buffer( &
          field_data, dims, stride_factors, &
          self%strided_buffer, strided_dims_local)

        call self%adios2_writer%write_data( &
          field_name, self%strided_buffer, &
          file, strided_shape, strided_start, strided_count &
          )
      else
        call self%adios2_writer%write_data( &
          field_name, field_data, &
          file, shape_dims, start_dims, count_dims &
          )
      end if
    end subroutine write_single_field

  end subroutine write_fields

  subroutine cleanup_strided_buffers(self)
    class(checkpoint_manager_adios2_t), intent(inout) :: self

    if (allocated(self%strided_buffer)) deallocate (self%strided_buffer)
  end subroutine cleanup_strided_buffers

  subroutine finalise(self)
    class(checkpoint_manager_adios2_t), intent(inout) :: self

    call self%cleanup_strided_buffers()
    call self%adios2_writer%finalise()
  end subroutine finalise

end module m_checkpoint_manager_impl

module m_checkpoint_manager
  !! Public facade and factory function for checkpoint manager
  use m_checkpoint_manager_impl, only: checkpoint_manager_adios2_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: checkpoint_manager_t, create_checkpoint_manager

  type :: checkpoint_manager_t
    type(checkpoint_manager_adios2_t) :: impl
  contains
    procedure :: init => cm_init
    procedure :: handle_restart => cm_handle_restart
    procedure :: handle_io_step => cm_handle_io_step
    procedure :: finalise => cm_finalise
    procedure :: is_restart => cm_is_restart
  end type checkpoint_manager_t

  interface checkpoint_manager_t
    module procedure create_checkpoint_manager
  end interface checkpoint_manager_t

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
