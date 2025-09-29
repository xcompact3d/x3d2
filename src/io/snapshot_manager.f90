module m_snapshot_manager
!! Snapshot manager for visualization output
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank
  use m_common, only: dp, i8, DIR_C, VERT, get_argument
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: allocate_io_writer, io_session_t
  use m_io_base, only: io_writer_t, io_file_t, io_mode_write
  use m_config, only: checkpoint_config_t
  use m_io_field_utils, only: field_buffer_map_t, field_ptr_t, &
                              setup_field_arrays, cleanup_field_arrays, &
                              stride_data_to_buffer, get_output_dimensions

  implicit none

  private
  public :: snapshot_manager_t

  type :: snapshot_manager_t
    class(io_writer_t), pointer :: writer => null()
    type(checkpoint_config_t) :: config
    integer, dimension(3) :: output_stride = [1, 1, 1]
    type(field_buffer_map_t), allocatable :: field_buffers(:)
    integer(i8), dimension(3) :: last_shape_dims = 0
    integer, dimension(3) :: last_stride_factors = 0
    integer(i8), dimension(3) :: last_output_shape = 0
    character(len=4096) :: vtk_xml = ""
  contains
    procedure :: init
    procedure :: handle_snapshot_step
    procedure :: finalise
    procedure, private :: write_snapshot
    procedure, private :: write_fields
    procedure, private :: cleanup_output_buffers
    procedure, private :: generate_vtk_xml
  end type snapshot_manager_t

contains

  subroutine init(self, comm)
    !! Initialize snapshot manager
    class(snapshot_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    call allocate_io_writer(self%writer)
    call self%writer%init(comm, "snapshot_writer")

    self%config = checkpoint_config_t()
    call self%config%read(nml_file=get_argument(1))

    if (self%config%snapshot_freq > 0) then
      call configure_output(self, comm)
    end if
  end subroutine init

  subroutine configure_output(self, comm)
    !! Configure snapshot output settings
    class(snapshot_manager_t), intent(inout) :: self
    integer, intent(in) :: comm

    integer :: myrank, ierr

    call MPI_Comm_rank(comm, myrank, ierr)

    self%output_stride = self%config%output_stride

    if (myrank == 0) then
      print *, 'Snapshot frequency: ', self%config%snapshot_freq
      print *, 'Snapshot prefix: ', trim(self%config%snapshot_prefix)
      print *, 'Output stride: ', self%output_stride
    end if
  end subroutine configure_output

  subroutine handle_snapshot_step(self, solver, timestep, comm)
    !! Handle snapshot writing at a given timestep
    class(snapshot_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in), optional :: comm

    integer :: comm_to_use

    comm_to_use = MPI_COMM_WORLD
    if (present(comm)) comm_to_use = comm

    call self%write_snapshot(solver, timestep, comm_to_use)
  end subroutine handle_snapshot_step

  subroutine write_snapshot(self, solver, timestep, comm)
    !! Write a snapshot file for visualization
    class(snapshot_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in) :: comm

    character(len=*), parameter :: field_names(*) = ["u", "v", "w"]
    integer :: myrank, ierr
    character(len=256) :: filename
    integer(i8), dimension(3) :: output_shape_dims
    integer, dimension(3) :: global_dims, output_dims
    type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
    real(dp), dimension(3) :: origin, original_spacing, output_spacing
    real(dp) :: simulation_time
    logical :: snapshot_uses_stride = .true.
    type(io_session_t) :: io_session
    integer :: i

    if (self%config%snapshot_freq <= 0) return
    if (mod(timestep, self%config%snapshot_freq) /= 0) return

    call MPI_Comm_rank(comm, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%config%snapshot_prefix), '_', timestep, '.bp'

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

    call io_session%open(filename, comm, io_mode_write)
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
      solver, io_session, solver%u%data_loc &
      )

    call io_session%close()

    call cleanup_field_arrays(solver, field_ptrs, host_fields)
  end subroutine write_snapshot

  subroutine generate_vtk_xml(self, dims, fields, origin, spacing)
    !! Generate VTK XML string for ImageData format for ParaView's ADIOS2VTXReader
    class(snapshot_manager_t), intent(inout) :: self
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

  subroutine write_fields( &
    self, field_names, host_fields, solver, io_session, data_loc &
    )
    !! Write field data with striding for snapshots
    class(snapshot_manager_t), intent(inout) :: self
    character(len=*), dimension(:), intent(in) :: field_names
    class(field_ptr_t), dimension(:), target, intent(in) :: host_fields
    class(solver_t), intent(in) :: solver
    type(io_session_t), intent(inout) :: io_session
    integer, intent(in) :: data_loc

    integer :: i_field

    ! Prepare buffers with striding for snapshots
    call prepare_field_buffers(solver, self%output_stride, field_names, data_loc)

    do i_field = 1, size(field_names)
      call write_single_field( &
        trim(field_names(i_field)), host_fields(i_field)%ptr, data_loc)
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

      integer, dimension(3) :: output_dims_local
      integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
      integer(i8), dimension(3) :: output_shape, output_start, output_count
      integer :: dims(3), buffer_idx
      logical :: buffer_found

      dims = solver%mesh%get_dims(data_loc)
      shape_dims = int(solver%mesh%get_global_dims(data_loc), i8)
      start_dims = int(solver%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)

      call get_output_dimensions( &
        shape_dims, start_dims, count_dims, self%output_stride, &
        output_shape, output_start, output_count, &
        output_dims_local, &
        self%last_shape_dims, self%last_stride_factors, self%last_output_shape &
        )

      ! Find the matching buffer for this field
      buffer_found = .false.
      do buffer_idx = 1, size(self%field_buffers)
        if (trim(self%field_buffers(buffer_idx)%field_name) == trim(field_name)) then
          buffer_found = .true.
          exit
        end if
      end do

      if (buffer_found) then
        call stride_data_to_buffer( &
          host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), dims, &
          self%output_stride, self%field_buffers(buffer_idx)%buffer, &
          output_dims_local &
          )

        call io_session%write_data( &
          field_name, self%field_buffers(buffer_idx)%buffer, &
          start_dims=output_start, count_dims=output_count &
          )
      else
        print *, 'INTERNAL ERROR: No buffer found for field: ', trim(field_name)
        error stop 'Missing field buffer'
      end if
    end subroutine write_single_field

  end subroutine write_fields

  subroutine cleanup_output_buffers(self)
    !! Clean up dynamic field buffers
    class(snapshot_manager_t), intent(inout) :: self
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
    !! Clean up snapshot manager
    class(snapshot_manager_t), intent(inout) :: self

    call self%cleanup_output_buffers()
    call self%writer%finalise()

    if (associated(self%writer)) deallocate(self%writer)
  end subroutine finalise

end module m_snapshot_manager