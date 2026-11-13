module m_inflow_planes
  !! Time-indexed y-z velocity planes for precursor-successor inflow.
  !!
  !! A precursor run samples u, v, w on one y-z plane every `sample_every`
  !! iterations and writes them in blocks of `ntimesteps` planes, one ADIOS2
  !! file per block (Incompact3d `append_outflow`/`write_outflow`). A successor
  !! reads the blocks back with `read_inflow_block` and imposes them as a
  !! time-varying inflow. The on-disk layout is documented under "Precursor
  !! inflow planes" in docs/source/user/input_file.rst; bump
  !! INFLOW_FORMAT_VERSION when it changes.
  use m_mpi, only: MPI_COMM_WORLD

  use m_base_backend, only: base_backend_t
  use m_common, only: dp, i8, DIR_C, VERT
  use m_config, only: inflow_config_t
  use m_field, only: field_t
  use m_io_session, only: writer_session_t, reader_session_t
  use m_logger, only: log_info, log_warn, abort_run, str
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: inflow_plane_writer_t, inflow_plane_reader_t, inflow_file_name, &
            read_inflow_block, INFLOW_FORMAT_VERSION

  integer, parameter :: INFLOW_FORMAT_VERSION = 1

  type :: inflow_plane_writer_t
    class(base_backend_t), pointer :: backend => null()
    type(mesh_t), pointer :: mesh => null()
    type(inflow_config_t) :: cfg
    !> Buffered planes, (ntimesteps, ny, nz) per component on local y-z
    real(dp), allocatable :: u(:, :, :), v(:, :, :), w(:, :, :)
    !> Simulation time of each plane; 3D only to fit the array writer
    real(dp), allocatable :: times(:, :, :)
    !> Padded Cartesian scratch for the field-to-host copy
    real(dp), allocatable :: host(:, :, :)
    real(dp) :: dt = 0._dp       !! solver time step
    integer :: x_index = 0       !! global x vertex index of the plane
    integer :: i_local = 0       !! local x index, 0 if not on this rank
    integer :: n_filled = 0      !! planes sampled into the current block
  contains
    procedure :: sample
    procedure, private :: write_block
  end type inflow_plane_writer_t

  interface inflow_plane_writer_t
    module procedure init
  end interface inflow_plane_writer_t

  type :: plane_block_t
    integer :: id = 0            !! file index held, 0 if empty
    real(dp), allocatable :: u(:, :, :), v(:, :, :), w(:, :, :)
  end type plane_block_t

  type :: inflow_plane_reader_t
    !! Replays the planes of `ninflows` files as a time-varying inflow.
    !! Plane n (1-based, counted across the files) is the inflow at time
    !! n*plane_dt; between planes it is interpolated linearly, and after
    !! the last plane the series repeats from the first.
    type(mesh_t), pointer :: mesh => null()
    type(inflow_config_t) :: cfg
    integer :: ntimesteps = 0    !! planes per file
    real(dp) :: plane_dt = 0._dp !! time between planes
    !> Two files are enough to interpolate across a file boundary
    type(plane_block_t) :: blocks(2)
  contains
    procedure :: planes_at
    procedure, private :: block_slot
  end type inflow_plane_reader_t

  interface inflow_plane_reader_t
    module procedure init_reader
  end interface inflow_plane_reader_t

contains

  function init(backend, mesh, cfg, dt) result(writer)
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(inflow_config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt
    type(inflow_plane_writer_t) :: writer

    integer :: dims(3), global_dims(3), padded_dims(3), nt

    writer%backend => backend
    writer%mesh => mesh
    writer%cfg = cfg
    writer%dt = dt

    if (.not. cfg%write_planes) return

#ifndef WITH_ADIOS2
    call abort_run('inflow_nml: write_planes requires a build WITH_ADIOS2.', &
                   collective=.true.)
#endif

    dims = mesh%get_dims(VERT)
    global_dims = mesh%get_global_dims(VERT)

    writer%x_index = cfg%plane_x_index
    if (writer%x_index == 0) writer%x_index = global_dims(1)
    if (writer%x_index > global_dims(1)) &
      call abort_run('inflow_nml: plane_x_index exceeds the x vertex &
                     &count.', collective=.true.)

    writer%i_local = writer%x_index - mesh%par%n_offset(1)
    if (writer%i_local < 1 .or. writer%i_local > dims(1)) writer%i_local = 0

    nt = cfg%ntimesteps
    allocate (writer%u(nt, dims(2), dims(3)), writer%v(nt, dims(2), dims(3)), &
              writer%w(nt, dims(2), dims(3)), source=0._dp)
    allocate (writer%times(nt, 1, 1), source=0._dp)
    if (writer%i_local > 0) then
      padded_dims = backend%allocator%get_padded_dims(DIR_C)
      allocate (writer%host(padded_dims(1), padded_dims(2), padded_dims(3)))
    end if
  end function init

  pure function inflow_file_name(prefix, block) result(filename)
    !! File holding block `block` (1-based) of a plane database
    character(*), intent(in) :: prefix
    integer, intent(in) :: block
    character(:), allocatable :: filename

    character(len=20) :: block_str

    write (block_str, '(I0)') block
    filename = trim(prefix)//'_'//trim(block_str)//'.bp'
  end function inflow_file_name

  subroutine sample(self, iter, t, u, v, w)
    !! Record the plane at iteration `iter` if it is a sampling step, and
    !! write the block once its last plane is in. A block only partly sampled
    !! by this run (a restart inside it) is skipped rather than written with
    !! stale planes.
    class(inflow_plane_writer_t), intent(inout) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t
    class(field_t), intent(in) :: u, v, w

    integer :: isample, slot, block, nt

    if (.not. self%cfg%write_planes) return
    if (iter <= self%cfg%start_iter) return
    if (mod(iter - self%cfg%start_iter, self%cfg%sample_every) /= 0) return

    nt = self%cfg%ntimesteps
    isample = (iter - self%cfg%start_iter)/self%cfg%sample_every
    slot = mod(isample - 1, nt) + 1
    if (slot == 1) self%n_filled = 0

    if (self%i_local > 0) then
      call extract_plane(self, self%u, slot, u)
      call extract_plane(self, self%v, slot, v)
      call extract_plane(self, self%w, slot, w)
    end if
    self%times(slot, 1, 1) = t
    self%n_filled = self%n_filled + 1

    if (slot /= nt) return

    block = isample/nt
    if (self%n_filled == nt) then
      call self%write_block(block, iter)
    else
      call log_warn('inflow block '//str(block)//' was only partly sampled &
                    &by this run and is not written')
    end if
  end subroutine sample

  subroutine extract_plane(self, planes, slot, f)
    class(inflow_plane_writer_t), intent(inout) :: self
    real(dp), intent(inout) :: planes(:, :, :)
    integer, intent(in) :: slot
    class(field_t), intent(in) :: f

    integer :: dims(3)

    dims = self%mesh%get_dims(VERT)
    call self%backend%get_field_data(self%host, f)
    planes(slot, :, :) = self%host(self%i_local, 1:dims(2), 1:dims(3))
  end subroutine extract_plane

  subroutine write_block(self, block, iter)
    class(inflow_plane_writer_t), intent(inout) :: self
    integer, intent(in) :: block, iter

    type(writer_session_t) :: writer
    character(:), allocatable :: filename
    integer(i8) :: shape_dims(3), start_dims(3), count_dims(3)
    integer :: dims(3), global_dims(3)
    integer :: format_version, ntimesteps, sample_every, first_iter, x_index
    real(dp) :: dt, x

    dims = self%mesh%get_dims(VERT)
    global_dims = self%mesh%get_global_dims(VERT)

    ! Scalars are put deferred, so they need storage that outlives the puts
    format_version = INFLOW_FORMAT_VERSION
    ntimesteps = self%cfg%ntimesteps
    sample_every = self%cfg%sample_every
    first_iter = iter - (ntimesteps - 1)*sample_every
    x_index = self%x_index
    dt = self%dt
    x = self%mesh%geo%origin(1) + (x_index - 1)*self%mesh%geo%d(1)

    filename = inflow_file_name(self%cfg%inflow_prefix, block)
    call log_info('Writing inflow planes to '//filename)

    call writer%open(filename, MPI_COMM_WORLD)
    if (.not. writer%is_session_functional()) &
      call abort_run('Writing inflow planes needs a functional ADIOS2 &
                     &session.')

    if (self%mesh%par%is_root()) then
      call writer%write_data('format_version', format_version)
      call writer%write_data('ntimesteps', ntimesteps)
      call writer%write_data('sample_every', sample_every)
      call writer%write_data('block', block)
      call writer%write_data('first_iter', first_iter)
      call writer%write_data('x_index', x_index)
      call writer%write_data('x', x)
      call writer%write_data('dt', dt)
      call writer%write_data('ny', global_dims(2))
      call writer%write_data('nz', global_dims(3))
      call writer%write_data( &
        'time', self%times, [int(ntimesteps, i8), 1_i8, 1_i8], &
        [0_i8, 0_i8, 0_i8], [int(ntimesteps, i8), 1_i8, 1_i8])
    end if

    if (self%i_local > 0) then
      shape_dims = int([ntimesteps, global_dims(2), global_dims(3)], i8)
      start_dims = int([0, self%mesh%par%n_offset(2), &
                        self%mesh%par%n_offset(3)], i8)
      count_dims = int([ntimesteps, dims(2), dims(3)], i8)
      call writer%write_data('u', self%u, shape_dims, start_dims, count_dims)
      call writer%write_data('v', self%v, shape_dims, start_dims, count_dims)
      call writer%write_data('w', self%w, shape_dims, start_dims, count_dims)
    end if

    call writer%close()
  end subroutine write_block

  subroutine read_inflow_block(filename, mesh, u, v, w, times, plane_dt)
    !! Read one block on this rank's y-z subdomain. The arrays are allocated
    !! here as (ntimesteps, ny, nz) on the local VERT dims, and `times` as
    !! (ntimesteps). The rank layout may differ from the precursor's, but the
    !! global y-z vertex counts must match.
    character(*), intent(in) :: filename
    type(mesh_t), intent(in) :: mesh
    real(dp), allocatable, intent(out) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(dp), allocatable, intent(out) :: times(:)
    !> Time between planes, sample_every*dt of the precursor
    real(dp), optional, intent(out) :: plane_dt

    type(reader_session_t) :: reader
    integer(i8) :: start_dims(3), count_dims(3)
    integer :: dims(3), global_dims(3)
    integer :: format_version, ntimesteps, ny, nz, sample_every
    real(dp) :: dt
    real(dp), allocatable :: times_3d(:, :, :)

    dims = mesh%get_dims(VERT)
    global_dims = mesh%get_global_dims(VERT)

    call reader%open(filename, MPI_COMM_WORLD)
    call reader%read_data('format_version', format_version)
    if (format_version /= INFLOW_FORMAT_VERSION) &
      call abort_run('Inflow planes: unsupported format_version.', &
                     collective=.true.)
    call reader%read_data('ntimesteps', ntimesteps)
    call reader%read_data('ny', ny)
    call reader%read_data('nz', nz)
    if (ny /= global_dims(2) .or. nz /= global_dims(3)) &
      call abort_run('Inflow planes: y-z grid does not match the mesh.', &
                     collective=.true.)
    if (present(plane_dt)) then
      call reader%read_data('sample_every', sample_every)
      call reader%read_data('dt', dt)
      plane_dt = sample_every*dt
    end if

    allocate (u(ntimesteps, dims(2), dims(3)))
    allocate (v(ntimesteps, dims(2), dims(3)), w(ntimesteps, dims(2), dims(3)))
    allocate (times_3d(ntimesteps, 1, 1))

    call reader%read_data('time', times_3d, [0_i8, 0_i8, 0_i8], &
                          [int(ntimesteps, i8), 1_i8, 1_i8])

    start_dims = int([0, mesh%par%n_offset(2), mesh%par%n_offset(3)], i8)
    count_dims = int([ntimesteps, dims(2), dims(3)], i8)
    call reader%read_data('u', u, start_dims, count_dims)
    call reader%read_data('v', v, start_dims, count_dims)
    call reader%read_data('w', w, start_dims, count_dims)
    call reader%close()

    times = times_3d(:, 1, 1)
  end subroutine read_inflow_block

  function init_reader(mesh, cfg) result(reader)
    !! Opens the first file to learn the plane count and spacing
    type(mesh_t), target, intent(inout) :: mesh
    type(inflow_config_t), intent(in) :: cfg
    type(inflow_plane_reader_t) :: reader

    type(reader_session_t) :: session
    integer :: format_version, sample_every
    real(dp) :: dt

    reader%mesh => mesh
    reader%cfg = cfg
    if (.not. cfg%read_planes) return

#ifndef WITH_ADIOS2
    call abort_run('inflow_nml: read_planes requires a build WITH_ADIOS2.', &
                   collective=.true.)
#endif

    ! Only the metadata here; the planes are read on first use, so the
    ! constructed reader carries no plane arrays through the copies into
    ! its owner (nvfortran mis-copies them).
    call session%open(inflow_file_name(cfg%inflow_read_prefix, 1), &
                      MPI_COMM_WORLD)
    call session%read_data('format_version', format_version)
    if (format_version /= INFLOW_FORMAT_VERSION) &
      call abort_run('Inflow planes: unsupported format_version.', &
                     collective=.true.)
    call session%read_data('ntimesteps', reader%ntimesteps)
    call session%read_data('sample_every', sample_every)
    call session%read_data('dt', dt)
    call session%close()
    reader%plane_dt = sample_every*dt
    if (reader%plane_dt <= 0._dp) &
      call abort_run('Inflow planes: the plane spacing must be positive.', &
                     collective=.true.)

    call log_info('Inflow planes: '//str(cfg%ninflows)//' file(s) of ' &
                  //str(reader%ntimesteps)//' planes, plane spacing ' &
                  //str(reader%plane_dt, '(es12.5)'))
  end function init_reader

  subroutine planes_at(self, t, u, v, w)
    !! Inflow plane at time `t` on this rank's local (ny, nz)
    class(inflow_plane_reader_t), intent(inout) :: self
    real(dp), intent(in) :: t
    real(dp), intent(out) :: u(:, :), v(:, :), w(:, :)

    integer :: n_planes, n0, n1, slot0, slot1, i0, i1, j, k
    real(dp) :: position, frac

    n_planes = self%ntimesteps*self%cfg%ninflows

    ! Fractional plane number, wrapped into [1, n_planes + 1). Planes
    ! n_planes and 1 are neighbours when the series repeats.
    position = modulo(t/self%plane_dt - 1._dp, real(n_planes, dp)) + 1._dp
    n0 = min(int(position), n_planes)
    frac = position - n0
    n1 = mod(n0, n_planes) + 1

    slot0 = self%block_slot((n0 - 1)/self%ntimesteps + 1, 0)
    slot1 = self%block_slot((n1 - 1)/self%ntimesteps + 1, slot0)
    i0 = mod(n0 - 1, self%ntimesteps) + 1
    i1 = mod(n1 - 1, self%ntimesteps) + 1

    ! Plane i of each block is a strided section; loop rather than form
    ! array temporaries
    do k = 1, size(u, 2)
      do j = 1, size(u, 1)
        u(j, k) = (1._dp - frac)*self%blocks(slot0)%u(i0, j, k) &
                  + frac*self%blocks(slot1)%u(i1, j, k)
        v(j, k) = (1._dp - frac)*self%blocks(slot0)%v(i0, j, k) &
                  + frac*self%blocks(slot1)%v(i1, j, k)
        w(j, k) = (1._dp - frac)*self%blocks(slot0)%w(i0, j, k) &
                  + frac*self%blocks(slot1)%w(i1, j, k)
      end do
    end do
  end subroutine planes_at

  integer function block_slot(self, id, keep) result(slot)
    !! Slot holding file `id`, reading it into the slot other than `keep`
    !! if it is not loaded yet
    class(inflow_plane_reader_t), intent(inout) :: self
    integer, intent(in) :: id, keep

    real(dp), allocatable :: times(:)
    real(dp) :: plane_dt

    do slot = 1, 2
      if (self%blocks(slot)%id == id) return
    end do

    slot = 1
    if (keep == 1) slot = 2

    call read_inflow_block( &
      inflow_file_name(self%cfg%inflow_read_prefix, id), self%mesh, &
      self%blocks(slot)%u, self%blocks(slot)%v, self%blocks(slot)%w, &
      times, plane_dt)
    self%blocks(slot)%id = id

    if (size(times) /= self%ntimesteps .or. &
        abs(plane_dt - self%plane_dt) > 1e-12_dp*self%plane_dt) then
      call abort_run('Inflow planes: all files need the same plane count &
                     &and spacing.', collective=.true.)
    end if
  end function block_slot

end module m_inflow_planes
