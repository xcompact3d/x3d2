module m_mesh
  !! Mesh module providing high-level mesh management and query functions.
  !!
  !! This module defines the `mesh_t` type which aggregates geometry, grid, and
  !! parallel decomposition information. It provides methods to query mesh
  !! dimensions, coordinates, and other mesh properties for both global and
  !! local (per MPI rank) domains.

  use iso_fortran_env, only: stderr => error_unit

  use mpi
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      CELL, VERT, X_FACE, Y_FACE, Z_FACE, &
                      X_EDGE, Y_EDGE, Z_EDGE, &
                      BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, BC_HALO
  use m_field, only: field_t
  use m_mesh_content

  implicit none

  type :: mesh_t
    !! Mesh type containing all mesh information for the simulation.
    !!
    !! This type aggregates three main components:
    !! - geo: Geometry information (coordinates, stretching)
    !! - grid: Grid dimensions and boundary conditions
    !! - par: Parallel domain decomposition information
    !!
    !! The mesh is initialised once and should be treated as read-only
    !! during the simulation.
    type(geo_t), allocatable :: geo     !! Geometry information
    class(grid_t), allocatable :: grid  !! Grid dimensions and boundary conditions
    class(par_t), allocatable :: par    !! Parallel decomposition information
  contains
    procedure :: get_dims         !! Get local dimensions for a data location
    procedure :: get_global_dims  !! Get global dimensions for a data location

    procedure :: get_n_dir        !! Get number of grid points in a direction
    procedure :: get_n_phi        !! Get number of grid points for a field
    generic :: get_n => get_n_dir, get_n_phi  !! Generic interface for get_n

    procedure :: get_coordinates  !! Get coordinate array for a direction
  end type mesh_t

  interface mesh_t
    module procedure mesh_init
  end interface mesh_t

contains

  function mesh_init(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                     stretching, beta, use_2decomp) result(mesh)
    !! Initialise the mesh object with global domain parameters.
    !!
    !! Creates and fully initialises a mesh object containing geometry, grid, and
    !! parallel decomposition information. The mesh should be treated as read-only
    !! after initialisation. Supports both uniform and stretched meshes, and can
    !! use either 2decomp or generic domain decomposition.
    use m_decomp, only: is_avail_2decomp, decomposition_2decomp
    integer, dimension(3), intent(in) :: dims_global  !! Global grid dimensions [nx, ny, nz]
    integer, dimension(3), intent(in) :: nproc_dir    !! Number of processors in each direction
    real(dp), dimension(3), intent(in) :: L_global    !! Physical domain lengths [Lx, Ly, Lz]
    character(len=*), dimension(2), intent(in) :: BC_x  !! Boundary conditions in x (lower, upper)
    character(len=*), dimension(2), intent(in) :: BC_y  !! Boundary conditions in y (lower, upper)
    character(len=*), dimension(2), intent(in) :: BC_z  !! Boundary conditions in z (lower, upper)
    character(len=*), dimension(3), optional, intent(in) :: stretching  !! Mesh stretching type per direction
    real(dp), dimension(3), optional, intent(in) :: beta  !! Stretching parameters per direction
    logical, optional, intent(in) :: use_2decomp      !! Flag to use 2decomp library
    class(mesh_t), allocatable :: mesh                !! Initialised mesh object

    character(len=20), dimension(3, 2) :: BC_all
    logical :: is_first_domain, is_last_domain
    integer :: dir, j
    integer :: ierr

    allocate (mesh)
    allocate (mesh%geo)
    allocate (mesh%grid)
    allocate (mesh%par)

    BC_all(1, 1) = BC_x(1); BC_all(1, 2) = BC_x(2)
    BC_all(2, 1) = BC_y(1); BC_all(2, 2) = BC_y(2)
    BC_all(3, 1) = BC_z(1); BC_all(3, 2) = BC_z(2)
    do dir = 1, 3
      do j = 1, 2
        select case (trim(BC_all(dir, j)))
        case ('periodic')
          mesh%grid%BCs_global(dir, j) = BC_PERIODIC
        case ('neumann')
          mesh%grid%BCs_global(dir, j) = BC_NEUMANN
        case ('dirichlet')
          mesh%grid%BCs_global(dir, j) = BC_DIRICHLET
        case default
          error stop 'Unknown BC'
        end select
      end do
    end do

    do dir = 1, 3
      if (any(mesh%grid%BCs_global(dir, :) == BC_PERIODIC) .and. &
          (.not. all(mesh%grid%BCs_global(dir, :) == BC_PERIODIC))) then
        error stop 'BCs are incompatible: in a direction make sure to have &
                    &either both sides periodic or none.'
      end if
      mesh%grid%periodic_BC(dir) = all(mesh%grid%BCs_global(dir, :) &
                                       == BC_PERIODIC)
    end do

    ! Set global vertex dims
    mesh%grid%global_vert_dims(:) = dims_global

    ! Set global cell dims
    do dir = 1, 3
      if (mesh%grid%periodic_BC(dir)) then
        mesh%grid%global_cell_dims(dir) = mesh%grid%global_vert_dims(dir)
      else
        mesh%grid%global_cell_dims(dir) = mesh%grid%global_vert_dims(dir) - 1
      end if
    end do

    ! Parallel domain decomposition
    mesh%par%nproc_dir(:) = nproc_dir
    mesh%par%nproc = product(nproc_dir(:))
    call MPI_Comm_rank(MPI_COMM_WORLD, mesh%par%nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mesh%par%nproc, ierr)

    ! Either use 2decomp or the generic decomposition
    if (present(use_2decomp)) then
      if (is_avail_2decomp() .and. use_2decomp) then
        call decomposition_2decomp(mesh%grid, mesh%par)
      else
        call decomposition_generic(mesh%grid, mesh%par)
      end if
    else
      call decomposition_generic(mesh%grid, mesh%par)
    end if

    ! Set subdomain BCs
    do dir = 1, 3
      is_first_domain = mesh%par%nrank_dir(dir) == 0
      is_last_domain = mesh%par%nrank_dir(dir) + 1 == mesh%par%nproc_dir(dir)
      ! subdomain-subdomain boundaries are identical to periodic BCs
      if (is_first_domain .and. is_last_domain) then
        mesh%grid%BCs(dir, 1) = mesh%grid%BCs_global(dir, 1)
        mesh%grid%BCs(dir, 2) = mesh%grid%BCs_global(dir, 2)
      else if (is_first_domain) then
        mesh%grid%BCs(dir, 1) = mesh%grid%BCs_global(dir, 1)
        mesh%grid%BCs(dir, 2) = BC_HALO
      else if (is_last_domain) then
        mesh%grid%BCs(dir, 1) = BC_HALO
        mesh%grid%BCs(dir, 2) = mesh%grid%BCs_global(dir, 2)
      else
        mesh%grid%BCs(dir, :) = BC_HALO
      end if
    end do

    ! Geometry
    mesh%geo%L = L_global
    mesh%geo%d = mesh%geo%L/mesh%grid%global_cell_dims

    if (present(stretching)) then
      mesh%geo%stretching = stretching
    else
      mesh%geo%stretching(:) = 'uniform'
    end if

    if (present(beta)) then
      mesh%geo%beta = beta
    else
      mesh%geo%beta(:) = 1
    end if

    call mesh%geo%obtain_coordinates( &
      mesh%grid%vert_dims, mesh%grid%cell_dims, mesh%par%n_offset &
      )

  end function mesh_init

  subroutine decomposition_generic(grid, par)
    ! Generic decomposition used when 2decomp isn't used

    use m_mesh_content, only: par_t, grid_t

    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par
    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer :: i, nproc_x, nproc_y, nproc_z

    if (par%is_root()) then
      print *, "Domain decomposition by x3d2 (generic)"
    end if

    ! Number of processes on a direction basis
    nproc_x = par%nproc_dir(1)
    nproc_y = par%nproc_dir(2)
    nproc_z = par%nproc_dir(3)

    ! Define number of cells and vertices in each direction
    grid%vert_dims = grid%global_vert_dims/par%nproc_dir

    ! A 3D array corresponding to each region in the global domain
    allocate (global_ranks(nproc_x, nproc_y, nproc_z))

    ! set the corresponding global rank for each sub-domain
    global_ranks = reshape([(i, i=0, par%nproc - 1)], &
                           shape=[nproc_x, nproc_y, nproc_z])

    call par%compute_rank_pos_from_global(global_ranks)
    call grid%copy_vert2cell_dims(par)

    par%n_offset(:) = grid%vert_dims(:)*par%nrank_dir(:)

  end subroutine

  pure function get_dims(self, data_loc) result(dims)
    !! Get local domain dimensions for a specific data location.
    !!
    !! Returns the dimensions of the local subdomain (on this MPI rank) for
    !! the specified data location (VERT, CELL, X_FACE, etc.).
    class(mesh_t), intent(in) :: self  !! Mesh object
    integer, intent(in) :: data_loc    !! Data location flag (VERT, CELL, etc.)
    integer, dimension(3) :: dims      !! Local dimensions [nx, ny, nz]

    dims = get_dims_dataloc(data_loc, self%grid%vert_dims, self%grid%cell_dims)
  end function

  pure function get_global_dims(self, data_loc) result(dims)
    !! Get global domain dimensions for a specific data location.
    !!
    !! Returns the dimensions of the entire global domain for the specified
    !! data location (VERT, CELL, X_FACE, etc.).
    class(mesh_t), intent(in) :: self  !! Mesh object
    integer, intent(in) :: data_loc    !! Data location flag (VERT, CELL, etc.)
    integer, dimension(3) :: dims      !! Global dimensions [nx, ny, nz]

    dims = get_dims_dataloc(data_loc, self%grid%global_vert_dims, &
                            self%grid%global_cell_dims)
  end function

  pure function get_dims_dataloc(data_loc, vert_dims, cell_dims) result(dims)
  !! Getter for domain dimensions
    integer, intent(in) :: data_loc
    integer, dimension(3), intent(in) :: vert_dims, cell_dims
    integer, dimension(3) :: dims

    select case (data_loc)
    case (VERT)
      dims = vert_dims
    case (CELL)
      dims = cell_dims
    case (X_FACE)
      dims(1) = vert_dims(1)
      dims(2:3) = cell_dims(2:3)
    case (Y_FACE)
      dims(1) = cell_dims(1)
      dims(2) = vert_dims(2)
      dims(3) = cell_dims(3)
    case (Z_FACE)
      dims(1:2) = cell_dims(1:2)
      dims(3) = vert_dims(3)
    case (X_EDGE)
      dims(1) = cell_dims(1)
      dims(2:3) = vert_dims(2:3)
    case (Y_EDGE)
      dims(1) = vert_dims(1)
      dims(2) = cell_dims(2)
      dims(3) = vert_dims(3)
    case (Z_EDGE)
      dims(1:2) = vert_dims(1:2)
      dims(3) = cell_dims(3)
    case default
      error stop "Unknown location in get_dims_dataloc"
    end select
  end function get_dims_dataloc

  pure function get_n_phi(self, phi) result(n)
    !! Get the main dimension (pencil length) for a field.
    !!
    !! Returns the number of grid points along the primary direction for the
    !! given field, accounting for both the field's orientation (dir) and
    !! data location on the staggered grid.
    class(mesh_t), intent(in) :: self   !! Mesh object
    class(field_t), intent(in) :: phi   !! Field to query
    integer :: n                        !! Number of grid points in main direction

    n = self%get_n(phi%dir, phi%data_loc)

  end function

  pure function get_n_dir(self, dir, data_loc) result(n)
    !! Get the main dimension for a field with given direction and data location.
    !!
    !! Returns the number of grid points along a specified direction for a field
    !! located at the given position on the staggered grid. Handles the different
    !! grid dimensions for vertex-centered vs cell-centered data.
    class(mesh_t), intent(in) :: self     !! Mesh object
    integer, intent(in) :: dir            !! Primary direction (DIR_X, DIR_Y, DIR_Z)
    integer, intent(in) :: data_loc       !! Data location (VERT, CELL, X_FACE, etc.)
    integer :: n                          !! Number of grid points in direction
    integer :: n_cell, n_vert

    n_cell = self%grid%cell_dims(dir)
    n_vert = self%grid%vert_dims(dir)

    ! default to n_vert
    n = n_vert

    select case (data_loc)
    case (CELL)
      n = n_cell
    case (VERT)
      n = n_vert
    case (X_FACE)
      if (dir /= DIR_X) then
        n = n_cell
      end if
    case (Y_FACE)
      if (dir /= DIR_Y) then
        n = n_cell
      end if
    case (Z_FACE)
      if (dir /= DIR_Z) then
        n = n_cell
      end if
    case (X_EDGE)
      if (dir == DIR_X) then
        n = n_cell
      end if
    case (Y_EDGE)
      if (dir == DIR_Y) then
        n = n_cell
      end if
    case (Z_EDGE)
      if (dir == DIR_Z) then
        n = n_cell
      end if
    case default
      error stop "Unknown direction in get_n_dir"
    end select
  end function get_n_dir

  pure function get_coordinates(self, i, j, k, data_loc_op) result(coords)
    !! Get physical coordinates for a grid point with given indices.
    !!
    !! Returns the physical (x, y, z) coordinates for a grid point specified by
    !! local Cartesian indices (i, j, k) at the given data location. Default
    !! location is vertex-centered (VERT). Note: Avoid calling this function in
    !! hot loops due to performance overhead.
    class(mesh_t), intent(in) :: self             !! Mesh object
    integer, intent(in) :: i, j, k                !! Local Cartesian indices
    integer, optional, intent(in) :: data_loc_op  !! Data location (default: VERT)
    integer :: data_loc
    real(dp), dimension(3) :: coords              !! Physical coordinates [x, y, z]

    if (present(data_loc_op)) then
      data_loc = data_loc_op
    else
      data_loc = VERT
    end if

    select case (data_loc)
    case (VERT)
      coords(1) = self%geo%vert_coords(i, 1)
      coords(2) = self%geo%vert_coords(j, 2)
      coords(3) = self%geo%vert_coords(k, 3)
    case (CELL)
      coords(1) = self%geo%midp_coords(i, 1)
      coords(2) = self%geo%midp_coords(j, 2)
      coords(3) = self%geo%midp_coords(k, 3)
    case (X_FACE)
      coords(1) = self%geo%vert_coords(i, 1)
      coords(2) = self%geo%midp_coords(j, 2)
      coords(3) = self%geo%midp_coords(k, 3)
    case (Y_FACE)
      coords(1) = self%geo%midp_coords(i, 1)
      coords(2) = self%geo%vert_coords(j, 2)
      coords(3) = self%geo%midp_coords(k, 3)
    case (Z_FACE)
      coords(1) = self%geo%midp_coords(i, 1)
      coords(2) = self%geo%midp_coords(j, 2)
      coords(3) = self%geo%vert_coords(k, 3)
    case (X_EDGE)
      coords(1) = self%geo%midp_coords(i, 1)
      coords(2) = self%geo%vert_coords(j, 2)
      coords(3) = self%geo%vert_coords(k, 3)
    case (Y_EDGE)
      coords(1) = self%geo%vert_coords(i, 1)
      coords(2) = self%geo%midp_coords(j, 2)
      coords(3) = self%geo%vert_coords(k, 3)
    case (Z_EDGE)
      coords(1) = self%geo%vert_coords(i, 1)
      coords(2) = self%geo%vert_coords(j, 2)
      coords(3) = self%geo%midp_coords(k, 3)
    case default
      error stop "Unknown data_loc in get_coordinates"
    end select
  end function

end module m_mesh
