module m_mesh
  use iso_fortran_env, only: stderr => error_unit

  use mpi
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      CELL, VERT, none, X_FACE, Y_FACE, Z_FACE, &
                      X_EDGE, Y_EDGE, Z_EDGE
  use m_field, only: field_t
  use m_par
  use m_grid

  implicit none

  ! Stores geometry information
  type :: geo_t
    real(dp), dimension(3) :: d ! size of a cell in each direction (=edge length, distance between centers, distance between vertices)
    real(dp), dimension(3) :: L ! Global dimensions of the domain in each direction
  end type

  ! The mesh class stores all the information about the global and local (due to domain decomposition) mesh
  ! It also includes getter functions to access some of its parameters
  type :: mesh_t
    integer, private :: sz
    type(geo_t), allocatable :: geo ! object containing geometry information
    class(grid_t), allocatable :: grid ! object containing grid information
    class(par_t), allocatable :: par ! object containing parallel domain decomposition information
  contains
    procedure :: get_SZ

    procedure :: get_dims
    procedure :: get_global_dims

    procedure :: get_n_groups_dir
    procedure :: get_n_groups_phi
    generic :: get_n_groups => get_n_groups_dir, get_n_groups_phi

    procedure :: get_field_dims_dir
    procedure :: get_field_dims_phi
    procedure :: get_field_dims_phi_dataloc
    generic :: get_field_dims => get_field_dims_dir, get_field_dims_phi, &
      get_field_dims_phi_dataloc

    procedure :: get_n_dir
    procedure :: get_n_phi
    generic :: get_n => get_n_dir, get_n_phi

    procedure :: get_padded_dims_phi
    procedure :: get_padded_dims_dir
    generic :: get_padded_dims => get_padded_dims_dir, get_padded_dims_phi

    procedure :: get_coordinates

    procedure :: set_sz
    procedure :: set_padded_dims
  end type mesh_t

  interface mesh_t
    module procedure mesh_init
  end interface mesh_t

contains

  function mesh_init(dims_global, nproc_dir, L_global, &
                     periodic_BC) result(mesh)

    use m_decomp, only: decomp_t
    !! Completely initialise the mesh object.
    !! Upon initialisation the mesh object can be read-only and shouldn't be edited
    !! Takes as argument global information about the mesh like its length, number of cells and decomposition in each direction
    integer, dimension(3), intent(in) :: dims_global
    integer, dimension(3), intent(in) :: nproc_dir ! Number of proc in each direction
    real(dp), dimension(3), intent(in) :: L_global
    logical, dimension(3), optional, intent(in) :: periodic_BC
    class(mesh_t), allocatable :: mesh
    class(decomp_t), allocatable :: decomp

    integer :: dir
    integer :: ierr

    allocate(mesh)
    allocate (mesh%geo)
    allocate (mesh%par)
    mesh%grid%global_vert_dims(:) = dims_global

    if (present(periodic_BC)) then
      mesh%grid%periodic_BC(:) = periodic_BC
    else
      ! Default to periodic BC
      mesh%grid%periodic_BC(:) = .true.
    end if

    do dir = 1, 3
      if (mesh%grid%periodic_BC(dir)) then
        mesh%grid%global_cell_dims(dir) = mesh%grid%global_vert_dims(dir)
      else
        mesh%grid%global_cell_dims(dir) = mesh%grid%global_vert_dims(dir) - 1
      end if
    end do

    ! Geometry
    mesh%geo%L = L_global
    mesh%geo%d = mesh%geo%L/mesh%grid%global_cell_dims

    ! Parallel domain decomposition
    mesh%par%nproc_dir(:) = nproc_dir
    mesh%par%nproc = product(nproc_dir(:))
    call MPI_Comm_rank(MPI_COMM_WORLD, mesh%par%nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mesh%par%nproc, ierr)

    call decomp%decomposition(mesh%grid, mesh%par)

    ! Define default values
    mesh%grid%vert_dims_padded = mesh%grid%vert_dims
    mesh%sz = 1

  end function mesh_init

  subroutine set_padded_dims(self, vert_dims)
    class(mesh_t), intent(inout) :: self
    integer, dimension(3), intent(in) :: vert_dims

    self%grid%vert_dims_padded = vert_dims

  end subroutine

  subroutine set_sz(self, sz)
    class(mesh_t), intent(inout) :: self
    integer, intent(in) :: sz

    self%sz = sz

  end subroutine

  pure function get_sz(self) result(sz)
  !! Getter for parameter SZ
    class(mesh_t), intent(in) :: self
    integer :: sz

    sz = self%sz

  end function

  pure function get_dims(self, data_loc) result(dims)
  !! Getter for local domain dimensions
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: data_loc
    integer, dimension(3) :: dims

    dims = get_dims_dataloc(data_loc, self%grid%vert_dims, self%grid%cell_dims)
  end function

  pure function get_global_dims(self, data_loc) result(dims)
  !! Getter for local domain dimensions
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: data_loc
    integer, dimension(3) :: dims

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
    case (none)
      error stop "Unknown location in get_dims_dataloc"
    end select
  end function get_dims_dataloc

  pure function get_padded_dims_dir(self, dir) result(dims_padded)
  !! Getter for padded dimensions with structure in `dir` direction
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: dir
    integer, dimension(3) :: dims_padded

    if (dir == DIR_C) then
      dims_padded = self%grid%vert_dims_padded
    else
      dims_padded(1) = self%sz
      dims_padded(2) = self%grid%vert_dims_padded(dir)
      dims_padded(3) = self%get_n_groups(dir)
    end if

  end function

  pure function get_padded_dims_phi(self, phi) result(dims_padded)
  !! Getter for padded dimensions for field phi
  !! Gets the field direction from the field itself
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer, dimension(3) :: dims_padded

    dims_padded = self%get_padded_dims(phi%dir)

  end function

  pure function get_n_groups_dir(self, dir) result(n_groups)
  !! Getter for the number of groups for fields in direction `dir`
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: dir
    integer :: n_groups

    n_groups = (product(self%grid%vert_dims_padded(:))/ &
                self%grid%vert_dims_padded(dir))/self%sz

  end function

  pure function get_n_groups_phi(self, phi) result(n_groups)
  !! Getter for the number of groups for fields phi
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer :: n_groups

    n_groups = self%get_n_groups(phi%dir)

  end function

  pure function get_field_dims_phi(self, phi) result(dims)
  !! Getter for the dimensions of field phi
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer, dimension(3) :: dims

    dims = self%get_field_dims(phi%dir, phi%data_loc)

  end function

  pure function get_field_dims_phi_dataloc(self, phi, data_loc) result(dims)
  !! Getter for the dimensions of field phi where data is located on `data_loc`
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer, intent(in) :: data_loc
    integer, dimension(3) :: dims

    dims = self%get_field_dims(phi%dir, data_loc)

  end function

  pure function get_field_dims_dir(self, dir, data_loc) result(dims)
  !! Getter for the dimensions of an array directed along `dir` where data would be located on `data_loc`
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: dir
    integer, intent(in) :: data_loc
    integer, dimension(3) :: dims

    if (dir == DIR_C) then
      dims(1) = self%get_n(DIR_X, data_loc)
      dims(2) = self%get_n(DIR_Y, data_loc)
      dims(3) = self%get_n(DIR_Z, data_loc)
    else
      dims(1) = self%sz
      dims(2) = self%get_n(dir, data_loc)
      dims(3) = self%get_n_groups(dir)
    end if

  end function

  pure function get_n_phi(self, phi) result(n)
  !! Getter for the main dimension of field phi
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer :: n

    n = self%get_n(phi%dir, phi%data_loc)

  end function

  pure function get_n_dir(self, dir, data_loc) result(n)
  !! Getter for the main dimension a field oriented along `dir` with data on `data_loc`
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: dir
    integer, intent(in) :: data_loc
    integer :: n, n_cell, n_vert

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
    case (none)
      error stop "Unknown direction in get_n_dir"
    end select
  end function get_n_dir

  pure function get_coordinates(self, i, j, k) result(xloc)
  !! Get the physical location of a cell center with i,j,k local indices
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: i, j, k
    real(dp), dimension(3) :: xloc

    xloc(1) = (i - 1 + self%par%n_offset(1))*self%geo%d(1)
    xloc(2) = (j - 1 + self%par%n_offset(2))*self%geo%d(2)
    xloc(3) = (k - 1 + self%par%n_offset(3))*self%geo%d(3)
  end function

end module m_mesh
