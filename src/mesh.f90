module m_mesh
  use iso_fortran_env, only: stderr => error_unit

  use mpi
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                     CELL, VERT, NONE, X_FACE, Y_FACE, Z_FACE, &
                     X_EDGE, Y_EDGE, Z_EDGE
  use m_field, only: field_t

  implicit none


  ! Stores geometry information
  type :: geo_t
    real(dp), dimension(3) :: d ! size of a cell in each direction (=edge length, distance between centers, distance between vertices)
    real(dp), dimension(3) :: L ! Global dimensions of the domain in each direction
  end type

  ! Stores parallel domain related information
  type :: parallel_t
    integer :: nrank ! local rank ID
    integer :: nproc ! total number of ranks/proc participating in the domain decomposition
    integer, dimension(3) :: nrank_dir ! local rank ID in each direction
    integer, dimension(3) :: nproc_dir ! total number of proc in each direction
    integer, dimension(3) :: n_offset  ! number of cells offset in each direction due to domain decomposition
    integer, dimension(3) :: pnext ! rank ID of the previous rank in each direction
    integer, dimension(3) :: pprev ! rank ID of the next rank in each direction
  contains
    procedure :: is_root ! returns if the current rank is the root rank
  end type


  ! The mesh class stores all the information about the global and local (due to domain decomposition) mesh
  ! It also includes getter functions to access some of its parameters
  type :: mesh_t
    integer, dimension(3), private :: dims_global ! global number of vertices in each direction

    integer, dimension(3), private :: vert_dims_padded ! local domain size including padding (cartesian structure)
    integer, dimension(3), private :: vert_dims ! local number of vertices in each direction without padding (cartesian structure)
    integer, dimension(3), private :: cell_dims ! local number of cells in each direction without padding (cartesian structure)
    logical, dimension(3), private :: periodic_dir ! Whether or not a direction has a periodic BC
    integer, private :: sz
    class(geo_t), allocatable :: geo ! object containing geometry information
    class(parallel_t), allocatable :: par ! object containing parallel domain decomposition information
  contains
    procedure :: get_SZ

    procedure :: get_n_groups_dir
    procedure :: get_n_groups_phi
    generic :: get_n_groups => get_n_groups_dir, get_n_groups_phi

    procedure :: get_dims_dir
    procedure :: get_dims_phi
    procedure :: get_dims_phi_dataloc
    generic :: get_dims => get_dims_dir, get_dims_phi, get_dims_phi_dataloc

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

  function mesh_init(dims_global, nproc_dir, L_global) result(mesh)
    !! Completely initialise the mesh object. 
    !! Upon initialisation the mesh object can be read-only and shouldn't be edited
    !! Takes as argument global information about the mesh like its length, number of cells and decomposition in each direction
    integer, dimension(3), intent(in) :: dims_global
    integer, dimension(3), intent(in) :: nproc_dir ! Number of proc in each direction
    real(dp), dimension(3), intent(in) :: L_global
    type(mesh_t) :: mesh

    integer :: nx, ny, nz, dir
    integer :: ierr
    logical :: is_last_domain

    allocate(mesh%geo)
    allocate(mesh%par)
    mesh%dims_global(:) = dims_global
    ! Geometry
    mesh%geo%L = L_global
    mesh%geo%d = mesh%geo%L(:) / mesh%dims_global(:)

    ! Parallel domain decomposition
    mesh%par%nproc_dir(:) = nproc_dir
    mesh%par%nproc = product(nproc_dir(:))
    call MPI_Comm_rank(MPI_COMM_WORLD, mesh%par%nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mesh%par%nproc, ierr)
    call domain_decomposition(mesh)

    ! Local mesh dimensions
    nx = mesh%dims_global(1)/mesh%par%nproc_dir(1)
    ny = mesh%dims_global(2)/mesh%par%nproc_dir(2)
    nz = mesh%dims_global(3)/mesh%par%nproc_dir(3)


    ! TODO: fixme: hard code no periodic boundary (for now)
    mesh%periodic_dir(:) = .false.

    ! Define number of cells and vertices in each direction
    mesh%vert_dims = [nx, ny, nz]
    mesh%cell_dims = mesh%vert_dims

    do dir=1, 3
      is_last_domain = (mesh%par%nrank_dir(dir) == mesh%par%nproc_dir(dir))
      if (is_last_domain) then
        if (mesh%periodic_dir(dir)) then
          mesh%vert_dims(dir) = mesh%vert_dims(dir) + 1
          mesh%cell_dims(dir) = mesh%cell_dims(dir) + 1
        end if

        mesh%cell_dims(dir) = mesh%cell_dims(dir) - 1
      end if
    end do

  end function mesh_init

  subroutine set_padded_dims(self, vert_dims)
    class(mesh_t), intent(inout) :: self
    integer, dimension(3), intent(in) :: vert_dims

    self%vert_dims_padded = vert_dims

  end subroutine

  subroutine set_sz(self, sz)
    class(mesh_t), intent(inout) :: self
    integer, intent(in) :: sz

    self%sz = sz

  end subroutine

  subroutine domain_decomposition(mesh)
    !! Supports 1D, 2D, and 3D domain decomposition.
    !!
    !! Current implementation allows only constant sub-domain size across a
    !! given direction.
    class(mesh_t), intent(inout) :: mesh

    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer :: i, nproc_x, nproc_y, nproc_z, nproc
    integer, dimension(3) :: subd_pos, subd_pos_prev, subd_pos_next
    integer :: dir

    ! Number of processes on a direction basis
    nproc_x = mesh%par%nproc_dir(1)
    nproc_y = mesh%par%nproc_dir(2)
    nproc_z = mesh%par%nproc_dir(3)

    ! A 3D array corresponding to each region in the global domain
    allocate (global_ranks(nproc_x, nproc_y, nproc_z))

    ! set the corresponding global rank for each sub-domain
    global_ranks = reshape([(i, i=0, mesh%par%nproc - 1)], &
                           shape=[nproc_x, nproc_y, nproc_z])

    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, mesh%par%nrank)

    ! local/directional position of the subdomain
    mesh%par%nrank_dir(:) = subd_pos(:) - 1

    mesh%par%n_offset(:) = (mesh%dims_global(:)/mesh%par%nproc_dir(:)) * mesh%par%nrank_dir(:)

    do dir=1, 3
      nproc = mesh%par%nproc_dir(dir)
      subd_pos_prev(:) = subd_pos(:)
      subd_pos_prev(dir) = modulo(subd_pos(dir) - 2, nproc) + 1
      mesh%par%pprev(dir) = global_ranks(subd_pos_prev(1), &
                                         subd_pos_prev(2), &
                                         subd_pos_prev(3))

      subd_pos_next(:) = subd_pos(:)
      subd_pos_next(dir) = modulo(subd_pos(dir) - nproc, nproc) + 1
      mesh%par%pnext(dir) = global_ranks(subd_pos_next(1), &
                                         subd_pos_next(2), &
                                         subd_pos_next(3))
    end do

  end subroutine domain_decomposition

  pure function get_sz(self) result(sz)
  !! Getter for parameter SZ
    class(mesh_t), intent(in) :: self
    integer :: sz

    sz = self%sz

  end function

  pure function get_padded_dims_dir(self, dir) result(dims_padded)
  !! Getter for padded dimensions with structure in `dir` direction
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: dir
    integer, dimension(3) :: dims_padded

    if (dir == DIR_C) then
      dims_padded = self%vert_dims_padded
    else
      dims_padded(1) = self%sz 
      dims_padded(2) = self%vert_dims_padded(dir)
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

      n_groups = (product(self%vert_dims_padded(:))/self%vert_dims_padded(dir)) /self%sz

  end function

  pure function get_n_groups_phi(self, phi) result(n_groups)
  !! Getter for the number of groups for fields phi
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer :: n_groups

      n_groups = self%get_n_groups(phi%dir)

  end function

  pure function get_dims_phi(self, phi) result(dims)
  !! Getter for the dimensions of field phi
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer, dimension(3) :: dims

    dims = self%get_dims(phi%dir, phi%data_loc)

  end function

  pure function get_dims_phi_dataloc(self, phi, data_loc) result(dims)
  !! Getter for the dimensions of field phi where data is located on `data_loc`
    class(mesh_t), intent(in) :: self
    class(field_t), intent(in) :: phi
    integer, intent(in) :: data_loc
    integer, dimension(3) :: dims

    dims = self%get_dims(phi%dir, data_loc)

  end function

  pure function get_dims_dir(self, dir, data_loc) result(dims)
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

    n_cell = self%cell_dims(dir)
    n_vert = self%vert_dims(dir)

    ! default to n_vert
    n = n_vert

    select case(data_loc)
      case(CELL)
        n = n_cell
      case(VERT)
        n = n_vert
      case(X_FACE)
        if (dir /= DIR_X) then
          n = n_cell
        end if
      case(Y_FACE)
        if (dir /= DIR_Y) then
          n = n_cell
        end if
      case(Z_FACE)
        if (dir /= DIR_Z) then
          n = n_cell
        end if
      case(X_EDGE)
        if (dir == DIR_X) then
          n = n_cell
        end if
      case(Y_EDGE)
        if (dir == DIR_Y) then
          n = n_cell
        end if
      case(Z_EDGE)
        if (dir == DIR_Z) then
          n = n_cell
        end if
      case(NONE)
        error stop
    end select
  end function

  pure function get_coordinates(self, i, j, k) result(xloc)
  !! Get the physical location of a cell center with i,j,k local indices
    class(mesh_t), intent(in) :: self
    integer, intent(in) :: i, j, k
    real(dp), dimension(3) :: xloc

    xloc(1) = (i - 1 + self%par%n_offset(1))*self%geo%d(1)
    xloc(2) = (j - 1 + self%par%n_offset(2))*self%geo%d(2)
    xloc(3) = (k - 1 + self%par%n_offset(3))*self%geo%d(3)
  end function

  pure function is_root(self) result(is_root_rank)
  !! Returns whether or not the current rank is the root rank
    class(parallel_t), intent(in) :: self
    logical :: is_root_rank

    is_root_rank = (self%nrank == 0)

end function


end module m_mesh