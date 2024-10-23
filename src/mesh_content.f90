module m_mesh_content

  use m_common, only: dp
  implicit none

  ! Stores geometry information
  type :: geo_t
    real(dp), dimension(3) :: d ! size of a cell in each direction (=edge length, distance between centers, distance between vertices)
    real(dp), dimension(3) :: L ! Global dimensions of the domain in each direction
  end type

 ! Stores grid information
  type :: grid_t
    integer, dimension(3) :: global_vert_dims ! global number of vertices in each direction without padding (cartesian structure)
    integer, dimension(3) :: global_cell_dims ! global number of cells in each direction without padding (cartesian structure)

    integer, dimension(3) :: vert_dims_padded ! local domain size including padding (cartesian structure)
    integer, dimension(3) :: vert_dims ! local number of vertices in each direction without padding (cartesian structure)
    integer, dimension(3) :: cell_dims ! local number of cells in each direction without padding (cartesian structure)
    logical, dimension(3) :: periodic_BC ! Whether or not a direction has a periodic BC
    integer, dimension(3, 2) :: BCs_global
    integer, dimension(3, 2) :: BCs
  end type

  ! Stores parallel domain related information
  type :: par_t
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

  contains

  pure function is_root(self) result(is_root_rank)
  !! Returns wether or not the current rank is the root rank
    class(par_t), intent(in) :: self
    logical :: is_root_rank

    is_root_rank = (self%nrank == 0)

  end function

end module m_mesh_content
