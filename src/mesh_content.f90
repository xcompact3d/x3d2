module m_mesh_content

  use m_common, only: dp
  implicit none

  type :: geo_t
    !! Stores geometry information
    real(dp), dimension(3) :: d ! size of a cell in each direction (=edge length, distance between centers, distance between vertices)
    real(dp), dimension(3) :: L ! Global dimensions of the domain in each direction
  end type

  type :: grid_t
    !! Stores grid information
    integer, dimension(3) :: global_vert_dims ! global number of vertices in each direction without padding (cartesian structure)
    integer, dimension(3) :: global_cell_dims ! global number of cells in each direction without padding (cartesian structure)

    integer, dimension(3) :: vert_dims_padded ! local domain size including padding (cartesian structure)
    integer, dimension(3) :: vert_dims ! local number of vertices in each direction without padding (cartesian structure)
    integer, dimension(3) :: cell_dims ! local number of cells in each direction without padding (cartesian structure)
    logical, dimension(3) :: periodic_BC ! Whether or not a direction has a periodic BC
    integer, dimension(3, 2) :: BCs_global
    integer, dimension(3, 2) :: BCs
  contains
    procedure :: copy_cell2vert_dims  ! Copies cell_dims to vert_dims taking periodicity into account
    procedure :: copy_vert2cell_dims  ! Copies vert_dims to cell_dims taking periodicity into account
  end type

  type :: par_t
    !! Stores parallel domain related information
    integer :: nrank ! local rank ID
    integer :: nproc ! total number of ranks/proc participating in the domain decomposition
    integer, dimension(3) :: nrank_dir ! local rank ID in each direction
    integer, dimension(3) :: nproc_dir ! total number of proc in each direction
    integer, dimension(3) :: n_offset  ! number of cells offset in each direction due to domain decomposition
    integer, dimension(3) :: pnext ! rank ID of the previous rank in each direction
    integer, dimension(3) :: pprev ! rank ID of the next rank in each direction
  contains
    procedure :: is_root ! returns if the current rank is the root rank
    procedure :: compute_global_rank_layout   ! determines the layout of ranks in a 3-D grid
    procedure :: compute_rank_pos_from_global ! fills in pnext, pprev and nrank_dir from global ranks map
  end type

contains

  pure function is_root(self) result(is_root_rank)
    !! Returns wether or not the current rank is the root rank
    class(par_t), intent(in) :: self
    logical :: is_root_rank

    is_root_rank = (self%nrank == 0)

  end function

  pure function compute_global_rank_layout(self) result(global_rank_layout)
    class(par_t), intent(in) :: self
    integer, dimension(self%nproc_dir(1), self%nproc_dir(2), self%nproc_dir(3)) :: global_rank_layout

    integer :: i

    global_rank_layout = reshape([(i, i=0, self%nproc - 1)], &
                                 shape=[self%nproc_dir(1), &
                                        self%nproc_dir(2), &
                                        self%nproc_dir(3)])

  end function compute_global_rank_layout
  
  pure subroutine compute_rank_pos_from_global(self, global_ranks)
    !! From the global rank maps, fills in the rank position as well
    !! as the previous and next rank in the `par` structure

    class(par_t), intent(inout) :: self
    integer, dimension(:, :, :), intent(in) :: global_ranks
    integer, dimension(3) :: subd_pos, subd_pos_prev, subd_pos_next
    integer :: dir, nproc

    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, self%nrank)

    ! local/directional position of the subdomain
    self%nrank_dir(:) = subd_pos(:) - 1

    do dir = 1, 3
      nproc = self%nproc_dir(dir)
      subd_pos_prev(:) = subd_pos(:)
      subd_pos_prev(dir) = modulo(subd_pos(dir) - 2, nproc) + 1
      self%pprev(dir) = global_ranks(subd_pos_prev(1), &
                                     subd_pos_prev(2), &
                                     subd_pos_prev(3))

      subd_pos_next(:) = subd_pos(:)
      subd_pos_next(dir) = modulo(subd_pos(dir) - nproc, nproc) + 1
      self%pnext(dir) = global_ranks(subd_pos_next(1), &
                                     subd_pos_next(2), &
                                     subd_pos_next(3))
    end do

  end subroutine

  pure subroutine copy_vert2cell_dims(self, par)
    !! Copies vert_dims information to cell_dims taking
    !! periodicity into account
    class(grid_t), intent(inout) :: self
    type(par_t), intent(in) :: par
    integer :: dir
    logical :: is_last_domain

    do dir = 1, 3
      is_last_domain = (par%nrank_dir(dir) + 1 == par%nproc_dir(dir))
      if (is_last_domain .and. (.not. self%periodic_BC(dir))) then
        self%cell_dims(dir) = self%vert_dims(dir) - 1
      else
        self%cell_dims(dir) = self%vert_dims(dir)
      end if
    end do

  end subroutine

  pure subroutine copy_cell2vert_dims(self, par)
    !! Copies cell_dims information to vert_dims taking
    !! periodicity into account
    class(grid_t), intent(inout) :: self
    type(par_t), intent(in) :: par
    integer :: dir
    logical :: is_last_domain

    do dir = 1, 3
      is_last_domain = (par%nrank_dir(dir) + 1 == par%nproc_dir(dir))
      if (is_last_domain .and. (.not. self%periodic_BC(dir))) then
        self%vert_dims(dir) = self%cell_dims(dir) + 1
      else
        self%vert_dims(dir) = self%cell_dims(dir)
      end if
    end do

  end subroutine

end module m_mesh_content
