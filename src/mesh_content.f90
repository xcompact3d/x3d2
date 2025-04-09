module m_mesh_content

  use m_common, only: dp, pi
  implicit none

  type :: geo_t
    !! Stores geometry information
    !> Origin: coordinates of vertex (1, 1, 1)
    real(dp) :: origin(3)
    !> size of a cell in each direction for a uniform mesh
    real(dp) :: d(3)
    !> Global dimensions of the domain in each direction
    real(dp) :: L(3)
    !> Global coordinates at vertices
    real(dp), allocatable, dimension(:, :) :: vert_coords
    !> Global coordinates at midpoints
    real(dp), allocatable, dimension(:, :) :: midp_coords
    !> Stretching type
    character(len=20), dimension(3) :: stretching
    !> Stretching
    logical :: stretched(3)
    !> Stretching parameter
    real(dp) :: beta(3)
    !> Stretching factors at vertices
    real(dp), allocatable, dimension(:, :) :: vert_ds, vert_ds2, vert_d2s
    !> Stretching factors at midpoints
    real(dp), allocatable, dimension(:, :) :: midp_ds, midp_ds2, midp_d2s
  contains
    procedure :: obtain_coordinates
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
    procedure :: compute_rank_pos_from_global ! fills in pnext, pprev and nrank_dir from global ranks map
  end type

contains

  pure function is_root(self) result(is_root_rank)
    !! Returns wether or not the current rank is the root rank
    class(par_t), intent(in) :: self
    logical :: is_root_rank

    is_root_rank = (self%nrank == 0)

  end function

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

  subroutine obtain_coordinates(self, vert_dims, cell_dims, n_offset)
    !! Obtains global coordinates for all the vertices and midpoints
    implicit none
    class(geo_t) :: self
    integer, intent(in) :: vert_dims(3), cell_dims(3), n_offset(3)

    integer :: dir, i, i_glob
    real(dp) :: L_inf, alpha, beta, r, const, s, yeta_vt, yeta_mp, coord

    allocate (self%vert_coords(maxval(vert_dims), 3))
    allocate (self%vert_ds(maxval(vert_dims), 3))
    allocate (self%vert_ds2(maxval(vert_dims), 3))
    allocate (self%vert_d2s(maxval(vert_dims), 3))

    allocate (self%midp_coords(maxval(cell_dims), 3))
    allocate (self%midp_ds(maxval(cell_dims), 3))
    allocate (self%midp_ds2(maxval(cell_dims), 3))
    allocate (self%midp_d2s(maxval(cell_dims), 3))

    ! vertex coordinates
    do dir = 1, 3
      if (trim(self%stretching(dir)) == 'uniform') then
        self%stretched(dir) = .false.
        self%vert_coords(1:vert_dims(dir), dir) = &
          [((n_offset(dir) + i - 1)*self%d(dir), i=1, vert_dims(dir))]
        self%vert_ds(:, dir) = 1._dp
        self%vert_ds2(:, dir) = 1._dp
        self%vert_d2s(:, dir) = 0._dp
        self%midp_coords(1:cell_dims(dir), dir) = &
          [((n_offset(dir) + i - 0.5_dp)*self%d(dir), i=1, cell_dims(dir))]
        self%midp_ds(:, dir) = 1._dp
        self%midp_ds2(:, dir) = 1._dp
        self%midp_d2s(:, dir) = 0._dp
        ! all sorted, move on to the next direction
        cycle
      else
        self%stretched(dir) = .true.
      end if

      L_inf = self%L(dir)/2
      beta = self%beta(dir)
      alpha = abs((L_inf - sqrt((pi*beta)**2 + L_inf**2))/(2*beta*L_inf))
      r = sqrt((alpha*beta + 1)/(alpha*beta))
      const = sqrt(beta)/(2*sqrt(alpha)*sqrt(alpha*beta + 1))
      s = self%d(dir)/self%L(dir)

      do i = 1, vert_dims(dir)
        i_glob = i + n_offset(dir)
        select case (trim(self%stretching(dir)))
        case ('centred')
          yeta_vt = (i_glob - 1)*s
          yeta_mp = (i_glob - 0.5_dp)*s
        case ('both-ends')
          yeta_vt = (i_glob - 1)*s - 0.5_dp
          yeta_mp = (i_glob - 0.5_dp)*s - 0.5_dp
        case ('bottom')
          yeta_vt = (i_glob - 1)*s/2 - 0.5_dp
          yeta_mp = (i_glob - 0.5_dp)*s/2 - 0.5_dp
        case default
          error stop 'Invalid stretching type'
        end select

        ! vertex coordinates
        coord = const*atan2(r*sin(pi*yeta_vt), cos(pi*yeta_vt)) &
                *(2*alpha*beta - cos(2*pi*yeta_vt) + 1) &
                /(sin(pi*yeta_vt)**2 + alpha*beta)
        self%vert_coords(i, dir) = coord + pi*const
        self%vert_ds(i, dir) = self%L(dir)*(alpha/pi &
                                            + sin(pi*yeta_vt)**2/(pi*beta))
        self%vert_ds2(i, dir) = self%vert_ds(i, dir)**2
        self%vert_d2s(i, dir) = 2*cos(pi*yeta_vt)*sin(pi*yeta_vt)/beta

        ! midpoint coordinates
        coord = const*atan2(r*sin(pi*yeta_mp), cos(pi*yeta_mp)) &
                *(2*alpha*beta - cos(2*pi*yeta_mp) + 1) &
                /(sin(pi*yeta_mp)**2 + alpha*beta)
        self%midp_coords(i, dir) = coord + pi*const
        self%midp_ds(i, dir) = self%L(dir)*(alpha/pi &
                                            + sin(pi*yeta_mp)**2/(pi*beta))
        self%midp_ds2(i, dir) = self%vert_ds(i, dir)**2
        self%midp_d2s(i, dir) = 2*cos(pi*yeta_mp)*sin(pi*yeta_mp)/beta
      end do

      ! final shifts/corrections
      select case (trim(self%stretching(dir)))
      case ('centred')
        self%vert_coords(:, dir) = self%vert_coords(:, dir) - L_inf
        self%midp_coords(:, dir) = self%midp_coords(:, dir) - L_inf
      case ('bottom')
        self%vert_coords(:, dir) = 2*(self%vert_coords(:, dir))
        self%vert_d2s(:, dir) = self%vert_d2s(:, dir)/2
        self%midp_coords(:, dir) = 2*(self%midp_coords(:, dir))
        self%midp_d2s(:, dir) = self%midp_d2s(:, dir)/2
      end select
    end do

  end subroutine obtain_coordinates

end module m_mesh_content
