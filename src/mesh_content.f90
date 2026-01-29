module m_mesh_content
  !! Module containing mesh content types for geometry, grid, and parallel decomposition.
  !!
  !! This module defines three main types:
  !!
  !! - `geo_t`: Geometry information including coordinates and mesh stretching
  !! - `grid_t`: Grid dimensions and boundary conditions
  !! - `par_t`: Parallel domain decomposition information

  use m_common, only: dp, pi
  implicit none

  type :: geo_t
    !! Geometry information type for domain coordinates and mesh stretching.
    !!
    !! This type stores physical domain dimensions, coordinates at grid points,
    !! and mesh stretching parameters. Coordinates and stretching factors are
    !! stored for both vertex-centered and cell-centered locations.
    real(dp) :: origin(3)  !! Coordinates of vertex (1, 1, 1)
    real(dp) :: d(3)       !! Cell size in each direction for uniform mesh
    real(dp) :: L(3)       !! Global domain dimensions in each direction
    real(dp), allocatable, dimension(:, :) :: vert_coords  !! Global coordinates at vertices
    real(dp), allocatable, dimension(:, :) :: midp_coords  !! Global coordinates at cell midpoints
    character(len=20), dimension(3) :: stretching  !! Stretching type in each direction
    logical :: stretched(3)     !! Whether each direction has stretching applied
    real(dp) :: alpha(3)        !! Stretching parameter \(\alpha\) in each direction
    real(dp) :: beta(3)         !! Stretching parameter \(\beta\) in each direction
    !> Stretching factors at vertices: \(\frac{ds}{d\xi}\), \(\frac{d^2s}{d\xi^2}\), \(\frac{d^2\xi}{ds^2}\)
    real(dp), allocatable, dimension(:, :) :: vert_ds, vert_ds2, vert_d2s
    !> Stretching factors at midpoints: \(\frac{ds}{d\xi}\), \(\frac{d^2s}{d\xi^2}\), \(\frac{d^2\xi}{ds^2}\)
    real(dp), allocatable, dimension(:, :) :: midp_ds, midp_ds2, midp_d2s
  contains
    procedure :: obtain_coordinates  !! Compute coordinates and stretching factors
  end type

  type :: grid_t
    !! Grid information type for mesh dimensions and boundary conditions.
    !!
    !! This type stores both global and local (per MPI rank) grid dimensions,
    !! accounting for both vertex-centered and cell-centered data. It also
    !! manages boundary condition information.
    integer, dimension(3) :: global_vert_dims  !! Global number of vertices in each direction
    integer, dimension(3) :: global_cell_dims  !! Global number of cells in each direction
    integer, dimension(3) :: vert_dims         !! Local number of vertices in each direction
    integer, dimension(3) :: cell_dims         !! Local number of cells in each direction
    logical, dimension(3) :: periodic_BC       !! Whether each direction has periodic BC
    integer, dimension(3, 2) :: BCs_global     !! Global boundary conditions (lower, upper) in each direction
    integer, dimension(3, 2) :: BCs            !! Local subdomain boundary conditions (lower, upper)
  contains
    procedure :: copy_cell2vert_dims  !! Copy cell_dims to vert_dims accounting for periodicity
    procedure :: copy_vert2cell_dims  !! Copy vert_dims to cell_dims accounting for periodicity
  end type

  type :: par_t
    !! Parallel domain decomposition information type.
    !!
    !! This type stores all information related to MPI domain decomposition,
    !! including rank IDs, processor grid layout, and neighbor communication
    !! information for halo exchanges.
    integer :: nrank                !! Local MPI rank ID (0-based)
    integer :: nproc                !! Total number of MPI ranks
    integer, dimension(3) :: nrank_dir   !! Local rank ID in each direction (0-based)
    integer, dimension(3) :: nproc_dir   !! Number of processors in each direction
    integer, dimension(3) :: n_offset    !! Cell offset in each direction due to decomposition
    integer, dimension(3) :: pnext       !! Rank ID of next neighbor in each direction
    integer, dimension(3) :: pprev       !! Rank ID of previous neighbor in each direction
  contains
    procedure :: is_root                        !! Check if current rank is root (rank 0)
    procedure :: compute_rank_pos_from_global   !! Compute rank position and neighbors from global map
  end type

contains

  pure function is_root(self) result(is_root_rank)
    !! Check whether the current MPI rank is the root rank.
    !!
    !! The root rank is defined as rank 0 in the MPI communicator.
    class(par_t), intent(in) :: self  !! Parallel decomposition object
    logical :: is_root_rank           !! True if this is rank 0

    is_root_rank = (self%nrank == 0)

  end function

  pure subroutine compute_rank_pos_from_global(self, global_ranks)
    !! Compute rank position and neighbor ranks from global rank map.
    !!
    !! From the 3D global rank map, this subroutine determines the position
    !! of the current rank in the processor grid and identifies the previous
    !! and next neighboring ranks in each direction for halo communication.
    !! Periodic wrapping is applied for neighbor identification.
    class(par_t), intent(inout) :: self                !! Parallel decomposition object to update
    integer, dimension(:, :, :), intent(in) :: global_ranks  !! 3D map of MPI ranks
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
    !! Copy vertex dimensions to cell dimensions accounting for periodicity.
    !!
    !! For periodic boundaries, vertex and cell dimensions are equal. For
    !! non-periodic boundaries on the last domain, cell dimensions are one
    !! less than vertex dimensions.
    class(grid_t), intent(inout) :: self  !! Grid object to update
    type(par_t), intent(in) :: par        !! Parallel decomposition info
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
    !! Copy cell dimensions to vertex dimensions accounting for periodicity.
    !!
    !! For periodic boundaries, vertex and cell dimensions are equal. For
    !! non-periodic boundaries on the last domain, vertex dimensions are one
    !! more than cell dimensions.
    class(grid_t), intent(inout) :: self  !! Grid object to update
    type(par_t), intent(in) :: par        !! Parallel decomposition info
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
    !! Compute global coordinates and stretching factors for grid points.
    !!
    !! This subroutine calculates coordinates at both vertex-centered and
    !! cell-centered locations, supporting both uniform and stretched meshes.
    !! For stretched meshes, it also computes the stretching factors
    !! \(\frac{ds}{d\xi}\), \(\frac{d^2s}{d\xi^2}\), and \(\frac{d^2\xi}{ds^2}\).
    implicit none
    class(geo_t) :: self                      !! Geometry object to populate
    integer, intent(in) :: vert_dims(3)       !! Local vertex dimensions
    integer, intent(in) :: cell_dims(3)       !! Local cell dimensions
    integer, intent(in) :: n_offset(3)        !! Cell offset due to domain decomposition

    integer :: dir, i, i_glob
    real(dp) :: L_inf, alpha, beta, r, const, s, yeta_vt, yeta_mp, coord
    real(dp), parameter :: beta_tolerance = epsilon(1._dp)

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
        self%alpha(dir) = 0._dp
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
      else
        self%stretched(dir) = .true.
        L_inf = self%L(dir)/2
        beta = self%beta(dir)
        if (beta <= beta_tolerance) then
          error stop 'Invalid beta in domain_settings'
        end if
        alpha = abs((L_inf - sqrt((pi*beta)**2 + L_inf**2))/(2*beta*L_inf))
        self%alpha(dir) = alpha

        r = sqrt((alpha*beta + 1)/(alpha*beta))
        const = sqrt(beta)/(2*sqrt(alpha)*sqrt(alpha*beta + 1))
        s = self%d(dir)/self%L(dir)

        do i = 1, vert_dims(dir)
          i_glob = i + n_offset(dir)
          select case (trim(self%stretching(dir)))
          case ('centred')
            yeta_vt = (i_glob - 1)*s
          case ('top-bottom')
            yeta_vt = (i_glob - 1)*s - 0.5_dp
          case ('bottom')
            yeta_vt = (i_glob - 1)*s/2 - 0.5_dp
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
        end do

        do i = 1, cell_dims(dir)
          i_glob = i + n_offset(dir)
          select case (trim(self%stretching(dir)))
          case ('centred')
            yeta_mp = (i_glob - 0.5_dp)*s
          case ('top-bottom')
            yeta_mp = (i_glob - 0.5_dp)*s - 0.5_dp
          case ('bottom')
            yeta_mp = (i_glob - 0.5_dp)*s/2 - 0.5_dp
          case default
            error stop 'Invalid stretching type'
          end select

          ! midpoint coordinates
          coord = const*atan2(r*sin(pi*yeta_mp), cos(pi*yeta_mp)) &
                  *(2*alpha*beta - cos(2*pi*yeta_mp) + 1) &
                  /(sin(pi*yeta_mp)**2 + alpha*beta)
          self%midp_coords(i, dir) = coord + pi*const
          self%midp_ds(i, dir) = self%L(dir)*(alpha/pi &
                                              + sin(pi*yeta_mp)**2/(pi*beta))
          self%midp_ds2(i, dir) = self%midp_ds(i, dir)**2
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
      end if
    end do

  end subroutine obtain_coordinates

end module m_mesh_content
