submodule(m_decomp) m_decomp_2decompfft

  use mpi
  use m_decomp, only: decomp_t
  implicit none

  type, extends(decomp_t) :: decomp_2decompfft_t
  contains
    procedure :: decomposition => decomposition_2decompfft
  end type

  contains

  module subroutine decomposition_2decompfft(self, grid, par)
    !! Supports 1D, 2D, and 3D domain decomposition.
    !!
    !! Current implementation allows only constant sub-domain size across a
    !! given direction.
    use m_grid, only: grid_t
    use m_par, only: par_t
    use decomp_2d, only: decomp_2d_init, DECOMP_2D_COMM_CART_X, xsize, xstart

    class(decomp_2decompfft_t) :: self
    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par 
    integer :: p_col, p_row
    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer, allocatable, dimension(:) :: global_ranks_lin
    integer :: nproc
    integer, dimension(3) :: subd_pos, subd_pos_prev, subd_pos_next
    logical, dimension(3) :: periodic_bc
    integer :: dir
    logical :: is_last_domain
    integer :: nx, ny, nz
    integer :: ierr
    integer :: cart_rank
    integer, dimension(2) :: coords

    if (par%is_root()) then
      print*, "Domain decomposition by 2decomp&fft"
    end if

    nx = grid%global_cell_dims(1)
    ny = grid%global_cell_dims(2)
    nz = grid%global_cell_dims(3)

    p_row = par%nproc_dir(2)
    p_col = par%nproc_dir(3)
    periodic_bc(:) = grid%periodic_BC(:)
    call decomp_2d_init(nx, ny, nz, p_row, p_col, periodic_bc)

    ! Get global_ranks
    allocate(global_ranks(1, p_row, p_col))
    allocate(global_ranks_lin(p_row*p_col))
    global_ranks_lin(:) = 0

    call MPI_Comm_rank(DECOMP_2D_COMM_CART_X, cart_rank, ierr)
    call MPI_Cart_coords(DECOMP_2D_COMM_CART_X, cart_rank, 2, coords, ierr)

    global_ranks_lin(coords(1)+1 + p_row*(coords(2))) = par%nrank

    call MPI_Allreduce(MPI_IN_PLACE, global_ranks_lin, p_row*p_col, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    global_ranks = reshape(global_ranks_lin, shape=[1, p_row, p_col])

    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, par%nrank)

    ! local/directional position of the subdomain
    par%nrank_dir(:) = subd_pos(:) - 1

    ! Get local domain size and offset from 2decomp
    grid%cell_dims(:) = xsize(:)
    par%n_offset(:) = xstart(:)

    ! compute vert_dims from cell_dims
    do dir = 1, 3
      is_last_domain = (par%nrank_dir(dir) + 1 == par%nproc_dir(dir))
      if (is_last_domain .and. (.not. grid%periodic_BC(dir))) then
        grid%vert_dims(dir) = grid%cell_dims(dir) +1
      else
        grid%vert_dims(dir) = grid%cell_dims(dir)
      end if
    end do

    ! Get neighbour ranks
    do dir = 1, 3
      nproc = par%nproc_dir(dir)
      subd_pos_prev(:) = subd_pos(:)
      subd_pos_prev(dir) = modulo(subd_pos(dir) - 2, nproc) + 1
      par%pprev(dir) = global_ranks(subd_pos_prev(1), &
                                         subd_pos_prev(2), &
                                         subd_pos_prev(3))

      subd_pos_next(:) = subd_pos(:)
      subd_pos_next(dir) = modulo(subd_pos(dir) - nproc, nproc) + 1
      par%pnext(dir) = global_ranks(subd_pos_next(1), &
                                         subd_pos_next(2), &
                                         subd_pos_next(3))
    end do

  end subroutine decomposition_2decompfft


end submodule