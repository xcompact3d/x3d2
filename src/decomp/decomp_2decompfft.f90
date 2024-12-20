submodule(m_decomp) m_decomp_2decompfft
!! Parallel decomposition provided by 2decomp&FFT

  use mpi
  implicit none

  type, extends(decomp_t) :: decomp_2decompfft_t
  contains
    procedure :: decomposition => decomposition_2decompfft
  end type

  contains

  module subroutine init_decomp(decomp)
    class(decomp_t), allocatable, intent(out):: decomp

    allocate(decomp_2decompfft_t :: decomp)
    decomp%is_avail_2decomp = .true.
  end subroutine

  module subroutine decomposition_2decompfft(self, grid, par)
    !! Performs 2D mesh decomposition using 2decomp&fft
    use m_mesh_content, only: par_t, grid_t
    use decomp_2d, only: decomp_2d_init, DECOMP_2D_COMM_CART_X, xsize, xstart
    use decomp_2d_mpi, only: nrank, nproc

    class(decomp_2decompfft_t) :: self
    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par 
    integer :: p_col, p_row
    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer, allocatable, dimension(:) :: global_ranks_lin
    integer :: nproc
    logical, dimension(3) :: periodic_bc
    integer :: nx, ny, nz
    integer :: ierr
    integer :: cart_rank
    integer, dimension(2) :: coords

    if (par%is_root()) then
      print*, "Domain decomposition by 2decomp&fft"
    end if
    nrank = par%nrank
    nproc = par%nproc

    nx = grid%global_cell_dims(1)
    ny = grid%global_cell_dims(2)
    nz = grid%global_cell_dims(3)

    p_row = par%nproc_dir(2)
    p_col = par%nproc_dir(3)
    if (p_row*p_col /= par%nproc) then
      error stop "Decomposition in X not supported by 2decomp&fft backend"
    end if
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

    ! Get local domain size and offset from 2decomp
    grid%cell_dims(:) = xsize(:)
    par%n_offset(:) = xstart(:)

    call par%compute_rank_pos_from_global(global_ranks)
    call grid%copy_cell2vert_dims(par)

  end subroutine decomposition_2decompfft

end submodule