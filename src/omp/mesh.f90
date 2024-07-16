module m_omp_mesh
    use mpi
    use decomp_2d, only: decomp_2d_init, DECOMP_2D_COMM_CART_X, xsize, xstart
    use m_mesh, only: mesh_t
    use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                        CELL, VERT, none, X_FACE, Y_FACE, Z_FACE, &
                        X_EDGE, Y_EDGE, Z_EDGE
    use m_field, only: field_t


    type, extends(mesh_t) :: omp_mesh_t
    contains
      procedure :: domain_decomposition => domain_decomposition_2decompfft
    end type omp_mesh_t

  interface omp_mesh_t
    module procedure init
  end interface omp_mesh_t

  private init

    contains


  function init(dims_global, nproc_dir, L_global, &
                     periodic_BC) result(mesh)
    integer, dimension(3), intent(in) :: dims_global
    integer, dimension(3), intent(in) :: nproc_dir ! Number of proc in each direction
    real(dp), dimension(3), intent(in) :: L_global
    logical, dimension(3), optional, intent(in) :: periodic_BC

    type(omp_mesh_t) :: mesh

    mesh%mesh_t = mesh_t(dims_global, nproc_dir, L_global, &
                     periodic_BC)

  end function


  subroutine domain_decomposition_2decompfft(mesh)
    !! Supports 1D, 2D, and 3D domain decomposition.
    !!
    !! Current implementation allows only constant sub-domain size across a
    !! given direction.
    class(omp_mesh_t), intent(inout) :: mesh
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

    nx = mesh%global_vert_dims(1)
    ny = mesh%global_vert_dims(2)
    nz = mesh%global_vert_dims(3)

    p_row = mesh%par%nproc_dir(2)
    p_col = mesh%par%nproc_dir(3)
    periodic_bc(:) = mesh%periodic_BC(:)
    call decomp_2d_init(nx, ny, nz, p_row, p_col, periodic_bc)

    mesh%vert_dims(:) = xsize(:)
    mesh%par%n_offset(:) = xstart(:)

    ! Get global_ranks
    allocate(global_ranks(1, p_row, p_col))
    allocate(global_ranks_lin(p_row*p_col))
    global_ranks_lin(:) = 0

    call MPI_Comm_rank(DECOMP_2D_COMM_CART_X, cart_rank, ierr)
    call MPI_Cart_coords(DECOMP_2D_COMM_CART_X, cart_rank, 2, coords, ierr)

    global_ranks_lin(coords(1) + p_row*(coords(2)-1)) = mesh%par%nrank

    call MPI_Allreduce(MPI_IN_PLACE, global_ranks_lin, p_row*p_col, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    global_ranks = reshape(global_ranks_lin, shape=[1, p_row, p_col])

    ! TODO: refactirise
    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, mesh%par%nrank)

    ! local/directional position of the subdomain
    mesh%par%nrank_dir(:) = subd_pos(:) - 1

    do dir = 1, 3
      is_last_domain = (mesh%par%nrank_dir(dir) + 1 == mesh%par%nproc_dir(dir))
      if (is_last_domain .and. (.not. mesh%periodic_BC(dir))) then
        mesh%cell_dims(dir) = mesh%vert_dims(dir) - 1
      else
        mesh%cell_dims(dir) = mesh%vert_dims(dir)
      end if
    end do

    mesh%par%n_offset(:) = mesh%vert_dims(:)*mesh%par%nrank_dir(:)

    do dir = 1, 3
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


  end subroutine domain_decomposition_2decompfft


end module m_omp_mesh