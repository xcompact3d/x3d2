module m_omp_mesh
    use decomp_2d_fft, only: decomp_2d_init

    type, extends(mesh_t) :: omp_mesh_t
    contains
      procedure :: domain_decomposition => domain_decomposition_2decompfft
    end type omp_mesh_t

    contains


  subroutine domain_decomposition_2decompfft(mesh)
    !! Supports 1D, 2D, and 3D domain decomposition.
    !!
    !! Current implementation allows only constant sub-domain size across a
    !! given direction.
    class(mesh_t), intent(inout) :: mesh

    nx = mesh%global_vert_dims(1)
    ny = mesh%global_vert_dims(2)
    nz = mesh%global_vert_dims(3)

    p_row = mesh%par%nproc_dir(2)
    p_col = mesh%par%nproc_dir(3)
    periodic_bc(:) = mesh%periodic_BC(:)
    call decomp_2d_init(nx, ny, nz, p_row, p_col, periodc_bc)

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