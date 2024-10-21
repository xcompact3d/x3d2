submodule(m_decomp) m_decomp_generic

  use m_decomp, only: decomp_t
  implicit none

  type, extends(decomp_t) :: decomp_generic_t
  contains
    procedure :: decomposition => decomposition_generic
  end type

  contains

  module subroutine decomposition_generic(self, grid, par)
    use m_mesh_content, only: par_t, grid_t

    class(decomp_generic_t) :: self
    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par 
    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer :: i, nproc_x, nproc_y, nproc_z, nproc
    integer, dimension(3) :: subd_pos, subd_pos_prev, subd_pos_next
    integer :: dir
    logical :: is_last_domain

    if (par%is_root()) then
      print*, "Domain decomposition by x3d2 (generic)"
    end if

    ! Number of processes on a direction basis
    nproc_x = par%nproc_dir(1)
    nproc_y = par%nproc_dir(2)
    nproc_z = par%nproc_dir(3)

    ! Define number of cells and vertices in each direction
    grid%vert_dims = grid%global_vert_dims/par%nproc_dir

    ! A 3D array corresponding to each region in the global domain
    allocate (global_ranks(nproc_x, nproc_y, nproc_z))

    ! set the corresponding global rank for each sub-domain
    global_ranks = reshape([(i, i=0, par%nproc - 1)], &
                           shape=[nproc_x, nproc_y, nproc_z])

    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, par%nrank)

    ! local/directional position of the subdomain
    par%nrank_dir(:) = subd_pos(:) - 1

    do dir = 1, 3
      is_last_domain = (par%nrank_dir(dir) + 1 == par%nproc_dir(dir))
      if (is_last_domain .and. (.not. grid%periodic_BC(dir))) then
        grid%cell_dims(dir) = grid%vert_dims(dir) - 1
      else
        grid%cell_dims(dir) = grid%vert_dims(dir)
      end if
    end do

    par%n_offset(:) = grid%vert_dims(:)*par%nrank_dir(:)

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

  end subroutine decomposition_generic


end submodule