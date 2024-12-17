module m_ordering

  use m_common, only: dp, get_dirs_from_rdr, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y

  use m_mesh, only: mesh_t

  implicit none
  interface get_index_reordering
    procedure get_index_reordering_rdr, get_index_reordering_dirs
  end interface

contains
   !!
   !! "Application storage" stores spatial data with a directionality for better cache locality
   !!  This set of functions converts indices from this application storage (_dir) to cartesian indices (_ijk)
   !!

  pure subroutine get_index_ijk(i, j, k, dir_i, dir_j, dir_k, dir, &
                                SZ, nx_padded, ny_padded, nz_padded)
      !! Get cartesian index from application storage directional one
    integer, intent(out) :: i, j, k                   ! cartesian indices
    integer, intent(in) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: dir                        ! direction of the applicatino storage indices
    integer, intent(in) :: SZ, nx_padded, ny_padded, nz_padded ! dimensions of the block

    select case (dir)
    case (DIR_X)
      i = dir_j
      j = mod(dir_k - 1, ny_padded/SZ)*SZ + dir_i
      k = 1 + (dir_k - 1)/(ny_padded/SZ)
    case (DIR_Y)
      i = mod(dir_k - 1, nx_padded/SZ)*SZ + dir_i
      j = dir_j
      k = 1 + (dir_k - 1)/(nx_padded/SZ)
    case (DIR_Z)
      i = mod(dir_k - 1, nx_padded/SZ)*SZ + dir_i
      j = 1 + (dir_k - 1)/(nx_padded/SZ)
      k = dir_j
    case (DIR_C)
      i = dir_i
      j = dir_j
      k = dir_k
    end select

  end subroutine get_index_ijk

  pure subroutine get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, &
                                SZ, nx_padded, ny_padded, nz_padded)
      !! Get application storage directional index from cartesian index
    integer, intent(out) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: i, j, k                     ! cartesian indices
    integer, intent(in) :: dir                        ! direction of the application storage indices
    integer, intent(in) :: SZ, nx_padded, ny_padded, nz_padded ! dimensions of the block

    select case (dir)
    case (DIR_X)
      dir_i = mod(j - 1, SZ) + 1
      dir_j = i
      dir_k = (ny_padded/SZ)*(k - 1) + 1 + (j - 1)/SZ
    case (DIR_Y)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = j
      dir_k = (nx_padded/SZ)*(k - 1) + 1 + (i - 1)/SZ
    case (DIR_Z)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = k
      dir_k = (nx_padded/SZ)*(j - 1) + 1 + (i - 1)/SZ
    case (DIR_C)
      dir_i = i
      dir_j = j
      dir_k = k
    end select

  end subroutine get_index_dir

  pure subroutine get_index_reordering_dirs( &
    out_i, out_j, out_k, in_i, in_j, in_k, dir_from, dir_to, mesh &
    )
      !! Converts a set of application storage directional index to an other direction.
      !! The two directions are defined by the reorder_dir variable, RDR_X2Y will go from storage in X to Y etc.
    integer, intent(out) :: out_i, out_j, out_k         ! new indices in the application storage
    integer, intent(in) :: in_i, in_j, in_k             ! original indices
    integer, intent(in) :: dir_from, dir_to
    type(mesh_t), intent(in) :: mesh
    integer :: i, j, k        ! Intermediary cartesian indices
    integer, dimension(3) :: dims_padded

    dims_padded = mesh%get_padded_dims(DIR_C)
    call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_from, mesh%get_sz(), &
                       dims_padded(1), dims_padded(2), dims_padded(3))
    call get_index_dir(out_i, out_j, out_k, i, j, k, dir_to, mesh%get_sz(), &
                       dims_padded(1), dims_padded(2), dims_padded(3))

  end subroutine get_index_reordering_dirs

  pure subroutine get_index_reordering_rdr(out_i, out_j, out_k, &
                                           in_i, in_j, in_k, reorder_dir, mesh)
    integer, intent(out) :: out_i, out_j, out_k         ! new indices in the application storage
    integer, intent(in) :: in_i, in_j, in_k             ! original indices
    integer, intent(in) :: reorder_dir
    type(mesh_t), intent(in) :: mesh
    integer :: dir_from, dir_to

    call get_dirs_from_rdr(dir_from, dir_to, reorder_dir)
    call get_index_reordering(out_i, out_j, out_k, in_i, in_j, in_k, &
                              dir_from, dir_to, mesh)

  end subroutine get_index_reordering_rdr

end module m_ordering
