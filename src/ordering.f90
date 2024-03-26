module m_ordering

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y

  implicit none
contains
   !!
   !! "Application storage" stores spatial data with a directionality for better cache locality
   !!  This set of functions converts indices from this application storage (_dir) to cartesian indices (_ijk)
   !!

   pure subroutine get_index_ijk(i, j, k, dir_i, dir_j, dir_k, dir, SZ, nx_loc, ny_loc, nz_loc)
      !! Get cartesian index from application storage directional one
    integer, intent(out) :: i, j, k                   ! cartesian indices
    integer, intent(in) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: dir                        ! direction of the applicatino storage indices
    integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc ! dimensions of the block

    select case (dir)
    case (DIR_X)
      i = dir_j
      j = mod(dir_k - 1, ny_loc/SZ)*SZ + dir_i
      k = 1 + (dir_k - 1)/(ny_loc/SZ)
    case (DIR_Y)
      i = mod(dir_k - 1, nx_loc/SZ)*SZ + dir_i
      j = dir_j
      k = 1 + (dir_k - 1)/(nx_loc/SZ)
    case (DIR_Z)
      i = mod(dir_k - 1, nx_loc/SZ)*SZ + dir_i
      j = 1 + (dir_k - 1)/(nx_loc/SZ)
      k = dir_j
    end select

  end subroutine get_index_ijk

   pure subroutine get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, SZ, nx_loc, ny_loc, nz_loc)
      !! Get application storage directional index from cartesian index
    integer, intent(out) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: i, j, k                     ! cartesian indices
    integer, intent(in) :: dir                        ! direction of the application storage indices
    integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc ! dimensions of the block

    select case (dir)
    case (DIR_X)
      dir_i = mod(j - 1, SZ) + 1
      dir_j = i
      dir_k = (ny_loc/SZ)*(k - 1) + 1 + (j - 1)/SZ
    case (DIR_Y)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = j
      dir_k = (nx_loc/SZ)*(k - 1) + 1 + (i - 1)/SZ
    case (DIR_Z)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = k
      dir_k = (nx_loc/SZ)*(j - 1) + 1 + (i - 1)/SZ
    end select

  end subroutine get_index_dir

   pure subroutine get_index_reordering(out_i, out_j, out_k, in_i, in_j, in_k, reorder_dir, SZ, nx_loc, ny_loc, nz_loc)
      !! Converts a set of application storage directional index to an other direction.
      !! The two directions are defined by the reorder_dir variable, RDR_X2Y will go from storage in X to Y etc.
    integer, intent(out) :: out_i, out_j, out_k         ! new indices in the application storage
    integer, intent(in) :: in_i, in_j, in_k             ! original indices
    integer, intent(in) :: reorder_dir
    integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc ! dimensions of the block
    integer :: i, j, k        ! Intermediary cartesian indices
    integer :: dir_in, dir_out

    select case (reorder_dir)
    case (RDR_X2Y)
      dir_in = DIR_X
      dir_out = DIR_Y
    case (RDR_X2Z)
      dir_in = DIR_X
      dir_out = DIR_Z
    case (RDR_Y2X)
      dir_in = DIR_Y
      dir_out = DIR_X
    case (RDR_Y2Z)
      dir_in = DIR_Y
      dir_out = DIR_Z
    case (RDR_Z2X)
      dir_in = DIR_Z
      dir_out = DIR_X
    case (RDR_Z2Y)
      dir_in = DIR_Z
      dir_out = DIR_Y
    end select

      call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_in, SZ, nx_loc, ny_loc, nz_loc)
      call get_index_dir(out_i, out_j, out_k, i, j, k, dir_out, SZ, nx_loc, ny_loc, nz_loc)

  end subroutine get_index_reordering

end module m_ordering
