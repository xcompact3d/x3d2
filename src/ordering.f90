module m_ordering

  use m_common, only: dp, get_dirs_from_rdr, DIR_X, DIR_Y, DIR_Z, DIR_C

  implicit none

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

  pure subroutine get_index_reordering( &
    out_i, out_j, out_k, in_i, in_j, in_k, dir_from, dir_to, sz, cart_padded &
    )
    !! Converts indices in between any two DIR_?
    integer, intent(out) :: out_i, out_j, out_k ! output indices
    integer, intent(in) :: in_i, in_j, in_k ! input indices
    integer, intent(in) :: dir_from, dir_to
    integer, intent(in) :: sz
    integer, intent(in) :: cart_padded(3) ! padded cartesian dimensions
    integer :: i, j, k        ! Intermediary cartesian indices

    call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_from, sz, &
                       cart_padded(1), cart_padded(2), cart_padded(3))
    call get_index_dir(out_i, out_j, out_k, i, j, k, dir_to, sz, &
                       cart_padded(1), cart_padded(2), cart_padded(3))

  end subroutine get_index_reordering

end module m_ordering
