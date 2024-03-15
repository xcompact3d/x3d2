module m_ordering

   use m_common, only: dp, dir_X, dir_Y, dir_Z, &
                       RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y

   implicit none
contains
   !!
   !! "Application storage" stores spatial data with a directionality for better cache locality
   !!  This set of functions converts indices from this application storage (_dir) to cartesian indices (_ijk)
   !! 

   ! Get cartesian index from application storage directional one
   pure subroutine get_index_ijk(i, j, k, dir_i, dir_j, dir_k, dir, SZ, nx_loc, ny_loc, nz_loc)
      integer, intent(out) :: i, j, k
      integer, intent(in) :: dir_i, dir_j, dir_k
      integer, intent(in) :: dir
      integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc
      
      select case (dir)
         case (dir_X)
            i = dir_j
            j = mod(dir_k - 1, ny_loc/SZ)*SZ + dir_i 
            k = ceiling(real(dir_k)/(ny_loc/SZ))
         case (dir_Y)
            i = mod(dir_k - 1, nx_loc/SZ)*SZ + dir_i 
            j = dir_j
            k = ceiling(real(dir_k)/(nx_loc/SZ))
         case (dir_Z)
            i = mod(dir_k - 1, nx_loc/SZ)*SZ + dir_i 
            j = ceiling(real(dir_k)/(nx_loc/SZ))
            k = dir_j
      end select

   end subroutine get_index_ijk

   ! Get application storage directional index from cartesian index
   pure subroutine get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, SZ, nx_loc, ny_loc, nz_loc)
      integer, intent(out) :: dir_i, dir_j, dir_k
      integer, intent(in) :: i, j, k
      integer, intent(in) :: dir
      integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc
      
      select case (dir)
         case (dir_X)
            dir_i = mod(j-1, SZ) + 1
            dir_j = i
            dir_k = int(ny_loc/SZ)*(k-1) + ceiling(real(j)/SZ)
         case (dir_Y)
            dir_i = mod(i-1, SZ) + 1
            dir_j = j
            dir_k = int(nx_loc/SZ)*(k-1) + ceiling(real(i)/SZ)
         case (dir_Z)
            dir_i = mod(i-1, SZ) + 1
            dir_j = k
            dir_k = int(nx_loc/SZ)*(j-1) + ceiling(real(i)/SZ)
      end select

   end subroutine get_index_dir

   ! Converts a set of application storage directional index to an other direction. 
   ! The two directions are defined by the reorder_dir variable, RDR_X2Y will go from storage in X to Y etc.
   pure subroutine get_index_reordering(out_i, out_j, out_k, in_i, in_j, in_k, reorder_dir, SZ, nx_loc, ny_loc, nz_loc)
      integer, intent(out) :: out_i, out_j, out_k
      integer, intent(in) :: in_i, in_j, in_k
      integer, intent(in) :: reorder_dir
      integer, intent(in) :: SZ, nx_loc, ny_loc, nz_loc
      integer :: i, j, k
      integer :: dir_in, dir_out

      select case (reorder_dir)
         case (RDR_X2Y)
            dir_in = dir_X
            dir_out = dir_Y
         case (RDR_X2Z)
            dir_in = dir_X
            dir_out = dir_Z
         case (RDR_Y2X)
            dir_in = dir_Y
            dir_out = dir_X
         case (RDR_Y2Z)
            dir_in = dir_Y
            dir_out = dir_Z
         case (RDR_Z2X)
            dir_in = dir_Z
            dir_out = dir_X
         case (RDR_Z2Y)
            dir_in = dir_Z
            dir_out = dir_Y
      end select

      call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_in, SZ, nx_loc, ny_loc, nz_loc)
      call get_index_dir(out_i, out_j, out_k, i, j, k, dir_out, SZ, nx_loc, ny_loc, nz_loc)

   end subroutine get_index_reordering


end module m_ordering
