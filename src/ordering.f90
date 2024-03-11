module m_ordering

   implicit none

   use m_omp_common, only: SZ

contains

   ! Get cartesian index from dir index
   pure elemental subroutine get_index_ijk(i, j, k, dir_i, dir_j, dir_k, dir)
      integer, intent(out) :: i, j, k
      integer, intent(in) :: dir_i, dir_j, dir_k
      integer, intent(in) :: dir
      
      select case (dir)
         case (dir_X)
            i = dir_j
            j = mod(dir_k, ny_loc/SZ)*SZ + dir_i 
            k = floor(dir_k/(ny_loc/SZ))
         case (dir_Y)
            i = mod(dir_k, nx_loc/SZ)*SZ + dir_i 
            j = dir_j
            k = floor(dir_k/(nx_loc/SZ))
         case (dir_Z)
            i = mod(dir_k, nx_loc/SZ)*SZ + dir_i 
            j = floor(dir_k/(nx_loc/SZ))
            k = dir_j
      end select

   end subroutine get_index_ijk

   ! Get directional index from cartesian index
   pure elemental subroutine get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir)
      integer, intent(out) :: dir_i, dir_j, dir_k
      integer, intent(in) :: i, j, k
      integer, intent(in) :: dir
      
      select case (dir)
         case (dir_X)
            dir_i = mod(j, SZ) + 1
            dir_j = i
            dir_k = (ny_loc/SZ)*(k-1) + floor(j/SZ)
         case (dir_Y)
            dir_i = mod(i, SZ) + 1
            dir_j = j
            dir_k = (nx_loc/SZ)*(k-1) + floor(i/SZ)
         case (dir_Z)
            dir_i = mod(i, SZ) + 1
            dir_j = k
            dir_k = (nx_loc/SZ)*(k-1) + floor(i/SZ)
      end select

   end subroutine get_index_dir

   ! Swap between two sets of directional indices following reorder_dir
   pure elemental subroutine get_index_reordering(out_i, out_j, out_k, in_i, in_j, in_k, reorder_dir)
      integer, intent(out) :: out_i, out_j, out_k
      integer, intent(in) :: in_i, in_j, in_k
      integer, intent(in) :: reorder_dir
      integer :: i, j, k
      integer :: dir_in, dir_out

      select case (dir)
         case (RDR_X2Y)
            dir_in = dir_X
            dir_out = dir_Z
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

      call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_in)
      call get_index_dir(out_i, out_j, out_k, i, j, k, dir_out)

   end subroutine get_index_reordering


end module m_ordering
