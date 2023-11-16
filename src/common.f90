module m_common
   implicit none

   integer, parameter :: dp=kind(0.0d0)
   real(dp), parameter :: pi = 4*atan(1.0_dp)

   type :: globs_t
      integer :: nx, ny, nz
      integer :: nx_loc, ny_loc, nz_loc
      integer :: n_groups_x, n_groups_y, n_groups_z
      real(dp) :: Lx, Ly, Lz
      real(dp) :: dx, dy, dz
      integer :: nproc_x = 1, nproc_y = 1, nproc_z = 1
      character(len=20) :: BC_x_s, BC_x_e, BC_y_s, BC_y_e, BC_z_s, BC_z_e
   end type globs_t

contains

   subroutine set_pprev_pnext(xprev, xnext, yprev, ynext, zprev, znext, &
                              xnproc, ynproc, znproc, nrank)
      implicit none

      integer, intent(out) :: xprev, xnext, yprev, ynext, zprev, znext
      integer, intent(in) :: xnproc, ynproc, znproc, nrank

      integer :: i, ix, iy, iz

      ! Set the pprev and pnext for each rank
      i = 0
      do iz = 1, znproc
         do iy = 1, ynproc
            do ix = 1, xnproc
               if (nrank == i) then
                  ! prev and next in x direction
                  if (ix == 1) then
                     xprev = i + (xnproc - 1)
                  else
                     xprev = i - 1
                  end if

                  if (ix == xnproc) then
                     xnext = i - (xnproc - 1)
                  else
                     xnext = i + 1
                  end if

                  ! prev and next in y direction
                  if (iy == 1) then
                     yprev = i + (xnproc*(ynproc - 1))
                  else
                     yprev = i - xnproc
                  end if

                  if (iy == ynproc) then
                     ynext = i - (xnproc*(ynproc - 1))
                  else
                     ynext = i + xnproc
                  end if

                  ! prev and next in z direction
                  if (iz == 1) then
                     zprev = i + (xnproc*ynproc*(znproc - 1))
                  else
                     zprev = i - xnproc*ynproc
                  end if

                  if (iz == znproc) then
                     znext = i - (xnproc*ynproc*(znproc - 1))
                  else
                     znext = i + xnproc*ynproc
                  end if
               end if

               ! increment rank number for the next one
               i = i + 1
            end do
         end do
      end do

   end subroutine set_pprev_pnext

end module m_common
