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

end module m_common
