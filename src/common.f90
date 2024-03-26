  module m_common
    implicit none

    integer, parameter :: dp=kind(0.0d0)
    real(dp), parameter :: pi = 4*atan(1.0_dp)

    integer, parameter :: RDR_X2Y = 12, RDR_X2Z = 13, RDR_Y2X = 21, &
      RDR_Y2Z = 23, RDR_Z2X = 31, RDR_Z2Y = 32, &
      RDR_C2X = 41, RDR_C2Y = 42, RDR_C2Z = 43, &
      RDR_X2C = 14, RDR_Y2C = 24, RDR_Z2C = 34
    integer, parameter :: DIR_X = 1, DIR_Y = 2, DIR_Z = 3, DIR_C = 4
    integer, parameter :: POISSON_SOLVER_FFT = 0, POISSON_SOLVER_CG = 1

    integer, protected :: &
      rdr_map(4, 4) = reshape([0, RDR_X2Y, RDR_X2Z, RDR_X2C, &
      RDR_Y2X, 0, RDR_Y2Z, RDR_Y2C, &
      RDR_Z2X, RDR_Z2Y, 0, RDR_Z2C, &
      RDR_C2X, RDR_C2Y, RDR_C2Z, 0], shape=[4, 4])

    type :: globs_t
      integer :: nx, ny, nz
      integer :: nx_loc, ny_loc, nz_loc
      integer :: n_groups_x, n_groups_y, n_groups_z
      real(dp) :: Lx, Ly, Lz
      real(dp) :: dx, dy, dz
      real(dp) :: nu, dt
      integer :: n_iters, n_output
      integer :: nproc_x = 1, nproc_y = 1, nproc_z = 1
      character(len=20) :: BC_x_s, BC_x_e, BC_y_s, BC_y_e, BC_z_s, BC_z_e
      integer :: poisson_solver_type
    end type globs_t

  contains

    subroutine set_pprev_pnext(xprev, xnext, yprev, ynext, zprev, znext, &
      xnproc, ynproc, znproc, nrank)
      implicit none

      integer, intent(out) :: xprev, xnext, yprev, ynext, zprev, znext
      integer, intent(in) :: xnproc, ynproc, znproc, nrank

      integer :: ix, iy, iz

      ix = modulo(nrank, xnproc)
      iy = modulo((nrank - ix)/xnproc, ynproc)
      iz = (nrank - ix - iy*xnproc)/(xnproc*ynproc)
      ! nrank == ix + iy*xnproc + iz*xnproc*ynproc

      ! prev and next in x direction
      if (ix == 0) then
        xprev = nrank + (xnproc - 1)
      else
        xprev = nrank - 1
      end if

      if (ix == xnproc - 1) then
        xnext = nrank - (xnproc - 1)
      else
        xnext = nrank + 1
      end if

      ! prev and next in y direction
      if (iy == 0) then
        yprev = nrank + (xnproc*(ynproc - 1))
      else
        yprev = nrank - xnproc
      end if

      if (iy == ynproc - 1) then
        ynext = nrank - (xnproc*(ynproc - 1))
      else
        ynext = nrank + xnproc
      end if

      ! prev and next in z direction
      if (iz == 0) then
        zprev = nrank + (xnproc*ynproc*(znproc - 1))
      else
        zprev = nrank - xnproc*ynproc
      end if

      if (iz == znproc - 1) then
        znext = nrank - (xnproc*ynproc*(znproc - 1))
      else
        znext = nrank + xnproc*ynproc
      end if

    end subroutine set_pprev_pnext

    integer function get_rdr_from_dirs(dir_from, dir_to) result(rdr_dir)
      !! Returns RDR_?2? value based on two direction inputs
      integer, intent(in) :: dir_from, dir_to

      rdr_dir = rdr_map(dir_from, dir_to)
    end function get_rdr_from_dirs

  end module m_common
