module m_common
  implicit none

  integer, parameter :: dp = kind(0.0d0)
  real(dp), parameter :: pi = 4*atan(1.0_dp)

  integer, parameter :: RDR_X2Y = 12, RDR_X2Z = 13, RDR_Y2X = 21, &
                        RDR_Y2Z = 23, RDR_Z2X = 31, RDR_Z2Y = 32, &
                        RDR_C2X = 41, RDR_C2Y = 42, RDR_C2Z = 43, &
                        RDR_X2C = 14, RDR_Y2C = 24, RDR_Z2C = 34
  integer, parameter :: DIR_X = 1, DIR_Y = 2, DIR_Z = 3, DIR_C = 4
  integer, parameter :: POISSON_SOLVER_FFT = 0, POISSON_SOLVER_CG = 1
  integer, parameter :: VERT = 1, & ! Vertex centered data
                        CELL = 2, & ! Cell centered data
                        X_FACE = 11, &  ! Data on faces normal to X
                        Y_FACE = 12, &  ! Data on faces normal to Y
                        Z_FACE = 13, &  ! Data on faces normal to Z
                        X_EDGE = 101, & ! Data on edges along X
                        Y_EDGE = 102, & ! Data on edges along Y
                        Z_EDGE = 103, & ! Data on edges along Z
                        none = -1 ! The location of data isn't specified
  integer, protected :: &
    rdr_map(4, 4) = reshape([0, RDR_Y2X, RDR_Z2X, RDR_C2X, &
                             RDR_X2Y, 0, RDR_Z2Y, RDR_C2Y, &
                             RDR_X2Z, RDR_Y2Z, 0, RDR_C2Z, &
                             RDR_X2C, RDR_Y2C, RDR_Z2C, 0], shape=[4, 4])

contains

  pure subroutine get_dirs_from_rdr(dir_from, dir_to, rdr_dir)
    integer, intent(out) :: dir_from, dir_to
    integer, intent(in) :: rdr_dir
    integer, dimension(2) :: dirs

    dirs = findloc(rdr_map, rdr_dir)
    dir_from = dirs(1)
    dir_to = dirs(2)

  end subroutine

  pure integer function get_rdr_from_dirs(dir_from, dir_to) result(rdr_dir)
      !! Returns RDR_?2? value based on two direction inputs
    integer, intent(in) :: dir_from, dir_to

    rdr_dir = rdr_map(dir_from, dir_to)
  end function get_rdr_from_dirs

end module m_common
