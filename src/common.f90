module m_common
  !! Common module containing global constants, parameters, and utility functions.
  !!
  !! This module provides:
  !! - Precision definitions (single or double precision based on compilation flags)
  !! - Mathematical constants (e.g., \(\pi\))
  !! - Direction and reordering constants for domain decomposition
  !! - Data location flags (vertex, cell, face, edge centered)
  !! - Boundary condition type constants
  !! - Utility functions for argument parsing and data manipulation
  use mpi

  implicit none

#ifdef SINGLE_PREC
  integer, parameter :: dp = kind(0.0e0)  !! Double precision kind parameter (single precision)
  integer, parameter :: nbytes = 4        !! Number of bytes for real numbers
  integer, parameter :: MPI_X3D2_DP = MPI_REAL  !! MPI datatype for real numbers
  logical, parameter :: is_sp = .true.    !! Flag indicating single precision
#else
  integer, parameter :: dp = kind(0.0d0)  !! Double precision kind parameter (double precision)
  integer, parameter :: nbytes = 8        !! Number of bytes for real numbers
  integer, parameter :: MPI_X3D2_DP = MPI_DOUBLE_PRECISION  !! MPI datatype for real numbers
  logical, parameter :: is_sp = .false.   !! Flag indicating double precision
#endif

  integer, parameter :: i8 = selected_int_kind(18)  !! Integer kind for 64-bit integers

  real(dp), parameter :: pi = 4*atan(1.0_dp)  !! Mathematical constant \(\pi\)

  !> Reordering constants for data layout transformations between directions.
  !! Format: RDR_<from_dir>2<to_dir> where directions are X, Y, Z, or C (complete/cell-centered)
  integer, parameter :: RDR_X2Y = 12, RDR_X2Z = 13, RDR_Y2X = 21, &
                        RDR_Y2Z = 23, RDR_Z2X = 31, RDR_Z2Y = 32, &
                        RDR_C2X = 41, RDR_C2Y = 42, RDR_C2Z = 43, &
                        RDR_X2C = 14, RDR_Y2C = 24, RDR_Z2C = 34
  integer, parameter :: DIR_X = 1  !! X direction index
  integer, parameter :: DIR_Y = 2  !! Y direction index
  integer, parameter :: DIR_Z = 3  !! Z direction index
  integer, parameter :: DIR_C = 4  !! Complete/cell-centered direction index
  integer, parameter :: POISSON_SOLVER_FFT = 0  !! FFT-based Poisson solver
  integer, parameter :: POISSON_SOLVER_CG = 1   !! Conjugate gradient Poisson solver
  integer, parameter :: VERT = 0000, & !! Vertex centered data
                        CELL = 1110, & !! Cell centered data
                        X_FACE = 1100, & !! Data on faces normal to X
                        Y_FACE = 1010, & !! Data on faces normal to Y
                        Z_FACE = 0110, & !! Data on faces normal to Z
                        X_EDGE = 0010, & !! Data on edges along X
                        Y_EDGE = 0100, & !! Data on edges along Y
                        Z_EDGE = 1000, & !! Data on edges along Z
                        NULL_LOC = -0001 !! The location of data isn't specified
  integer, parameter :: BC_PERIODIC = 0   !! Periodic boundary condition
  integer, parameter :: BC_NEUMANN = 1    !! Neumann boundary condition
  integer, parameter :: BC_DIRICHLET = 2  !! Dirichlet boundary condition
  integer, parameter :: BC_HALO = -1      !! Halo/ghost cell boundary condition
  !> Reordering map matrix for direction transformations.
  !! Maps from direction (row) to direction (column), yielding the reordering constant.
  integer, protected :: &
    rdr_map(4, 4) = reshape([0, RDR_Y2X, RDR_Z2X, RDR_C2X, &
                             RDR_X2Y, 0, RDR_Z2Y, RDR_C2Y, &
                             RDR_X2Z, RDR_Y2Z, 0, RDR_C2Z, &
                             RDR_X2C, RDR_Y2C, RDR_Z2C, 0], shape=[4, 4])

contains

  pure subroutine get_dirs_from_rdr(dir_from, dir_to, rdr_dir)
    !! Extract source and destination directions from a reordering constant.
    !!
    !! Given a reordering constant (e.g., RDR_X2Y), this subroutine determines
    !! the source direction and destination direction.
    integer, intent(out) :: dir_from  !! Source direction (DIR_X, DIR_Y, DIR_Z, or DIR_C)
    integer, intent(out) :: dir_to    !! Destination direction (DIR_X, DIR_Y, DIR_Z, or DIR_C)
    integer, intent(in) :: rdr_dir    !! Reordering constant (e.g., RDR_X2Y)
    integer, dimension(2) :: dirs

    dirs = findloc(rdr_map, rdr_dir)
    dir_from = dirs(1)
    dir_to = dirs(2)

  end subroutine

  pure integer function get_rdr_from_dirs(dir_from, dir_to) result(rdr_dir)
    !! Returns reordering constant based on two direction inputs.
    !!
    !! Given a source and destination direction, this function returns the
    !! corresponding reordering constant (e.g., RDR_X2Y for X to Y).
    integer, intent(in) :: dir_from  !! Source direction (DIR_X, DIR_Y, DIR_Z, or DIR_C)
    integer, intent(in) :: dir_to    !! Destination direction (DIR_X, DIR_Y, DIR_Z, or DIR_C)

    rdr_dir = rdr_map(dir_from, dir_to)
  end function get_rdr_from_dirs

  function get_argument(pos) result(arg)
    !! Retrieve a command-line argument at the specified position.
    !!
    !! This function wraps the intrinsic get_command_argument with error checking
    !! and automatic string trimming.
    integer, intent(in) :: pos  !! Position of the command-line argument (1-indexed)
    character(:), allocatable :: arg  !! The retrieved command-line argument

    character(len=200) :: temp
    integer :: stat

    call get_command_argument(pos, temp, status=stat)

    if (stat > 0) then
      error stop 'Argument retrieval failed!'
    else if (stat == -1) then
      error stop 'Argument is truncated!'
    end if

    arg = trim(temp)
  end function get_argument

  integer function move_data_loc(in_data_loc, dir, move) result(out_data_loc)
    !! Update data location by shifting along a specified direction.
    !!
    !! This function modifies a data location flag by moving it along one direction
    !! (X, Y, or Z) by a specified amount. The data location encoding uses powers of 10
    !! to represent positions in each direction.
    integer, intent(in) :: in_data_loc  !! Input data location flag
    integer, intent(in) :: dir          !! Direction to move (DIR_X, DIR_Y, or DIR_Z)
    integer, intent(in) :: move         !! Amount to move (typically -1, 0, or 1)

    out_data_loc = in_data_loc + move*(10**dir)
  end function move_data_loc

end module m_common
