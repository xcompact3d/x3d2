module m_dist_reduced
  !! Exact solve of the reduced system of the distributed tridiagonal solver.
  !!
  !! The first phase of the distributed solve (der_univ_dist) reduces each
  !! rank's block of a line so that its interior depends only on the block's
  !! own first and last unknowns, x_s and x_e. What couples the ranks is then
  !! a reduced system in the x_s and x_e of every rank along the line, in the
  !! order x_s(0), x_e(0), x_s(1), x_e(1), ...:
  !!
  !!   sa_1(p)*x_e(p-1) + x_s(p) + sc_1(p)*x_e(p)   = d_s(p)
  !!   sa_n(p)*x_s(p)   + x_e(p) + sc_n(p)*x_s(p+1) = d_e(p)
  !!
  !! where sa_1, sc_1, sa_n and sc_n are rank p's dist_sa(1), dist_sc(1),
  !! dist_sa(n) and dist_sc(n), and d_s and d_e its first and last values
  !! after the first phase. The system is cyclic: its corners couple the last
  !! rank to the first. On a non-periodic line they vanish, because the
  !! boundary rows have no outside neighbour.
  !!
  !! The default second phase (der_univ_subs) solves this system only
  !! approximately: it drops sc_1 and sa_n and pairs each block boundary with
  !! its neighbour alone. The dropped entries decay geometrically with the
  !! points per rank, at a rate set by how diagonally dominant the system is.
  !! That is invisible for the derivative operators, but not for the low-pass
  !! filter near its sharp limit filter_alpha -> 0.5, whose gain it then gets
  !! wrong by about sa_n even at zero wavenumber, every time it is applied.
  !!
  !! Here every rank gathers the reduced rows of its lines and solves the
  !! 2*nproc system exactly, keeping its own x_s and x_e. The system is small,
  !! so solving it on every rank costs less than more rounds of messages.
  use m_common, only: dp, MPI_X3D2_DP
  use m_logger, only: abort_run
  use m_mesh_content, only: par_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Allgather, &
                   MPI_Comm_split

  implicit none

  private
  public :: solve_reduced_exact

  !> Communicator along each direction, made on first use
  integer, save :: dir_comm(3) = MPI_COMM_NULL
  !> The rank layout each communicator was made for
  integer, save :: comm_layout(6, 3) = -1

contains

  subroutine solve_reduced_exact(x_s, x_e, d_s, d_e, sa_1, sc_1, sa_n, &
                                 sc_n, par, dir)
    !! Solve the reduced system of every line exactly. Line l of d_s and d_e
    !! must be the same physical line on every rank along dir, as it is for
    !! the neighbour exchange of the default second phase.
    real(dp), intent(out) :: x_s(:), x_e(:)
    real(dp), intent(in) :: d_s(:), d_e(:)
    !> This rank's dist_sa(1), dist_sc(1), dist_sa(n) and dist_sc(n)
    real(dp), intent(in) :: sa_1, sc_1, sa_n, sc_n
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir

    real(dp), allocatable :: coeffs(:), rows(:), all_rows(:)
    real(dp), allocatable :: lo(:), up(:), cp(:), den(:), z(:), y(:)
    real(dp) :: corner_tr, corner_bl, gamma, fact, zfac
    integer :: np, me, nl, m, l, i, p, comm, ierr

    np = par%nproc_dir(dir)
    me = par%nrank_dir(dir)
    nl = size(d_s)
    m = 2*np
    if (np < 2) &
      call abort_run('solve_reduced_exact needs a decomposed direction.')

    comm = direction_comm(par, dir)

    ! The coefficients can differ between ranks (boundary closures), so
    ! every rank's are gathered, then the rows of every line.
    allocate (coeffs(4*np), rows(2*nl), all_rows(2*nl*np))
    call MPI_Allgather([sa_1, sc_1, sa_n, sc_n], 4, MPI_X3D2_DP, &
                       coeffs, 4, MPI_X3D2_DP, comm, ierr)
    rows(1:nl) = d_s
    rows(nl + 1:2*nl) = d_e
    call MPI_Allgather(rows, 2*nl, MPI_X3D2_DP, all_rows, 2*nl, &
                       MPI_X3D2_DP, comm, ierr)

    ! Row i: lo(i)*x(i-1) + x(i) + up(i)*x(i+1), with x(0) = x(m) and
    ! x(m+1) = x(1) closing the cycle
    allocate (lo(m), up(m))
    do p = 0, np - 1
      lo(2*p + 1) = coeffs(4*p + 1)
      up(2*p + 1) = coeffs(4*p + 2)
      lo(2*p + 2) = coeffs(4*p + 3)
      up(2*p + 2) = coeffs(4*p + 4)
    end do
    corner_tr = lo(1)
    corner_bl = up(m)

    ! Cyclic tridiagonal by Sherman-Morrison: solve the system with the
    ! corners folded into the first and last diagonal entries, then correct
    ! along the solution z for u = (gamma, 0, ..., 0, corner_bl). The unit
    ! diagonal makes gamma = -1 safe. The factorisation is shared by every
    ! line.
    gamma = -1._dp
    allocate (cp(m), den(m), z(m), y(m))
    den(1) = 1._dp - gamma
    cp(1) = up(1)/den(1)
    do i = 2, m
      den(i) = 1._dp - lo(i)*cp(i - 1)
      if (i == m) den(i) = den(i) - corner_bl*corner_tr/gamma
      cp(i) = up(i)/den(i)
    end do
    z = 0._dp
    z(1) = gamma
    z(m) = corner_bl
    call thomas(z)
    zfac = 1._dp + z(1) + corner_tr*z(m)/gamma

    do l = 1, nl
      do p = 0, np - 1
        y(2*p + 1) = all_rows(2*nl*p + l)
        y(2*p + 2) = all_rows(2*nl*p + nl + l)
      end do
      call thomas(y)
      fact = (y(1) + corner_tr*y(m)/gamma)/zfac
      x_s(l) = y(2*me + 1) - fact*z(2*me + 1)
      x_e(l) = y(2*me + 2) - fact*z(2*me + 2)
    end do

  contains

    subroutine thomas(r)
      !! In place solve with the shared factorisation
      real(dp), intent(inout) :: r(:)
      integer :: k

      r(1) = r(1)/den(1)
      do k = 2, m
        r(k) = (r(k) - lo(k)*r(k - 1))/den(k)
      end do
      do k = m - 1, 1, -1
        r(k) = r(k) - cp(k)*r(k + 1)
      end do
    end subroutine thomas

  end subroutine solve_reduced_exact

  integer function direction_comm(par, dir) result(comm)
    !! The ranks sharing this rank's position in the other two directions,
    !! ordered along dir. Made collectively on first use, and again if the
    !! rank layout changes (a test building a second mesh, say).
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir

    integer :: layout(6), d1, d2, color, ierr

    layout = [par%nrank_dir, par%nproc_dir]
    if (dir_comm(dir) == MPI_COMM_NULL &
        .or. any(comm_layout(:, dir) /= layout)) then
      d1 = modulo(dir, 3) + 1
      d2 = modulo(dir + 1, 3) + 1
      color = par%nrank_dir(d1) + par%nproc_dir(d1)*par%nrank_dir(d2)
      ! A communicator made for an earlier layout is left as it is: freeing
      ! it is not worth a collective here.
      call MPI_Comm_split(MPI_COMM_WORLD, color, par%nrank_dir(dir), &
                          dir_comm(dir), ierr)
      comm_layout(:, dir) = layout
    end if
    comm = dir_comm(dir)
  end function direction_comm

end module m_dist_reduced
