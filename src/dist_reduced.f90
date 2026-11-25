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
  !! factorise_reduced gathers the coefficients and factorises the system
  !! once per operator, when it is set up; solve_reduced_exact gathers the
  !! rows and solves them on the host.
  use m_common, only: dp, MPI_X3D2_DP
  use m_logger, only: abort_run
  use m_mesh_content, only: par_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Allgather, &
                   MPI_Comm_free, MPI_Comm_split

  implicit none

  private
  public :: reduced_lu_t, factorise_reduced, solve_reduced_exact, &
            direction_comm

  type :: reduced_lu_t
    !! Factorisation of one operator's reduced system along a decomposed
    !! direction. It depends only on the coefficients of the ranks along the
    !! direction, so it is made once and shared by every line and every call.
    !!
    !! The cyclic system is solved by Sherman-Morrison: with the corners
    !! folded into the first and last diagonal entries it is tridiagonal,
    !! with sub-diagonal lo and Thomas factors den and cp, and its solution y
    !! of a right-hand side is corrected along z:
    !!
    !!   x = y - (y(1) + corner_g*y(m))/zfac*z
    integer :: m = 0 !! Size of the system, 2*nproc along the direction
    integer :: me = -1 !! This rank's position along the direction
    real(dp), allocatable :: lo(:), cp(:), den(:), z(:)
    real(dp) :: corner_g = 0._dp, zfac = 1._dp
  end type reduced_lu_t

  !> Communicator along each direction, made on first use
  integer, save :: dir_comm(3) = MPI_COMM_NULL
  !> The rank layout each communicator was made for
  integer, save :: comm_layout(6, 3) = -1

contains

  subroutine factorise_reduced(lu, sa_1, sc_1, sa_n, sc_n, par, dir)
    !! Gather the coefficients of every rank along dir and factorise the
    !! reduced system. Collective over the ranks along dir.
    type(reduced_lu_t), intent(out) :: lu
    !> This rank's dist_sa(1), dist_sc(1), dist_sa(n) and dist_sc(n)
    real(dp), intent(in) :: sa_1, sc_1, sa_n, sc_n
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir

    !> The unit diagonal makes gamma = -1 safe
    real(dp), parameter :: gamma = -1._dp
    real(dp), allocatable :: coeffs(:), up(:)
    real(dp) :: corner_tr, corner_bl
    integer :: np, m, i, p, ierr

    np = par%nproc_dir(dir)
    m = 2*np
    if (np < 2) &
      call abort_run('factorise_reduced needs a decomposed direction.')

    ! The coefficients can differ between ranks (boundary closures), so
    ! every rank's are gathered.
    allocate (coeffs(4*np))
    call MPI_Allgather([sa_1, sc_1, sa_n, sc_n], 4, MPI_X3D2_DP, &
                       coeffs, 4, MPI_X3D2_DP, direction_comm(par, dir), &
                       ierr)

    ! Row i: lo(i)*x(i-1) + x(i) + up(i)*x(i+1), with x(0) = x(m) and
    ! x(m+1) = x(1) closing the cycle
    lu%m = m
    lu%me = par%nrank_dir(dir)
    allocate (lu%lo(m), up(m))
    do p = 0, np - 1
      lu%lo(2*p + 1) = coeffs(4*p + 1)
      up(2*p + 1) = coeffs(4*p + 2)
      lu%lo(2*p + 2) = coeffs(4*p + 3)
      up(2*p + 2) = coeffs(4*p + 4)
    end do
    corner_tr = lu%lo(1)
    corner_bl = up(m)

    ! Sherman-Morrison with u = (gamma, 0, ..., 0, corner_bl): fold the
    ! corners into the first and last diagonal entries, then z solves the
    ! folded system for u.
    allocate (lu%cp(m), lu%den(m), lu%z(m))
    lu%den(1) = 1._dp - gamma
    lu%cp(1) = up(1)/lu%den(1)
    ! The corner correction of the last entry stays out of the loop: with it
    ! inside, nvfortran 25.3 at -O2 breaks the recurrence through cp(i - 1),
    ! and double precision builds solved a different system.
    do i = 2, m
      lu%den(i) = 1._dp - lu%lo(i)*lu%cp(i - 1)
      lu%cp(i) = up(i)/lu%den(i)
    end do
    lu%den(m) = lu%den(m) - corner_bl*corner_tr/gamma
    lu%cp(m) = up(m)/lu%den(m)
    lu%z = 0._dp
    lu%z(1) = gamma
    lu%z(m) = corner_bl
    call thomas_lu(lu, lu%z)
    lu%corner_g = corner_tr/gamma
    lu%zfac = 1._dp + lu%z(1) + lu%corner_g*lu%z(m)

  end subroutine factorise_reduced

  subroutine solve_reduced_exact(x_s, x_e, d_s, d_e, lu, par, dir)
    !! Solve the reduced system of every line exactly on the host. Line l of
    !! d_s and d_e must be the same physical line on every rank along dir, as
    !! it is for the neighbour exchange of the default second phase.
    real(dp), intent(out) :: x_s(:), x_e(:)
    real(dp), intent(in) :: d_s(:), d_e(:)
    !> The factorisation from factorise_reduced for this operator
    type(reduced_lu_t), intent(in) :: lu
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir

    real(dp), allocatable :: rows(:), all_rows(:), y(:)
    real(dp) :: fact
    integer :: np, me, nl, m, l, p, ierr

    np = par%nproc_dir(dir)
    me = lu%me
    m = lu%m
    nl = size(d_s)
    if (m /= 2*np .or. me /= par%nrank_dir(dir)) &
      call abort_run('solve_reduced_exact: the reduced system was not '// &
                     'factorised for this rank layout.')

    allocate (rows(2*nl), all_rows(2*nl*np))
    rows(1:nl) = d_s
    rows(nl + 1:2*nl) = d_e
    call MPI_Allgather(rows, 2*nl, MPI_X3D2_DP, all_rows, 2*nl, &
                       MPI_X3D2_DP, direction_comm(par, dir), ierr)

    allocate (y(m))
    do l = 1, nl
      do p = 0, np - 1
        y(2*p + 1) = all_rows(2*nl*p + l)
        y(2*p + 2) = all_rows(2*nl*p + nl + l)
      end do
      call thomas_lu(lu, y)
      fact = (y(1) + lu%corner_g*y(m))/lu%zfac
      x_s(l) = y(2*me + 1) - fact*lu%z(2*me + 1)
      x_e(l) = y(2*me + 2) - fact*lu%z(2*me + 2)
    end do

  end subroutine solve_reduced_exact

  subroutine thomas_lu(lu, r)
    !! In place solve of the folded tridiagonal system
    type(reduced_lu_t), intent(in) :: lu
    real(dp), intent(inout) :: r(:)
    integer :: k

    r(1) = r(1)/lu%den(1)
    do k = 2, lu%m
      r(k) = (r(k) - lu%lo(k)*r(k - 1))/lu%den(k)
    end do
    do k = lu%m - 1, 1, -1
      r(k) = r(k) - lu%cp(k)*r(k + 1)
    end do
  end subroutine thomas_lu

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
      ! Every rank reaches here together, as the split is collective, so the
      ! communicator made for an earlier layout can be freed first.
      if (dir_comm(dir) /= MPI_COMM_NULL) &
        call MPI_Comm_free(dir_comm(dir), ierr)
      call MPI_Comm_split(MPI_COMM_WORLD, color, par%nrank_dir(dir), &
                          dir_comm(dir), ierr)
      comm_layout(:, dir) = layout
    end if
    comm = dir_comm(dir)
  end function direction_comm

end module m_dist_reduced
