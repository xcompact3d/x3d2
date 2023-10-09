module m_tridiagsolv
  !! This module provides access to derived types containing data and
  !! functionality to solve triadiagonal or quasi-triadiagonal system
  !! of equations.
  !!
  !! Derived type [[m_tridiagsolv(module):tridiagsolv(type)]] is used
  !! to solve a pure triadiagonal system for a specific set of
  !! coefficients.  The type is instantiated by passing the lower,
  !! upper and central diagonals.  The system can then been solved for
  !! arbitrary right-hand vectors by colling the
  !! [[m_tridiagsolv(module):solve(subroutine)]] type-bound procedure.
  !!
  !! Derived type [[m_tridiagsolv(module):periodic_tridiagsolv(type)]]
  !! plays a role similar to
  !! [[m_tridiagsolv(module):tridiagsolv(type)]] for
  !! quasi-triadiagonal system typical of configurations with periodic
  !! boundary conditions.
   implicit none

   type :: tridiagsolv
     !! Tridiagonal solver for a particular tridiagonal coefficient
     !! matrix.  The solution is computed using the TMDA (or Thomas)
     !! algorithm, see
     !! [(Wikipedia)Tridiagonal matrix algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).
     !!
     !! A solver can be created by passing the lower, central and
     !! upper diagonals to the type constructor, see
     !! [[m_tridiagsolv(module):construct(function)]].
     !!
     !! ```f90
     !! real :: low(7), up(7)
     !! real :: diag(8)
     !!
     !! diag = 1.
     !! call random_number(low)
     !! call random_number(up)
     !!
     !! thomas_solver = tridiagsolv(low, up, diag)
     !! ```
     !!
     !! The tridiagonal system can then be solved for a collection of
     !! M arbitrary right-hand-side vectors using the `solve`
     !! type-bound procedure.  This collection is passed as a
     !! two-dimensional array of shape `M * size(diag)`, that is the
     !! right-hand side vectors represented along the second
     !! dimension.
     !!
     !! ```f90
     !! ! Solve for M right-hand side vectors at once.
     !! integer, parameter :: M = 4
     !! real :: rhs(M, 8), solution(M, 8)
     !! call random_number(rhs)
     !!
     !! thomas_solver%solve(rhs, solution)
     !! ```
      private
      real, allocatable :: fwd(:), bwd(:), updiag(:)
   contains
      procedure, public :: solve
   end type tridiagsolv

   type, extends(tridiagsolv) :: periodic_tridiagsolv
     !! Solver for a quasi-tridiagonal tridiagonal coefficient matrix
     !! originating from periodic boundary conditions.
     !!
     !! ```
     !! b1  c1  0   0    0   a1
     !! a2  b2  c2  0    0   0
     !! 0   a3  b3  c3   0   0
     !! 0   0   a4  b4  c4   0
     !! 0   0   0   a5  b5   c5
     !! c6   0   0  0   a6   b6
     !! ```
     !!
     !! The solution is computed using a modified TMDA (or Thomas)
     !! algorithm using the Shermann-Morrison formula, see
     !! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variant
     !!
     !! A solver is created in the same way as the full tridiagonal
     !! case, by passing the lower, central and upper diagonals to the
     !! type constructor,
     !! see [[m_tridiagsolv(module):construct_periodic(function)]].
     !!
     !!  ```
     !!  real :: low(7), up(7)
     !!  real :: diag(8)
     !!
     !!  diag = 1.
     !!  call random_number(low)
     !!  call random_number(up)
     !!
     !!  thomas_solver = tridiagsolv(low, up, diag)
     !! ```
      real :: gamma = -1.
      !> Helper vector, solution to intermediate system A' * q = u
      real, allocatable :: q(:)
   contains
      procedure, public :: solve => solve_periodic
   end type periodic_tridiagsolv

   interface tridiagsolv
      module procedure construct
   end interface tridiagsolv

   interface periodic_tridiagsolv
      module procedure construct_periodic
   end interface periodic_tridiagsolv

contains

   function construct(low, up, diag) result(solver)
      type(tridiagsolv) :: solver
      real, intent(in) :: low(:), up(:)
      real, optional, intent(in) :: diag(:)
      real, allocatable :: fwd(:), i_bwd(:)
      integer :: i, n

      n = size(low) + 1
      ! Allocate solver's coefficient arrays.
      allocate (solver%fwd(n), solver%bwd(n))
      solver%updiag = up
      ! Allocate extra temp arrays for readability.  These only live
      ! within this the function scop.
      allocate (fwd(n), i_bwd(n))

      ! Initialise bwd array to diagonal.  The default behavior is to
      ! work with a tridiag array with 1s on the main diagonal.
      if (present(diag)) then
         i_bwd = diag
      else
         i_bwd = 1.
      end if

      !i_bwd is the /inverse/ of the coefficients for the bwd step.
      do i = 2, n
         fwd(i) = low(i - 1) / i_bwd(i - 1)
         i_bwd(i) = i_bwd(i) - up(i - 1) * fwd(i)
      end do
      solver%fwd = fwd
      solver%bwd = 1./i_bwd
   end function construct

   function construct_periodic(low, up) result(solver)
      real, intent(in) :: low(:), up(:)
      type(periodic_tridiagsolv) :: solver

      real, allocatable :: u(:), u_vec(:, :), q_vec(:, :)
      real, allocatable :: diag(:)
      integer :: i, n
      n = size(low) + 1
      ! Runs non-periodic (base type) constructor which sets up forward
      ! and backward coefficient arrays.
      diag = [2., &
           & (1., i=2, n - 1), &
           & 1.+up(1) * up(1) &
           & ]
      solver%tridiagsolv = tridiagsolv(low, up, diag=diag)
      ! q member array is used in the Sermann-Morrison formula

      ! solve() method expects a rank 2 array of shape (SZ X n) as
      ! input, where SZ is the number of pencils in the pencil group.
      ! Therefore artificially work a temp array u_vec and q_vec of rank
      ! 2.
      allocate (u(n), u_vec(1, n), q_vec(1, n))
      u = [solver%gamma, (0., i=2, n - 1), up(1)]
      u_vec = transpose(reshape(u, [n, 1], pad=u))
      call solver%tridiagsolv%solve(u_vec, q_vec)
      solver%q = q_vec(1, :)
   end function construct_periodic

   pure subroutine solve(self, f, df)
    !! Solve the tridiagonal system associated to the
    !! [[m_tridiagsolv(module):tridiagsolv(type)]] instance.
    !!
    !! The system is solved for a collection of M right-hand side
    !! vectors, arranged into a two-dimensional array of shape M X N,
    !! where N is the size of the system.
      class(tridiagsolv), intent(in) :: self
      real, intent(in) :: f(:, :)
      !> Solutions
      real, intent(out) :: df(:, :)
      integer :: i, j, n

      ! TODO: Overwrite f in place and get rid of df array.
      df = f
      n = size(self%fwd)
      if (size(f, 2) /= n) error stop
      do j = 2, n
         do i = 1, size(f, 1)
            df(i, j) = df(i, j) - df(i, j - 1) * self%fwd(j)
         end do
      end do
      do i = 1, size(f, 1)
         df(i, n) = df(i, n) * self%bwd(n)
      end do

      do j = n - 1, 1, -1
         do i = 1, size(f, 1)
            df(i, j) = (df(i, j) - df(i, j + 1) * self%updiag(j)) &
                      & * self%bwd(j)
         end do
      end do
   end subroutine solve

   pure subroutine solve_periodic(self, f, df)
    !! Solve the tridiagonal system associated to the
    !! [[m_tridiagsolv(module):periodic_tridiagsolv(type)]]
    !! instance.
    !!
    !! The system is solved for a collection of M right-hand side
    !! vectors, arranged into a two-dimensional array of shape M X N,
    !! where N is the size of the system.
      class(periodic_tridiagsolv), intent(in) :: self
      real, intent(in) :: f(:, :)
      real, intent(out) :: df(:, :)
      real :: y(size(f, 1), size(f, 2) - 1)
      integer :: i, m
      real :: alpha

      ! In the periodic case the first element of the differentiation
      ! pencil (f(1)) is the same as the last (f(n)).  Alhtough we have
      ! n data, the effective domain size if n - 1.
      m = size(f, 2) - 1
      call self%tridiagsolv%solve(f(:, 1:m), y)

      alpha = self%updiag(1)
      select type (self)
      type is (periodic_tridiagsolv)
         do i = 1, size(f, 1)
            df(i, 1:m) = y(i, :) - ((y(i, 1) - alpha * y(i, m)) &
                 & / (1.+self%q(1) - alpha * self%q(m))) * self%q
         end do
      class default
         error stop
      end select
    !! And don't forget to enforce periocity
      df(:, size(df, 2)) = df(:, 1)
   end subroutine solve_periodic
end module m_tridiagsolv
