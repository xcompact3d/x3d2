module m_poisson_cg
  !! Module defining a Poisson solver based on the (preconditioned) Conjugate
  !! Gradient method.

  use m_base_poisson_cg, only: poisson_solver_t
  use m_poisson_cg_backend, only: init_solver

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t

  implicit none

  private

  type, public :: poisson_cg_t
    !! Conjugate Gradient based Poisson solver.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.

    ! Prevent default access to components of type.
    private
    class(poisson_solver_t), allocatable :: solver

  contains
    private
    ! Solves the Poisson problem
    procedure, public :: solve
  end type poisson_cg_t

  interface poisson_cg_t
    !! Public constructor for the poisson_cg_t type.
    module procedure init_cg
  end interface poisson_cg_t

contains

  subroutine solve(self, p, f, backend)
    ! Solves the Poisson problem
    class(poisson_cg_t) :: self
    class(field_t), intent(inout) :: p ! Pressure solution
    class(field_t), intent(in) :: f    ! Poisson RHS
    class(base_backend_t), intent(in) :: backend

    call self%solver%solve(p, f, backend)
  end subroutine solve

  function init_cg(backend, mesh) result(solver)
    !! Initialises the conjugate gradient (CG) solver.
    !! XXX: The solver implementation is responsible for initialising the
    !!      preconditioner, i.e. it should at some point during its initialisation
    !!      do the equivalent of: `call self%precon = poisson_precon_t(backend)`.
    class(base_backend_t), target, intent(in) :: backend
    type(poisson_cg_t) :: solver
    type(mesh_t), intent(in) :: mesh

    call init_solver(solver%solver, backend, mesh)

  end function init_cg

end module m_poisson_cg
