module m_poisson_cg_backend

  use m_base_poisson_cg, only: poisson_solver_t, poisson_precon_t, laplace_operator_t

  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: init_solver
  public :: init_precon

  type, extends(poisson_solver_t), public :: poisson_solver_impl
    private
    type(laplace_operator_t) :: lapl ! The high-order Laplacian operator
    class(poisson_precon_t), allocatable :: precon ! The (low-order) preconditioner
  contains
    private
    ! Solves the Poisson problem
    procedure, public :: solve => solve
  end type poisson_solver_impl

  type, extends(poisson_precon_t), public :: poisson_precon_impl
    !! Wrapper definition of the Poisson preconditioner.
    !! User code should use this class which will instantiate the
    !! backend-specific code (determined at compile time).
    private
  contains
    private
    ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
    procedure, public :: apply => apply_precon
  end type poisson_precon_impl

  interface poisson_precon_t
    !! Public constructor for the Poisson preconditioner object.
    module procedure init_precon
  end interface poisson_precon_t

contains

  module subroutine init_solver(solver, backend, mesh)
    class(poisson_solver_t), allocatable, intent(out) :: solver
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), intent(in) :: mesh

    allocate (poisson_solver_impl :: solver)

    select type (solver)
    type is (poisson_solver_impl)
      solver%precon = init_precon(backend)
      solver%lapl = laplace_operator_t(backend, mesh)
    class default
      error stop "Dummy CG solver init failed"
    end select

    error stop "This dummy module does not implement CG, recompile the code with PETSc"
  end subroutine init_solver

  subroutine solve(self, p, f, backend)
    class(poisson_solver_impl) :: self
    class(field_t), intent(inout) :: p
    class(field_t), intent(in) :: f
    class(base_backend_t), intent(in) :: backend

    ! Call the backend-specific subroutine
    error stop "Dummy CG backend called"
  end subroutine solve

  function init_precon(backend) result(precon)
    !! Public constructor for the Poisson preconditioner object.
    class(base_backend_t), target, intent(in) :: backend
    class(poisson_precon_t), allocatable :: precon

    allocate (poisson_precon_impl :: precon)

  end function init_precon

  subroutine apply_precon(self, p, b, backend)
    ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
    class(poisson_precon_impl) :: self
    class(field_t), intent(in) :: p    ! Pressure solution
    class(field_t), intent(inout) :: b ! The evaluated matrix-vector product
    class(base_backend_t), intent(in) :: backend

    ! Call the backend-specific subroutine
    error stop "Dummy CG backend called"
  end subroutine apply_precon

end module m_poisson_cg_backend
