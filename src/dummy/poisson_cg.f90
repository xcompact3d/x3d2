!! Dummy implementatino of the iterative Poisson solver

submodule(m_poisson_cg) m_dummy_poisson_cg

contains
  
  module subroutine init_precon_impl(precon, backend)
    class(poisson_precon_impl_t), allocatable, intent(out) :: precon
    class(base_backend_t), intent(in) :: backend

    error stop "This dummy module does not implement CG, recompile the code with PETSc"

  end subroutine init_precon_impl

  module subroutine init_solver(solver, backend)
    class(poisson_solver_t), allocatable, intent(out) :: solver
    class(base_backend_t), target, intent(in) :: backend

    error stop "This dummy module does not implement CG, recompile the code with PETSc"

  end subroutine init_solver
  
end submodule m_dummy_poisson_cg
