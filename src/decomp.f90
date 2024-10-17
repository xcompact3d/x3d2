module m_decomp
implicit none

  type, abstract :: decomp_t 
  contains
  procedure(decomposition), public, deferred :: decomposition
  end type decomp_t

  interface
    subroutine decomposition(self, grid, par)
      use m_grid, only: grid_t
      use m_par, only: par_t
      import :: decomp_t
      class(decomp_t) :: self
      class(grid_t), intent(inout) :: grid
      class(par_t), intent(inout) :: par
    end subroutine
  end interface
  
  contains

  module subroutine test()
    print *, "test"
  end subroutine

end module m_decomp