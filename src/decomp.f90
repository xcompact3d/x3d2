module m_decomp
  use m_mesh_content, only: par_t, grid_t
  implicit none

  ! Object implementing the actual decomposition
  type, abstract :: decomp_t 
  contains
    procedure(decomposition), public, deferred :: decomposition
  end type decomp_t

  ! Abstraction layer to allow different version of the decomposition computation
  type, public :: decomp_mod_t
    class(decomp_t), allocatable :: decomp
    contains
    procedure, public :: decomp_grid
  end type

  interface decomp_mod_t
    module procedure init_decomp_mod
  end interface
 
  interface
    module subroutine init_decomp(decomp)
      class(decomp_t), allocatable, intent(out) :: decomp
    end subroutine

    module subroutine decomposition(self, grid, par)
      use m_mesh_content, only: par_t, grid_t
      class(decomp_t) :: self
      class(grid_t), intent(inout) :: grid
      class(par_t), intent(inout) :: par
    end subroutine

   end interface
 
  contains

  function init_decomp_mod() result(decomp)
    type(decomp_mod_t) :: decomp

    call init_decomp(decomp%decomp)

  end function
 
  module subroutine decomp_grid(self, grid, par)
      class(decomp_mod_t) :: self
      class(grid_t), intent(inout) :: grid
      class(par_t), intent(inout) :: par

      call self%decomp%decomposition(grid, par)
  end subroutine


end module m_decomp