submodule(m_decomp) m_decomp_dummy

  use mpi
  implicit none

  type, extends(decomp_t) :: decomp_dummy_t
  contains
    procedure :: decomposition => decomposition_dummy
  end type

  contains

  module subroutine init_decomp(decomp)
    class(decomp_t), allocatable, intent(out):: decomp

    allocate(decomp_dummy_t :: decomp)
    decomp%is_avail_2decomp = .false.
  end subroutine

  module subroutine decomposition_dummy(self, grid, par)
    use m_mesh_content, only: par_t, grid_t
    class(decomp_generic_t) :: self
    class(grid_t), intent(inout) :: grid

     error stop "This build doesn't support 2decomp decomposition"
  end subroutine decomposition_dummy

end submodule