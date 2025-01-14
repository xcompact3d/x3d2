module m_decomp

  implicit none

  !public :: is_avail_decomp
  !public :: decomposition

  contains

  module function is_avail_2decomp() result(avail)
    logical :: avail

    avail = .false.
  end function

  module subroutine decomposition_2decomp(grid, par)
    use m_mesh_content, only: par_t, grid_t
    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par 

     error stop "This build doesn't support 2decomp decomposition"
  end subroutine decomposition_dummy

end module