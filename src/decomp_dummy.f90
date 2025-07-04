module m_decomp
  !! Dummy implementation of the decomposition to be linked against when 2decomp&fft isn't available.

  implicit none

contains

  function is_avail_2decomp() result(avail)
    logical :: avail

    avail = .false.
  end function

  subroutine decomposition_2decomp(grid, par)
    use m_mesh_content, only: par_t, grid_t
    class(grid_t), intent(inout) :: grid
    class(par_t), intent(inout) :: par

    error stop "This build doesn't support 2decomp decomposition"
  end subroutine decomposition_2decomp

end module
