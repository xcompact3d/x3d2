module m_poisson_cg
  !! Module defining a Poisson solver based on the (preconditioned) Conjugate Gradient method.

  use m_allocator, only: field_t
  use m_tdsops, only: dirps_t
  
  implicit none

  private

  public :: poissmult
  
  type, abstract, public :: poisson_cg_t
     !! Conjugate Gradient based Poisson solver.
     !! Supports any decomposition that is also supported by the underlying finite difference
     !! schemes.

     ! Prevent default access to components of type.
     private 
   contains
  end type poisson_cg_t
  
  interface
     module function init_cg(xdirps, ydirps, zdirps) result(poisson_cg)
       class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
       class(poisson_cg_t), allocatable :: poisson_cg
     end function init_cg
  end interface
  
  interface poisson_cg_t
     !! Public constructor for the poisson_cg_t type.
     module procedure init_cg
  end interface poisson_cg_t
  
contains

  subroutine poissmult(x, f)
    !! Computes the action of the Poisson operator, i.e. `f = Mx` where `M` is the discrete
    !! Laplacian.
    class(field_t), intent(in) :: x    ! The input field
    class(field_t), intent(inout) :: f ! The output field

    ! Compute lapl(x) => f
  end subroutine poissmult
  
end module m_poisson_cg
