module m_vector3d
  use m_stencil, only: stencil
  use m_tridiagsolv, only: tridiagsolv

  implicit none


  type, abstract :: vector3d
     character(len=:), allocatable :: name
     integer :: rankid, nranks
     type(stencil) :: bulk_stencil(2)
     type(stencil) :: left_stencils(2, 2)
     type(stencil) :: right_stencils(2, 2)
     class(tridiagsolv), pointer :: thomas_solver
   contains
     procedure(field_op), deferred :: transport
     procedure(field_op), deferred :: div
     procedure, public :: u => get_component_ptr
  end type vector3d

  abstract interface
     subroutine field_op(self, rslt)
       import :: vector3d
       class(vector3d), intent(in) :: self
       class(vector3d), intent(inout) :: rslt
     end subroutine field_op
  end interface

contains

    function get_component_ptr(self, i) result(ptr)
    class(vector3d), intent(in) :: self
    integer, intent(in) :: i
    real, pointer :: ptr(:, :, :)

    select case(i)
    case (1)
       ptr => self%u1%data
    case (2)
       ptr => self%u2%data
    case (3)
       ptr => self%u3%data
    end select
  end function get_component_ptr
end module m_vector3d
