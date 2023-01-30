module m_vector3d
  type, abstract :: vector3d
     character(len=:), allocatable :: name
     integer :: rankid, nranks
   contains
     procedure(field_op), deferred :: transport
     procedure(field_op), deferred :: div
  end type vector3d

  abstract interface
     subroutine field_op(self, rslt)
       import :: vector3d
       class(vector3d), intent(in) :: self
       class(vector3d), intent(inout) :: rslt
     end subroutine field_op
  end interface
end module m_vector3d
