module m_scalar
  type, abstract :: scalar
     character(len=:), allocatable :: name
     integer :: dim(3)
     integer :: nproc, rankid
   contains
     procedure, deferred :: transport
     procedure, deferred :: grad
  end type scalar

  abstract interface
     subroutine field_op(self, rslt)
       import :: scalar
       class(scalar), intent(in) :: self
       class(scalar), intent(inout) :: rslt
     end subroutine field_op
  end interface
end module m_scalar
