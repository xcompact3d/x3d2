module m_vector3d
   use m_stencil, only: stencil
   use m_tridiagsolv, only: tridiagsolv
   use m_allocator, only: allocator_t, memblock_t

   implicit none

   type, abstract :: vector3d
      character(len=:), allocatable :: name
      integer :: rankid, nranks
      class(allocator_t), pointer :: allocator
      class(memblock_t), pointer :: u1, u2, u3
      real :: xnu
   contains
      procedure(field_op), deferred :: transport, div
      procedure(transport_op), deferred :: transport_dir
      procedure, public :: u => get_component_ptr
   end type vector3d

   abstract interface
      function field_op(self)
         import :: vector3d
         class(vector3d), intent(in) :: self
         class(vector3d), allocatable :: field_op
      end function field_op
      subroutine transport_op(self, dim, rslt)
         import :: vector3d
         class(vector3d), intent(in) :: self
         integer, intent(in) :: dim
         real, intent(out) :: rslt(:, :, :)
      end subroutine transport_op
   end interface

contains

   function get_component_ptr(self, i) result(ptr)
      class(vector3d), intent(in) :: self
      integer, intent(in) :: i
      real, pointer :: ptr(:, :, :)

      select case (i)
      case (1)
         ptr => self%u1%data
      case (2)
         ptr => self%u2%data
      case (3)
         ptr => self%u3%data
      end select
   end function get_component_ptr
end module m_vector3d
