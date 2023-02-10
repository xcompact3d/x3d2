program test_allocator_main
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t

  implicit none

  interface
     subroutine test_allocator(alloc, allp)
       import allocator_t
       class(allocator_t), intent(in) :: alloc
       logical, intent(in) :: allp
     end subroutine test_allocator
  end interface

  logical :: allpass
  integer, parameter :: dims(3) = [8, 8, 8]
  class(allocator_t), allocatable :: allocator


  allocator = allocator_t(dims)

  call test_allocator(allocator, allpass)

contains



end program test_allocator_main
