module m_common
   implicit none

   integer, parameter :: dp=kind(0.0d0)
   real(dp), parameter :: pi = 4*atan(1.0_dp)

   type :: derps_t
      integer :: n
   end type derps_t

end module m_common
