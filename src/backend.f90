module m_base_backend
   use m_allocator, only: allocator_t, field_t
   use m_common, only: dp
   use m_tdsops, only: tdsops_t, dirps_t

   implicit none

   type, abstract :: base_backend_t
      !! base_backend class defines all the abstract operations that the
      !! solver class requires.
      !!
      !! For example, transport equation in solver class evaluates the
      !! derivatives in x, y, and z directions, and reorders the input
      !! fields as required. Then finally, combines all the directional
      !! derivatives to obtain the divergence of U*.
      !!
      !! All these high level operations solver class executes are
      !! defined here using the abstract interfaces. Every backend
      !! implementation extends the present abstact backend class to
      !! define the specifics of these operations based on the target
      !! architecture.

      real(dp) :: nu
      integer :: nx_loc, ny_loc, nz_loc
      class(allocator_t), pointer :: allocator
      class(dirps_t), pointer :: xdirps, ydirps, zdirps
   contains
      procedure(transeq_ders), deferred :: transeq_x
      procedure(transeq_ders), deferred :: transeq_y
      procedure(transeq_ders), deferred :: transeq_z
      procedure(tds_solve), deferred :: tds_solve
      procedure(reorder), deferred :: reorder
      procedure(sum_intox), deferred :: sum_yintox
      procedure(sum_intox), deferred :: sum_zintox
      procedure(vecadd), deferred :: vecadd
      procedure(get_fields), deferred :: get_fields
      procedure(set_fields), deferred :: set_fields
      procedure(alloc_tdsops), deferred :: alloc_tdsops
   end type base_backend_t

   abstract interface
      subroutine transeq_ders(self, du, dv, dw, u, v, w, dirps)
         !! transeq equation obtains the derivatives direction by
         !! direction, and the exact algorithm used to obtain these
         !! derivatives are decided at runtime. Backend implementations
         !! are responsible from directing calls to transeq_ders into
         !! the correct algorithm.
         import :: base_backend_t
         import :: field_t
         import :: dirps_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: du, dv, dw
         class(field_t), intent(in) :: u, v, w
         type(dirps_t), intent(in) :: dirps
      end subroutine transeq_ders
   end interface

   abstract interface
      subroutine tds_solve(self, du, u, dirps, tdsops)
         !! transeq equation obtains the derivatives direction by
         !! direction, and the exact algorithm used to obtain these
         !! derivatives are decided at runtime. Backend implementations
         !! are responsible from directing calls to transeq_ders into
         !! the correct algorithm.
         import :: base_backend_t
         import :: field_t
         import :: dirps_t
         import :: tdsops_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: du
         class(field_t), intent(in) :: u
         type(dirps_t), intent(in) :: dirps
         class(tdsops_t), intent(in) :: tdsops
      end subroutine tds_solve
   end interface

   abstract interface
      subroutine reorder(self, u_, u, direction)
         !! reorder subroutines are straightforward, they rearrange
         !! data into our specialist data structure so that regardless
         !! of the direction tridiagonal systems are solved efficiently
         !! and fast.
         import :: base_backend_t
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: u_
         class(field_t), intent(in) :: u
         integer, intent(in) :: direction
      end subroutine reorder
   end interface

   abstract interface
      subroutine sum_intox(self, u, u_)
         !! sum9into3 subroutine combines all the directional velocity
         !! derivatives into the corresponding x directional fields.
         import :: base_backend_t
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: u
         class(field_t), intent(in) :: u_
      end subroutine sum_intox
   end interface

   abstract interface
      subroutine vecadd(self, a, x, b, y)
         !! adds two vectors together: y = a*x + b*y
         import :: base_backend_t
         import :: dp
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         real(dp), intent(in) :: a
         class(field_t), intent(in) :: x
         real(dp), intent(in) :: b
         class(field_t), intent(inout) :: y
      end subroutine vecadd
   end interface

   abstract interface
      subroutine get_fields(self, u_out, v_out, w_out, u, v, w)
         !! copy the specialist data structure from device or host back
         !! to a regular 3D data structure.
         import :: base_backend_t
         import :: dp
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         real(dp), dimension(:, :, :), intent(out) :: u_out, v_out, w_out
         class(field_t), intent(in) :: u, v, w
      end subroutine get_fields

      subroutine set_fields(self, u, v, w, u_in, v_in, w_in)
         !! copy the initial condition stored in a regular 3D data
         !! structure into the specialist data structure arrays on the
         !! device or host.
         import :: base_backend_t
         import :: dp
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: u, v, w
         real(dp), dimension(:, :, :), intent(in) :: u_in, v_in, w_in
      end subroutine set_fields
   end interface

   abstract interface
      subroutine alloc_tdsops(self, tdsops, n, dx, operation, scheme, n_halo, &
                              from_to, bc_start, bc_end, sym, c_nu, nu0_nu)
         import :: base_backend_t
         import :: dp
         import :: tdsops_t
         implicit none

         class(base_backend_t) :: self
         class(tdsops_t), allocatable, intent(inout) :: tdsops
         integer, intent(in) :: n
         real(dp), intent(in) :: dx
         character(*), intent(in) :: operation, scheme
         integer, optional, intent(in) :: n_halo
         character(*), optional, intent(in) :: from_to, bc_start, bc_end
         logical, optional, intent(in) :: sym
         real(dp), optional, intent(in) :: c_nu, nu0_nu
      end subroutine alloc_tdsops
   end interface

end module m_base_backend
