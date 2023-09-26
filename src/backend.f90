module m_base_backend
   use m_allocator, only: allocator_t, field_t
   use m_common, only: derps_t

   implicit none

   type, abstract :: base_backend_t
      !! This is the base backend type, or class. Here we define all the
      !! operations we need at a very high level. Then, we have specific
      !! operations that will be implemented on each backend, and these are
      !! defined here using abstract interfaces.
      !!
      !! transeq is a good one to start with. It executes the derivations in x,
      !! y, and z directions using the corresponding procedure pointers. There
      !! are two different algorithms implemented for this operation,
      !! a distributed algorithm and the Thomas algorithm. At this base backend
      !! level we don't know which algorithm will be executed, that is decided
      !! at run time and therefore backend implementations are responsible from
      !! assigning these pointers to the right subroutines.
      !!
      !! trans_?2? subroutines are straightforward, they rearrange data into
      !! our specialist data structure so that regardless of the direction
      !! tridiagonal systems are solved efficiently and fast.
      !!
      !! Allocator is responsible from giving us a field sized array when
      !! requested. For example, when the derivations in x direction are
      !! completed and we are ready for the y directional derivatives, we need
      !! three fields to reorder and store the velocities in y direction. Also,
      !! we need three more fields for storing the results, and the get_block
      !! method of the allocator is used to arrange all these memory
      !! assignments. Later, when a field is no more required, release_block
      !! method of the allocator can be used to make this field available
      !! for later use.
      real :: nu
      class(allocator_t), pointer :: allocator
      class(derps_t), pointer :: xderps, yderps, zderps
   contains
      procedure :: transeq
      procedure :: run
      procedure(transeq_ders), deferred :: transeq_x
      procedure(transeq_ders), deferred :: transeq_y
      procedure(transeq_ders), deferred :: transeq_z
      procedure(transposer), deferred :: trans_x2y
      procedure(transposer), deferred :: trans_x2z
      procedure(sum9into3), deferred :: sum_yzintox
   end type base_backend_t

   abstract interface
      subroutine transeq_ders(self, du, duu, d2u, u, v, w, conv, derps)
         import :: base_backend_t
         import :: field_t
         import :: derps_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(out) :: du, duu, d2u
         class(field_t), intent(in) :: u, v, w, conv
         type(derps_t), intent(in) :: derps
      end subroutine transeq_ders
   end interface

   abstract interface
      subroutine transposer(self, u_, v_, w_, u, v, w)
         import :: base_backend_t
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(out) :: u_, v_, w_
         class(field_t), intent(in) :: u, v, w
      end subroutine transposer
   end interface

   abstract interface
      subroutine sum9into3(self, du, dv, dw, du_y, dv_y, dw_y, du_z, dv_z, dw_z)
         import :: base_backend_t
         import :: field_t
         implicit none

         class(base_backend_t) :: self
         class(field_t), intent(inout) :: du, dv, dw
         class(field_t), intent(in) :: du_y, dv_y, dw_y, du_z, dv_z, dw_z
      end subroutine sum9into3
   end interface

contains

   subroutine transeq(self, du, dv, dw, u, v, w)
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(out) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w

      class(field_t), pointer :: u_y, v_y, w_y, u_z, v_z, w_z, du_y, dv_y, dw_y, du_z, dv_z, dw_z

      !1/2(nabla u curl u + u nabla u) + nu nablasq u

      ! call derivatives in x direction. Based on the run time arguments this
      ! executes a distributed algorithm or the Thomas algorithm.
      call transeq_x(du, dv, dw, u, v, w, u, self%xderps)

      ! request fields from the allocator
      u_y => self%allocator%get_block()
      v_y => self%allocator%get_block()
      w_y => self%allocator%get_block()
      du_y => self%allocator%get_block()
      dv_y => self%allocator%get_block()
      dw_y => self%allocator%get_block()

      ! reorder data from x orientation to y orientation
      call trans_x2y(u_y, v_y, w_y, u, v, w)
      ! similar to the x direction, obtain derivatives in y.
      call transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, v_y, self%yderps)

      ! we don't need the velocities in y orientation any more, so release
      ! them to open up space.
      ! It is important that this doesn't actually deallocate any memory,
      ! it just makes the corresponding memory space available for use.
      call self%allocator%release_block(u_y)
      call self%allocator%release_block(v_y)
      call self%allocator%release_block(w_y)

      ! just like in y direction, get some fields for the z derivatives.
      u_z => self%allocator%get_block()
      v_z => self%allocator%get_block()
      w_z => self%allocator%get_block()
      du_z => self%allocator%get_block()
      dv_z => self%allocator%get_block()
      dw_z => self%allocator%get_block()

      ! reorder from x to z
      call trans_x2z(u_z, v_z, w_z, u, v, w)
      ! get the derivatives in z
      call transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, w_z, self%zderps)

      ! there is no need to keep velocities in z orientation around, so release
      call self%allocator%release_block(u_z)
      call self%allocator%release_block(v_z)
      call self%allocator%release_block(w_z)

      ! gather all the contributions into the x result array
      ! this function does the data reordering and summations at once.
      call sumall(du, dv, dw, du_y, dv_y, dw_y, du_z, dv_z, dw_z)

      ! release all the unnecessary blocks.
      call self%allocator%release_block(du_y)
      call self%allocator%release_block(dv_y)
      call self%allocator%release_block(dw_y)
      call self%allocator%release_block(du_z)
      call self%allocator%release_block(dv_z)
      call self%allocator%release_block(dw_z)

   end subroutine transeq

   subroutine run(self)
      !! This is the main subroutine that will start execution. This may be
      !! moved to the main program but is here for now.
      implicit none

      class(base_backend_t) :: self
      class(field_t), pointer :: u, v, w, du, dv, dw, udiv, p, px, py, pz

      real :: dt, alpha

      ! transport equation
      call self%transeq(du, dv, dw, u, v, w)
      ! time integration

      !! pressure stuff
      !call self%divergence(u, v, w, udiv)
      !call self%poisson(udiv, p)
      !call self%gradient(p, px, py, pz)
      !! velocity correction

   end subroutine run

end module m_base_backend

