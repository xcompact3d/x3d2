module m_base_backend
   use m_allocator, only: allocator_t, field_t
   use m_common, only derps_t

   implicit none

   type, abstract :: base_backend_t
      real :: nu
      class(allocator_t), pointer :: allocator
      class(derps_t), pointer :: xderps, yderps, zderps
   contains
      procedure :: transeq
      procedure :: run
      procedure(transeq_ders), pointer :: transeq_x, transeq_y, transeq_z
      procedure(transposer), deferred :: trans_x2y, trans_x2z
      procedure(sum9into3), deferred :: sum_yzintox
   end type base_backend_t

   abstract interface
      subroutine transeq_ders(self, du, duu, d2u, u, v, w, conv, derps)
         implicit none
         import :: base_backend_t
         import :: field_t

         class(base_backend_t) :: self
         class(field_t), intent(out) :: du, duu, d2u
         class(field_t), intent(in) :: u, v, w, conv
         type(derps_t), intent(in) :: derps
      end subroutine transeq_ders
   end interface

   abstract interface
      subroutine transposer(self, u_, v_, w_, u, v, w)
         implicit none
         import :: base_backend_t
         import :: field_t

         class(base_backend_t) :: self
         class(field_t), intent(out) :: u_, v_, w_
         class(field_t), intent(in) :: u, v, w
      end subroutine transposer
   end interface

   abstract interface
      subroutine sum9into3(self, du, dv, dw, du_y, dv_y, dw_y, du_z, dv_z, dw_z)
         import :: base_backend_t
         import :: field_t

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

      call transeq_x(du, dv, dw, u, v, w, u, self%xderps)

      u_y => self%allocator%get_block()
      v_y => self%allocator%get_block()
      w_y => self%allocator%get_block()
      du_y => self%allocator%get_block()
      dv_y => self%allocator%get_block()
      dw_y => self%allocator%get_block()

      call trans_x2y(u_y, v_y, w_y, u, v, w)
      call transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, v_y, self%yderps)

      call self%allocator%release_block(u_y)
      call self%allocator%release_block(v_y)
      call self%allocator%release_block(w_y)

      u_z => self%allocator%get_block()
      v_z => self%allocator%get_block()
      w_z => self%allocator%get_block()
      du_z => self%allocator%get_block()
      dv_z => self%allocator%get_block()
      dw_z => self%allocator%get_block()

      call trans_x2z(u_z, v_z, w_z, u, v, w)
      call transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, w_z, self%zderps)

      call self%allocator%release_block(u_z)
      call self%allocator%release_block(v_z)
      call self%allocator%release_block(w_z)

      call sumall(du, dv, dw, du_y, dv_y, dw_y, du_z, dv_z, dw_z)

      call self%allocator%release_block(du_y)
      call self%allocator%release_block(dv_y)
      call self%allocator%release_block(dw_y)
      call self%allocator%release_block(du_z)
      call self%allocator%release_block(dv_z)
      call self%allocator%release_block(dw_z)

   end subroutine transeq

   subroutine run(self)
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

