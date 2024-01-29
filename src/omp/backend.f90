module m_omp_backend
   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t
   use m_tdsops, only: dirps_t, tdsops_t

   use m_omp_common, only: SZ

   implicit none

   type, extends(base_backend_t) :: omp_backend_t
      !character(len=*), parameter :: name = 'omp'
      integer :: MPI_FP_PREC = dp
      real(dp), allocatable, dimension(:, :, :) :: &
         u_recv_s, u_recv_e, u_send_s, u_send_e, &
         v_recv_s, v_recv_e, v_send_s, v_send_e, &
         w_recv_s, w_recv_e, w_send_s, w_send_e, &
         du_send_s, du_send_e, du_recv_s, du_recv_e, &
         dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
         d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e
    contains
      procedure :: alloc_tdsops => alloc_omp_tdsops
      procedure :: transeq_x => transeq_x_omp
      procedure :: transeq_y => transeq_y_omp
      procedure :: transeq_z => transeq_z_omp
      procedure :: tds_solve => tds_solve_omp
      procedure :: trans_x2y => trans_x2y_omp
      procedure :: trans_x2z => trans_x2z_omp
      procedure :: trans_y2z => trans_y2z_omp
      procedure :: trans_z2y => trans_z2y_omp
      procedure :: trans_y2x => trans_y2x_omp
      procedure :: sum_yzintox => sum_yzintox_omp
      procedure :: vecadd => vecadd_omp
      procedure :: set_fields => set_fields_omp
      procedure :: get_fields => get_fields_omp
   end type omp_backend_t

   interface omp_backend_t
      module procedure init
   end interface omp_backend_t

 contains

   function init(globs, allocator) result(backend)
      implicit none

      class(globs_t) :: globs
      class(allocator_t), target, intent(inout) :: allocator
      type(omp_backend_t) :: backend

      integer :: n_halo, n_block

      select type(allocator)
      type is (allocator_t)
         ! class level access to the allocator
         backend%allocator => allocator
      end select

      n_halo = 4
      n_block = globs%n_groups_x

      allocate(backend%u_send_s(SZ, n_halo, n_block))
      allocate(backend%u_send_e(SZ, n_halo, n_block))
      allocate(backend%u_recv_s(SZ, n_halo, n_block))
      allocate(backend%u_recv_e(SZ, n_halo, n_block))
      allocate(backend%v_send_s(SZ, n_halo, n_block))
      allocate(backend%v_send_e(SZ, n_halo, n_block))
      allocate(backend%v_recv_s(SZ, n_halo, n_block))
      allocate(backend%v_recv_e(SZ, n_halo, n_block))
      allocate(backend%w_send_s(SZ, n_halo, n_block))
      allocate(backend%w_send_e(SZ, n_halo, n_block))
      allocate(backend%w_recv_s(SZ, n_halo, n_block))
      allocate(backend%w_recv_e(SZ, n_halo, n_block))

      allocate(backend%du_send_s(SZ, 1, n_block))
      allocate(backend%du_send_e(SZ, 1, n_block))
      allocate(backend%du_recv_s(SZ, 1, n_block))
      allocate(backend%du_recv_e(SZ, 1, n_block))
      allocate(backend%dud_send_s(SZ, 1, n_block))
      allocate(backend%dud_send_e(SZ, 1, n_block))
      allocate(backend%dud_recv_s(SZ, 1, n_block))
      allocate(backend%dud_recv_e(SZ, 1, n_block))
      allocate(backend%d2u_send_s(SZ, 1, n_block))
      allocate(backend%d2u_send_e(SZ, 1, n_block))
      allocate(backend%d2u_recv_s(SZ, 1, n_block))
      allocate(backend%d2u_recv_e(SZ, 1, n_block))

   end function init

   subroutine alloc_omp_tdsops( &
      self, tdsops, n, dx, operation, scheme, &
      n_halo, from_to, bc_start, bc_end, sym, c_nu, nu0_nu &
      )
      implicit none

      class(omp_backend_t) :: self
      class(tdsops_t), allocatable, intent(inout) :: tdsops
      integer, intent(in) :: n
      real(dp), intent(in) :: dx
      character(*), intent(in) :: operation, scheme
      integer, optional, intent(in) :: n_halo
      character(*), optional, intent(in) :: from_to, bc_start, bc_end
      logical, optional, intent(in) :: sym
      real(dp), optional, intent(in) :: c_nu, nu0_nu

      allocate(tdsops_t :: tdsops)

      select type (tdsops)
      type is (tdsops_t)
         tdsops = tdsops_t(n, dx, operation, scheme, n_halo, from_to, &
                           bc_start, bc_end, sym, c_nu, nu0_nu)
      end select

   end subroutine alloc_omp_tdsops

   subroutine transeq_x_omp(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      !call self%transeq_omp_dist(du, dv, dw, u, v, w, dirps)

   end subroutine transeq_x_omp

   subroutine transeq_y_omp(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass v, u, w
      !call self%transeq_omp_dist(dv, du, dw, v, u, w, dirps)

   end subroutine transeq_y_omp

   subroutine transeq_z_omp(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass w, u, v
      !call self%transeq_omp_dist(dw, du, dv, w, u, v, dirps)

   end subroutine transeq_z_omp

   subroutine tds_solve_omp(self, du, u, dirps, tdsops)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops

      !call self%tds_solve_dist(self, du, u, dirps, tdsops)

   end subroutine tds_solve_omp

   subroutine trans_x2y_omp(self, u_, u)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u

   end subroutine trans_x2y_omp

   subroutine trans_x2z_omp(self, u_, u)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u

   end subroutine trans_x2z_omp

   subroutine trans_y2z_omp(self, u_, u)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u

   end subroutine trans_y2z_omp

   subroutine trans_z2y_omp(self, u_, u)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u

   end subroutine trans_z2y_omp

   subroutine trans_y2x_omp(self, u_, u)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u

   end subroutine trans_y2x_omp

   subroutine sum_yzintox_omp(self, du, dv, dw, &
                               du_y, dv_y, dw_y, du_z, dv_z, dw_z)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: du_y, dv_y, dw_y, du_z, dv_z, dw_z

   end subroutine sum_yzintox_omp

   subroutine vecadd_omp(self, a, x, b, y)
      implicit none

      class(omp_backend_t) :: self
      real(dp), intent(in) :: a
      class(field_t), intent(in) :: x
      real(dp), intent(in) :: b
      class(field_t), intent(inout) :: y

   end subroutine vecadd_omp

   subroutine copy_into_buffers(u_send_s, u_send_e, u, n)
      implicit none

      real(dp), dimension(:, :, :), intent(out) :: u_send_s, u_send_e
      real(dp), dimension(:, :, :), intent(in) :: u
      integer, intent(in) :: n

      u_send_s(:, :, :) = u(:, 1:4, :)
      u_send_e(:, :, :) = u(:, n - 3:n, :)

   end subroutine copy_into_buffers

   subroutine set_fields_omp(self, u, v, w, u_in, v_in, w_in)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u, v, w
      real(dp), dimension(:, :, :), intent(in) :: u_in, v_in, w_in

      u%data = u_in
      v%data = v_in
      w%data = w_in

   end subroutine set_fields_omp

   subroutine get_fields_omp(self, u_out, v_out, w_out, u, v, w)
      implicit none

      class(omp_backend_t) :: self
      real(dp), dimension(:, :, :), intent(out) :: u_out, v_out, w_out
      class(field_t), intent(in) :: u, v, w

      u_out = u%data
      v_out = v%data
      w_out = w%data

   end subroutine get_fields_omp

end module m_omp_backend

