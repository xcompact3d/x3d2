module m_omp_backend
   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t
   use m_tdsops, only: dirps_t, tdsops_t
   use m_omp_exec_dist, only: exec_dist_tds_compact, exec_dist_transeq_compact
   use m_omp_sendrecv, only: sendrecv_fields

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
      procedure :: reorder => reorder_omp
      procedure :: sum_yintox => sum_yintox_omp
      procedure :: sum_zintox => sum_zintox_omp
      procedure :: vecadd => vecadd_omp
      procedure :: set_fields => set_fields_omp
      procedure :: get_fields => get_fields_omp
      procedure :: transeq_omp_dist
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

      call self%transeq_omp_dist(du, dv, dw, u, v, w, dirps)

   end subroutine transeq_x_omp

   subroutine transeq_y_omp(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass v, u, w
      call self%transeq_omp_dist(dv, du, dw, v, u, w, dirps)

   end subroutine transeq_y_omp

   subroutine transeq_z_omp(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass w, u, v
      call self%transeq_omp_dist(dw, du, dv, w, u, v, dirps)

   end subroutine transeq_z_omp

   subroutine transeq_omp_dist(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      call transeq_halo_exchange(self, u, v, w, dirps)

      call transeq_dist_component(self, du, u, u, &
            dirps%der1st, dirps%der1st_sym, dirps%der2nd, dirps)
      call transeq_dist_component(self, dv, v, u, &
            dirps%der1st_sym, dirps%der1st, dirps%der2nd_sym, dirps)
      call transeq_dist_component(self, dw, w, u, &
            dirps%der1st_sym, dirps%der1st, dirps%der2nd_sym, dirps)

   end subroutine transeq_omp_dist


   subroutine transeq_halo_exchange(self, u, v, w, dirps)
      class(omp_backend_t) :: self
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps
      integer :: n_halo

      ! TODO: don't hardcode n_halo
      n_halo = 4

      call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, dirps%n, dirps%n_blocks)
      call copy_into_buffers(self%v_send_s, self%v_send_e, v%data, dirps%n, dirps%n_blocks)
      call copy_into_buffers(self%w_send_s, self%w_send_e, w%data, dirps%n, dirps%n_blocks)

      call sendrecv_fields(self%u_recv_s, self%u_recv_e, self%u_send_s, self%u_send_e, &
                           SZ*n_halo*dirps%n_blocks, dirps%nproc, dirps%pprev, dirps%pnext)
      call sendrecv_fields(self%v_recv_s, self%v_recv_e, self%v_send_s, self%v_send_e, &
                           SZ*n_halo*dirps%n_blocks, dirps%nproc, dirps%pprev, dirps%pnext)
      call sendrecv_fields(self%w_recv_s, self%w_recv_e, self%w_send_s, self%w_send_e, &
                           SZ*n_halo*dirps%n_blocks, dirps%nproc, dirps%pprev, dirps%pnext)

   end subroutine transeq_halo_exchange

   !> Computes RHS_x^v following:
   ! rhs_x^v = -0.5*(u*dv/dx + duv/dx) + nu*d2v/dx2
   subroutine transeq_dist_component(self, rhs, v, u, tdsops_du, tdsops_dud, tdsops_d2u, dirps)

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: rhs
      class(field_t), intent(in) :: u, v
      class(tdsops_t), intent(in) :: tdsops_du
      class(tdsops_t), intent(in) :: tdsops_dud
      class(tdsops_t), intent(in) :: tdsops_d2u
      type(dirps_t), intent(in) :: dirps
      class(field_t), pointer :: du, d2u, dud

      du => self%allocator%get_block()
      dud => self%allocator%get_block()
      d2u => self%allocator%get_block()

      call exec_dist_transeq_compact(&
         rhs%data, du%data, dud%data, d2u%data, &
         self%du_send_s,  self%du_send_e,  self%du_recv_s,  self%du_recv_e, &
         self%dud_send_s, self%dud_send_e, self%dud_recv_s, self%dud_recv_e, &
         self%d2u_send_s, self%d2u_send_e, self%d2u_recv_s, self%d2u_recv_e, &
         u%data, self%u_recv_s, self%u_recv_e, &
         v%data, self%v_recv_s, self%v_recv_e, &
         tdsops_du, tdsops_dud, tdsops_d2u, self%nu, &
         dirps%nproc, dirps%pprev, dirps%pnext, dirps%n_blocks)

      call self%allocator%release_block(du)
      call self%allocator%release_block(dud)
      call self%allocator%release_block(d2u)

   end subroutine transeq_dist_component

   subroutine tds_solve_omp(self, du, u, dirps, tdsops)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops

      call tds_solve_dist(self, du, u, dirps, tdsops)

   end subroutine tds_solve_omp

   subroutine tds_solve_dist(self, du, u, dirps, tdsops)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops
      integer :: n_halo

      ! TODO: don't hardcode n_halo
      n_halo = 4
      call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, dirps%n, dirps%n_blocks)

      ! halo exchange
      call sendrecv_fields(self%u_recv_s, self%u_recv_e, self%u_send_s, self%u_send_e, &
                           SZ*n_halo*dirps%n_blocks, dirps%nproc, dirps%pprev, dirps%pnext)


      call exec_dist_tds_compact( &
         du%data, u%data, self%u_recv_s, self%u_recv_e, self%du_send_s, self%du_send_e, &
         self%du_recv_s, self%du_recv_e, &
         tdsops, dirps%nproc, dirps%pprev, dirps%pnext, dirps%n_blocks)

   end subroutine tds_solve_dist

   subroutine reorder_omp(self, u_, u, direction)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u
      integer, intent(in) :: direction

      do k=1,
         do j=1, 
            do i=1, 
               call get_index_reordering(out_i, out_j, out_k, i, j, k, direction)
               u_(out_i, out_j, out_k) = u(i,j,k)
            end do
         end do
      end do

   end subroutine reorder_omp

   subroutine sum_yintox_omp(self, u, u_)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u
      class(field_t), intent(in) :: u_

   end subroutine sum_yintox_omp

   subroutine sum_zintox_omp(self, u, u_)
      implicit none

      class(omp_backend_t) :: self
      class(field_t), intent(inout) :: u
      class(field_t), intent(in) :: u_

   end subroutine sum_zintox_omp

   subroutine vecadd_omp(self, a, x, b, y)
      implicit none

      class(omp_backend_t) :: self
      real(dp), intent(in) :: a
      class(field_t), intent(in) :: x
      real(dp), intent(in) :: b
      class(field_t), intent(inout) :: y

   end subroutine vecadd_omp

   subroutine copy_into_buffers(u_send_s, u_send_e, u, n, n_blocks)
      implicit none

      real(dp), dimension(:, :, :), intent(out) :: u_send_s, u_send_e
      real(dp), dimension(:, :, :), intent(in) :: u
      integer, intent(in) :: n
      integer, intent(in) :: n_blocks
      integer :: i, j, k
      integer :: n_halo = 4

      !$omp parallel do
      do k=1, n_blocks
         do j=1, n_halo
            !$omp simd
            do i=1, SZ
               u_send_s(i, j, k) = u(i, j, k)
               u_send_e(i, j, k) = u(i, n - n_halo + j, k)
            end do
            !$omp end simd
         end do
      end do
      !$omp end parallel do

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

