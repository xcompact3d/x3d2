module m_cuda_backend
   use iso_fortran_env, only: stderr => error_unit
   use cudafor

   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t, &
                       RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y
   use m_poisson_fft, only: poisson_fft_t
   use m_tdsops, only: dirps_t, tdsops_t

   use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
   use m_cuda_common, only: SZ
   use m_cuda_exec_dist, only: exec_dist_transeq_3fused, exec_dist_tds_compact
   use m_cuda_poisson_fft, only: cuda_poisson_fft_t
   use m_cuda_sendrecv, only: sendrecv_fields, sendrecv_3fields
   use m_cuda_tdsops, only: cuda_tdsops_t
   use m_cuda_kernels_dist, only: transeq_3fused_dist, transeq_3fused_subs
   use m_cuda_kernels_reorder, only: &
       reorder_x2y, reorder_x2z, reorder_y2x, reorder_y2z, reorder_z2x, &
       reorder_z2y, sum_yintox, sum_zintox, axpby, buffer_copy

   implicit none

   type, extends(base_backend_t) :: cuda_backend_t
      !character(len=*), parameter :: name = 'cuda'
      integer :: MPI_FP_PREC = dp
      real(dp), device, allocatable, dimension(:, :, :) :: &
         u_recv_s_dev, u_recv_e_dev, u_send_s_dev, u_send_e_dev, &
         v_recv_s_dev, v_recv_e_dev, v_send_s_dev, v_send_e_dev, &
         w_recv_s_dev, w_recv_e_dev, w_send_s_dev, w_send_e_dev, &
         du_send_s_dev, du_send_e_dev, du_recv_s_dev, du_recv_e_dev, &
         dud_send_s_dev, dud_send_e_dev, dud_recv_s_dev, dud_recv_e_dev, &
         d2u_send_s_dev, d2u_send_e_dev, d2u_recv_s_dev, d2u_recv_e_dev
      type(dim3) :: xblocks, xthreads, yblocks, ythreads, zblocks, zthreads
    contains
      procedure :: alloc_tdsops => alloc_cuda_tdsops
      procedure :: transeq_x => transeq_x_cuda
      procedure :: transeq_y => transeq_y_cuda
      procedure :: transeq_z => transeq_z_cuda
      procedure :: tds_solve => tds_solve_cuda
      procedure :: reorder => reorder_cuda
      procedure :: sum_yintox => sum_yintox_cuda
      procedure :: sum_zintox => sum_zintox_cuda
      procedure :: vecadd => vecadd_cuda
      procedure :: set_fields => set_fields_cuda
      procedure :: get_fields => get_fields_cuda
      procedure :: transeq_cuda_dist
      procedure :: transeq_cuda_thom
      procedure :: tds_solve_dist
   end type cuda_backend_t

   interface cuda_backend_t
      module procedure init
   end interface cuda_backend_t

 contains

   function init(globs, allocator, xdirps, ydirps, zdirps) result(backend)
      implicit none

      class(globs_t) :: globs
      class(allocator_t), target, intent(inout) :: allocator
      class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
      type(cuda_backend_t) :: backend

      type(cuda_poisson_fft_t) :: cuda_poisson_fft
      integer :: n_halo, n_block

      select type(allocator)
      type is (cuda_allocator_t)
         ! class level access to the allocator
         backend%allocator => allocator
      end select

      backend%xthreads = dim3(SZ, 1, 1)
      backend%xblocks = dim3(globs%n_groups_x, 1, 1)
      backend%ythreads = dim3(SZ, 1, 1)
      backend%yblocks = dim3(globs%n_groups_y, 1, 1)
      backend%zthreads = dim3(SZ, 1, 1)
      backend%zblocks = dim3(globs%n_groups_z, 1, 1)

      backend%nx_loc = globs%nx_loc
      backend%ny_loc = globs%ny_loc
      backend%nz_loc = globs%nz_loc

      n_halo = 4
      n_block = globs%n_groups_x

      allocate(backend%u_send_s_dev(SZ, n_halo, n_block))
      allocate(backend%u_send_e_dev(SZ, n_halo, n_block))
      allocate(backend%u_recv_s_dev(SZ, n_halo, n_block))
      allocate(backend%u_recv_e_dev(SZ, n_halo, n_block))
      allocate(backend%v_send_s_dev(SZ, n_halo, n_block))
      allocate(backend%v_send_e_dev(SZ, n_halo, n_block))
      allocate(backend%v_recv_s_dev(SZ, n_halo, n_block))
      allocate(backend%v_recv_e_dev(SZ, n_halo, n_block))
      allocate(backend%w_send_s_dev(SZ, n_halo, n_block))
      allocate(backend%w_send_e_dev(SZ, n_halo, n_block))
      allocate(backend%w_recv_s_dev(SZ, n_halo, n_block))
      allocate(backend%w_recv_e_dev(SZ, n_halo, n_block))

      allocate(backend%du_send_s_dev(SZ, 1, n_block))
      allocate(backend%du_send_e_dev(SZ, 1, n_block))
      allocate(backend%du_recv_s_dev(SZ, 1, n_block))
      allocate(backend%du_recv_e_dev(SZ, 1, n_block))
      allocate(backend%dud_send_s_dev(SZ, 1, n_block))
      allocate(backend%dud_send_e_dev(SZ, 1, n_block))
      allocate(backend%dud_recv_s_dev(SZ, 1, n_block))
      allocate(backend%dud_recv_e_dev(SZ, 1, n_block))
      allocate(backend%d2u_send_s_dev(SZ, 1, n_block))
      allocate(backend%d2u_send_e_dev(SZ, 1, n_block))
      allocate(backend%d2u_recv_s_dev(SZ, 1, n_block))
      allocate(backend%d2u_recv_e_dev(SZ, 1, n_block))

      if (globs%use_fft) then
         allocate(cuda_poisson_fft_t :: backend%poisson_fft)

         select type (poisson_fft => backend%poisson_fft)
         type is (cuda_poisson_fft_t)
            poisson_fft = cuda_poisson_fft_t(xdirps, ydirps, zdirps)
         end select
      end if

   end function init

   subroutine alloc_cuda_tdsops( &
      self, tdsops, n, dx, operation, scheme, &
      n_halo, from_to, bc_start, bc_end, sym, c_nu, nu0_nu &
      )
      implicit none

      class(cuda_backend_t) :: self
      class(tdsops_t), allocatable, intent(inout) :: tdsops
      integer, intent(in) :: n
      real(dp), intent(in) :: dx
      character(*), intent(in) :: operation, scheme
      integer, optional, intent(in) :: n_halo
      character(*), optional, intent(in) :: from_to, bc_start, bc_end
      logical, optional, intent(in) :: sym
      real(dp), optional, intent(in) :: c_nu, nu0_nu

      allocate(cuda_tdsops_t :: tdsops)

      select type (tdsops)
      type is (cuda_tdsops_t)
         tdsops = cuda_tdsops_t(n, dx, operation, scheme, n_halo, from_to, &
                                bc_start, bc_end, sym, c_nu, nu0_nu)
      end select

   end subroutine alloc_cuda_tdsops

   subroutine transeq_x_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      call self%transeq_cuda_dist(du, dv, dw, u, v, w, dirps, &
                                  self%xblocks, self%xthreads)

   end subroutine transeq_x_cuda

   subroutine transeq_y_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass v, u, w
      call self%transeq_cuda_dist(dv, du, dw, v, u, w, dirps, &
                                  self%yblocks, self%ythreads)

   end subroutine transeq_y_cuda

   subroutine transeq_z_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass w, u, v
      call self%transeq_cuda_dist(dw, du, dv, w, u, v, dirps, &
                                  self%zblocks, self%zthreads)

   end subroutine transeq_z_cuda

   subroutine transeq_cuda_dist(self, du, dv, dw, u, v, w, dirps, &
                                blocks, threads)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps
      type(dim3), intent(in) :: blocks, threads

      class(field_t), pointer :: temp_du, temp_duu, temp_d2u, &
                                 temp_dv, temp_dvu, temp_d2v, &
                                 temp_dw, temp_dwu, temp_d2w

      real(dp), device, pointer, dimension(:, :, :) :: &
         du_dev, duu_dev, d2u_dev, &
         dv_dev, dvu_dev, d2v_dev, &
         dw_dev, dwu_dev, d2w_dev

      real(dp), device, pointer, dimension(:, :, :) :: u_dev, v_dev, w_dev, &
                                                       ru_dev, rv_dev, rw_dev

      type(cuda_tdsops_t), pointer :: der1st, der1st_sym, der2nd, der2nd_sym

      select type(u); type is (cuda_field_t); u_dev => u%data_d; end select
      select type(v); type is (cuda_field_t); v_dev => v%data_d; end select
      select type(w); type is (cuda_field_t); w_dev => w%data_d; end select

      select type(du); type is (cuda_field_t); ru_dev => du%data_d; end select
      select type(dv); type is (cuda_field_t); rv_dev => dv%data_d; end select
      select type(dw); type is (cuda_field_t); rw_dev => dw%data_d; end select

      select type (tdsops => dirps%der1st)
      type is (cuda_tdsops_t); der1st => tdsops
      end select
      select type (tdsops => dirps%der1st_sym)
      type is (cuda_tdsops_t); der1st_sym => tdsops
      end select
      select type (tdsops => dirps%der2nd)
      type is (cuda_tdsops_t); der2nd => tdsops
      end select
      select type (tdsops => dirps%der2nd_sym)
      type is (cuda_tdsops_t); der2nd_sym => tdsops
      end select

      ! Copy halo data into buffer arrays
      call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, &
                             dirps%n)
      call copy_into_buffers(self%v_send_s_dev, self%v_send_e_dev, v_dev, &
                             dirps%n)
      call copy_into_buffers(self%w_send_s_dev, self%w_send_e_dev, w_dev, &
                             dirps%n)

      ! halo exchange
      call sendrecv_3fields( &
         self%u_recv_s_dev, self%u_recv_e_dev, &
         self%v_recv_s_dev, self%v_recv_e_dev, &
         self%w_recv_s_dev, self%w_recv_e_dev, &
         self%u_send_s_dev, self%u_send_e_dev, &
         self%v_send_s_dev, self%v_send_e_dev, &
         self%w_send_s_dev, self%w_send_e_dev, &
         SZ*4*blocks%x, dirps%nproc, dirps%pprev, dirps%pnext &
      )

      ! get some fields for storing the result
      temp_du => self%allocator%get_block()
      temp_duu => self%allocator%get_block()
      temp_d2u => self%allocator%get_block()

      select type(temp_du)
      type is (cuda_field_t); du_dev => temp_du%data_d
      end select
      select type(temp_duu)
      type is (cuda_field_t); duu_dev => temp_duu%data_d
      end select
      select type(temp_d2u)
      type is (cuda_field_t); d2u_dev => temp_d2u%data_d
      end select

      call exec_dist_transeq_3fused( &
         ru_dev, &
         u_dev, self%u_recv_s_dev, self%u_recv_e_dev, &
         u_dev, self%u_recv_s_dev, self%u_recv_e_dev, &
         du_dev, duu_dev, d2u_dev, &
         self%du_send_s_dev, self%du_send_e_dev, &
         self%du_recv_s_dev, self%du_recv_e_dev, &
         self%dud_send_s_dev, self%dud_send_e_dev, &
         self%dud_recv_s_dev, self%dud_recv_e_dev, &
         self%d2u_send_s_dev, self%d2u_send_e_dev, &
         self%d2u_recv_s_dev, self%d2u_recv_e_dev, &
         der1st, der2nd, self%nu, &
         dirps%nproc, dirps%pprev, dirps%pnext, &
         blocks, threads &
      )

      ! Release temporary blocks
      call self%allocator%release_block(temp_du)
      call self%allocator%release_block(temp_duu)
      call self%allocator%release_block(temp_d2u)

      temp_dv => self%allocator%get_block()
      temp_dvu => self%allocator%get_block()
      temp_d2v => self%allocator%get_block()

      select type(temp_dv)
      type is (cuda_field_t); dv_dev => temp_dv%data_d
      end select
      select type(temp_dvu)
      type is (cuda_field_t); dvu_dev => temp_dvu%data_d
      end select
      select type(temp_d2v)
      type is (cuda_field_t); d2v_dev => temp_d2v%data_d
      end select

      call exec_dist_transeq_3fused( &
         rv_dev, &
         v_dev, self%v_recv_s_dev, self%v_recv_e_dev, &
         u_dev, self%u_recv_s_dev, self%u_recv_e_dev, &
         dv_dev, dvu_dev, d2v_dev, &
         self%du_send_s_dev, self%du_send_e_dev, &
         self%du_recv_s_dev, self%du_recv_e_dev, &
         self%dud_send_s_dev, self%dud_send_e_dev, &
         self%dud_recv_s_dev, self%dud_recv_e_dev, &
         self%d2u_send_s_dev, self%d2u_send_e_dev, &
         self%d2u_recv_s_dev, self%d2u_recv_e_dev, &
         der1st_sym, der2nd_sym, self%nu, &
         dirps%nproc, dirps%pprev, dirps%pnext, &
         blocks, threads &
      )

      ! Release temporary blocks
      call self%allocator%release_block(temp_dv)
      call self%allocator%release_block(temp_dvu)
      call self%allocator%release_block(temp_d2v)

      temp_dw => self%allocator%get_block()
      temp_dwu => self%allocator%get_block()
      temp_d2w => self%allocator%get_block()

      select type(temp_dw)
      type is (cuda_field_t); dw_dev => temp_dw%data_d
      end select
      select type(temp_dwu)
      type is (cuda_field_t); dwu_dev => temp_dwu%data_d
      end select
      select type(temp_d2w)
      type is (cuda_field_t); d2w_dev => temp_d2w%data_d
      end select

      call exec_dist_transeq_3fused( &
         rw_dev, &
         w_dev, self%w_recv_s_dev, self%w_recv_e_dev, &
         u_dev, self%u_recv_s_dev, self%u_recv_e_dev, &
         dw_dev, dwu_dev, d2w_dev, &
         self%du_send_s_dev, self%du_send_e_dev, &
         self%du_recv_s_dev, self%du_recv_e_dev, &
         self%dud_send_s_dev, self%dud_send_e_dev, &
         self%dud_recv_s_dev, self%dud_recv_e_dev, &
         self%d2u_send_s_dev, self%d2u_send_e_dev, &
         self%d2u_recv_s_dev, self%d2u_recv_e_dev, &
         der1st_sym, der2nd_sym, self%nu, &
         dirps%nproc, dirps%pprev, dirps%pnext, &
         blocks, threads &
      )

      ! Release temporary blocks
      call self%allocator%release_block(temp_dw)
      call self%allocator%release_block(temp_dwu)
      call self%allocator%release_block(temp_d2w)

   end subroutine transeq_cuda_dist

   subroutine transeq_cuda_thom(self, du, dv, dw, u, v, w, dirps)
      !! Thomas algorithm implementation. So much more easier than the
      !! distributed algorithm. It is intended to work only on a single rank
      !! so there is no MPI communication.
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

   end subroutine transeq_cuda_thom

   subroutine tds_solve_cuda(self, du, u, dirps, tdsops)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops

      type(dim3) :: blocks, threads

      blocks = dim3(dirps%n_blocks, 1, 1); threads = dim3(SZ, 1, 1)

      call tds_solve_dist(self, du, u, dirps, tdsops, blocks, threads)

   end subroutine tds_solve_cuda

   subroutine tds_solve_dist(self, du, u, dirps, tdsops, blocks, threads)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops
      type(dim3), intent(in) :: blocks, threads

      real(dp), device, pointer, dimension(:, :, :) :: du_dev, u_dev

      type(cuda_tdsops_t), pointer :: tdsops_dev

      select type(du); type is (cuda_field_t); du_dev => du%data_d; end select
      select type(u); type is (cuda_field_t); u_dev => u%data_d; end select

      select type (tdsops)
      type is (cuda_tdsops_t); tdsops_dev => tdsops
      end select

      call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, &
                             tdsops_dev%n)

      call sendrecv_fields(self%u_recv_s_dev, self%u_recv_e_dev, &
                           self%u_send_s_dev, self%u_send_e_dev, &
                           SZ*4*blocks%x, dirps%nproc, &
                           dirps%pprev, dirps%pnext)

      ! call exec_dist
      call exec_dist_tds_compact( &
         du_dev, u_dev, &
         self%u_recv_s_dev, self%u_recv_e_dev, &
         self%du_send_s_dev, self%du_send_e_dev, &
         self%du_recv_s_dev, self%du_recv_e_dev, &
         tdsops_dev, dirps%nproc, dirps%pprev, dirps%pnext, &
         blocks, threads &
      )

   end subroutine tds_solve_dist

   subroutine reorder_cuda(self, u_o, u_i, direction)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u_o
      class(field_t), intent(in) :: u_i
      integer, intent(in) :: direction

      real(dp), device, pointer, dimension(:, :, :) :: u_o_d, u_i_d
      type(dim3) :: blocks, threads

      select type(u_o); type is (cuda_field_t); u_o_d => u_o%data_d; end select
      select type(u_i); type is (cuda_field_t); u_i_d => u_i%data_d; end select

      select case (direction)
      case (RDR_X2Y) ! x2y
         blocks = dim3(self%nx_loc/SZ, self%nz_loc, self%ny_loc/SZ)
         threads = dim3(SZ, SZ, 1)
         call reorder_x2y<<<blocks, threads>>>(u_o_d, u_i_d, self%nz_loc)
      case (RDR_X2Z) ! x2z
         blocks = dim3(self%nx_loc, self%ny_loc/SZ, 1)
         threads = dim3(SZ, 1, 1)
         call reorder_x2z<<<blocks, threads>>>(u_o_d, u_i_d, self%nz_loc)
      case (RDR_Y2X) ! y2x
         blocks = dim3(self%nx_loc/SZ, self%ny_loc/SZ, self%nz_loc)
         threads = dim3(SZ, SZ, 1)
         call reorder_y2x<<<blocks, threads>>>(u_o_d, u_i_d, self%nz_loc)
      case (RDR_Y2Z) ! y2z
         blocks = dim3(self%nx_loc/SZ, self%ny_loc/SZ, self%nz_loc)
         threads = dim3(SZ, SZ, 1)
         call reorder_y2z<<<blocks, threads>>>(u_o_d, u_i_d, &
                                               self%nx_loc, self%nz_loc)
      case (RDR_Z2X) ! z2x
         blocks = dim3(self%nx_loc, self%ny_loc/SZ, 1)
         threads = dim3(SZ, 1, 1)
         call reorder_z2x<<<blocks, threads>>>(u_o_d, u_i_d, self%nz_loc)
      case (RDR_Z2Y) ! z2y
         blocks = dim3(self%nx_loc/SZ, self%ny_loc/SZ, self%nz_loc)
         threads = dim3(SZ, SZ, 1)
         call reorder_z2y<<<blocks, threads>>>(u_o_d, u_i_d, &
                                               self%nx_loc, self%nz_loc)
      case default
         error stop 'Reorder direction is undefined.'
      end select

   end subroutine reorder_cuda

   subroutine sum_yintox_cuda(self, u, u_y)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u
      class(field_t), intent(in) :: u_y

      real(dp), device, pointer, dimension(:, :, :) :: u_d, u_y_d
      type(dim3) :: blocks, threads

      select type(u); type is (cuda_field_t); u_d => u%data_d; end select
      select type(u_y); type is (cuda_field_t); u_y_d => u_y%data_d; end select

      blocks = dim3(self%nx_loc/SZ, self%ny_loc/SZ, self%nz_loc)
      threads = dim3(SZ, SZ, 1)
      call sum_yintox<<<blocks, threads>>>(u_d, u_y_d, self%nz_loc)

   end subroutine sum_yintox_cuda

   subroutine sum_zintox_cuda(self, u, u_z)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u
      class(field_t), intent(in) :: u_z

      real(dp), device, pointer, dimension(:, :, :) :: u_d, u_z_d
      type(dim3) :: blocks, threads

      select type(u); type is (cuda_field_t); u_d => u%data_d; end select
      select type(u_z); type is (cuda_field_t); u_z_d => u_z%data_d; end select

      blocks = dim3(self%nx_loc, self%ny_loc/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call sum_zintox<<<blocks, threads>>>(u_d, u_z_d, self%nz_loc)

   end subroutine sum_zintox_cuda

   subroutine vecadd_cuda(self, a, x, b, y)
      implicit none

      class(cuda_backend_t) :: self
      real(dp), intent(in) :: a
      class(field_t), intent(in) :: x
      real(dp), intent(in) :: b
      class(field_t), intent(inout) :: y

      real(dp), device, pointer, dimension(:, :, :) :: x_d, y_d
      type(dim3) :: blocks, threads
      integer :: nx

      select type(x); type is (cuda_field_t); x_d => x%data_d; end select
      select type(y); type is (cuda_field_t); y_d => y%data_d; end select

      nx = size(x_d, dim = 2)
      blocks = dim3(size(x_d, dim = 3), 1, 1)
      threads = dim3(SZ, 1, 1)
      call axpby<<<blocks, threads>>>(nx, a, x_d, b, y_d)

   end subroutine vecadd_cuda

   subroutine copy_into_buffers(u_send_s_dev, u_send_e_dev, u_dev, n)
      implicit none

      real(dp), device, dimension(:, :, :), intent(out) :: u_send_s_dev, &
                                                           u_send_e_dev
      real(dp), device, dimension(:, :, :), intent(in) :: u_dev
      integer, intent(in) :: n

      type(dim3) :: blocks, threads
      integer :: n_halo = 4

      blocks = dim3(size(u_dev, dim = 3), 1, 1)
      threads = dim3(SZ, 1, 1)
      call buffer_copy<<<blocks, threads>>>(u_send_s_dev, u_send_e_dev, &
                                            u_dev, n, n_halo)

   end subroutine copy_into_buffers

   subroutine set_fields_cuda(self, u, v, w, u_in, v_in, w_in)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u, v, w
      real(dp), dimension(:, :, :), intent(in) :: u_in, v_in, w_in

      select type(u); type is (cuda_field_t); u%data_d = u_in; end select
      select type(v); type is (cuda_field_t); v%data_d = v_in; end select
      select type(w); type is (cuda_field_t); w%data_d = w_in; end select

   end subroutine set_fields_cuda

   subroutine get_fields_cuda(self, u_out, v_out, w_out, u, v, w)
      implicit none

      class(cuda_backend_t) :: self
      real(dp), dimension(:, :, :), intent(out) :: u_out, v_out, w_out
      class(field_t), intent(in) :: u, v, w

      select type(u); type is (cuda_field_t); u_out = u%data_d; end select
      select type(v); type is (cuda_field_t); v_out = v%data_d; end select
      select type(w); type is (cuda_field_t); w_out = w%data_d; end select

   end subroutine get_fields_cuda

end module m_cuda_backend

