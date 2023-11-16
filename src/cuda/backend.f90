module m_cuda_backend
   use cudafor

   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t
   use m_tdsops, only: dirps_t

   use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
   use m_cuda_common, only: SZ
   use m_cuda_exec_dist, only: exec_dist_transeq_3fused
   use m_cuda_sendrecv, only: sendrecv_3fields
   use m_cuda_tdsops, only: cuda_tdsops_t
   use m_cuda_kernels_dist, only: transeq_3fused_dist, transeq_3fused_subs

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
      procedure :: transeq_x => transeq_x_cuda
      procedure :: transeq_y => transeq_y_cuda
      procedure :: transeq_z => transeq_z_cuda
      procedure :: trans_x2y => trans_x2y_cuda
      procedure :: trans_x2z => trans_x2z_cuda
      procedure :: sum_yzintox => sum_yzintox_cuda
      procedure :: transeq_cuda_dist
      procedure :: transeq_cuda_thom
   end type cuda_backend_t

   interface cuda_backend_t
      module procedure constructor
   end interface cuda_backend_t

 contains

   function constructor(globs, allocator, xdirps, ydirps, zdirps) &
      result(backend)
      implicit none

      class(globs_t) :: globs
      class(allocator_t), target, intent(inout) :: allocator
      class(dirps_t), target, intent(inout) :: xdirps, ydirps, zdirps
      type(cuda_backend_t) :: backend

      integer :: n_halo, n_block

      select type(allocator)
      type is (cuda_allocator_t)
         ! class level access to the allocator
         backend%allocator => allocator
      end select

      ! class level access to derivative parameters
      backend%xdirps => xdirps
      backend%ydirps => ydirps
      backend%zdirps => zdirps

      backend%xthreads = dim3(SZ, 1, 1)
      backend%xblocks = dim3(globs%n_groups_x, 1, 1)

      allocate(cuda_tdsops_t :: backend%xdirps%der1st)
      allocate(cuda_tdsops_t :: backend%ydirps%der1st)
      allocate(cuda_tdsops_t :: backend%zdirps%der1st)
      allocate(cuda_tdsops_t :: backend%xdirps%der1st_sym)
      allocate(cuda_tdsops_t :: backend%ydirps%der1st_sym)
      allocate(cuda_tdsops_t :: backend%zdirps%der1st_sym)
      allocate(cuda_tdsops_t :: backend%xdirps%der2nd)
      allocate(cuda_tdsops_t :: backend%ydirps%der2nd)
      allocate(cuda_tdsops_t :: backend%zdirps%der2nd)
      allocate(cuda_tdsops_t :: backend%xdirps%der2nd_sym)
      allocate(cuda_tdsops_t :: backend%ydirps%der2nd_sym)
      allocate(cuda_tdsops_t :: backend%zdirps%der2nd_sym)

      select type (der1st => backend%xdirps%der1st)
      type is (cuda_tdsops_t)
         der1st = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                'first-deriv', 'compact6')
      end select
      select type (der1st => backend%ydirps%der1st)
      type is (cuda_tdsops_t)
         der1st = cuda_tdsops_t(globs%ny_loc, globs%dy, &
                                'first-deriv', 'compact6')
      end select
      select type (der1st => backend%zdirps%der1st)
      type is (cuda_tdsops_t)
         der1st = cuda_tdsops_t(globs%nz_loc, globs%dz, &
                                'first-deriv', 'compact6')
      end select
      select type (der1st_sym => backend%xdirps%der1st_sym)
      type is (cuda_tdsops_t)
         der1st_sym = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                    'first-deriv', 'compact6')
      end select
      select type (der1st_sym => backend%ydirps%der1st_sym)
      type is (cuda_tdsops_t)
         der1st_sym = cuda_tdsops_t(globs%ny_loc, globs%dy, &
                                    'first-deriv', 'compact6')
      end select
      select type (der1st_sym => backend%zdirps%der1st_sym)
      type is (cuda_tdsops_t)
         der1st_sym = cuda_tdsops_t(globs%nz_loc, globs%dz, &
                                    'first-deriv', 'compact6')
      end select
      select type (der2nd => backend%xdirps%der2nd)
      type is (cuda_tdsops_t)
         der2nd = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                'second-deriv', 'compact6')
      end select
      select type (der2nd => backend%ydirps%der2nd)
      type is (cuda_tdsops_t)
         der2nd = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                'second-deriv', 'compact6')
      end select
      select type (der2nd => backend%zdirps%der2nd)
      type is (cuda_tdsops_t)
         der2nd = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                'second-deriv', 'compact6')
      end select
      select type (der2nd_sym => backend%xdirps%der2nd_sym)
      type is (cuda_tdsops_t)
         der2nd_sym = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                    'second-deriv', 'compact6')
      end select
      select type (der2nd_sym => backend%ydirps%der2nd_sym)
      type is (cuda_tdsops_t)
         der2nd_sym = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                    'second-deriv', 'compact6')
      end select
      select type (der2nd_sym => backend%zdirps%der2nd_sym)
      type is (cuda_tdsops_t)
         der2nd_sym = cuda_tdsops_t(globs%nx_loc, globs%dx, &
                                    'second-deriv', 'compact6')
      end select

      n_halo = 4
      n_block = ydirps%n*zdirps%n/SZ

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

      ! Assign transeq_? into right functions
      ! The idea is that these assignments will be conditional
      !backend%transeq_x => transeq_cuda_dist
      !backend%transeq_x => transeq_cuda_thom
      !backend%transeq_y => transeq_cuda_dist
      !backend%transeq_z => transeq_cuda_dist

   end function constructor

   subroutine transeq_x_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      call self%transeq_cuda_dist(du, dv, dw, u, v, w, dirps, &
                                  self%xthreads, self%xblocks)

   end subroutine transeq_x_cuda

   subroutine transeq_y_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass v, u, w
      call self%transeq_cuda_dist(dv, du, dw, v, u, w, dirps, &
                                  self%ythreads, self%yblocks)

   end subroutine transeq_y_cuda

   subroutine transeq_z_cuda(self, du, dv, dw, u, v, w, dirps)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps

      ! u, v, w is reordered so that we pass w, u, v
      call self%transeq_cuda_dist(dw, du, dv, w, u, v, dirps, &
                                  self%zthreads, self%zblocks)

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
      !self%u_send_s_dev(:, :, :) = u_dev(:, 1:4, :)
      !self%u_send_e_dev(:, :, :) = u_dev(:, dirps%n - 3:dirps%n, :)
      !self%v_send_s_dev(:, :, :) = v_dev(:, 1:4, :)
      !self%v_send_e_dev(:, :, :) = v_dev(:, dirps%n - 3:dirps%n, :)
      !self%w_send_s_dev(:, :, :) = w_dev(:, 1:4, :)
      !self%w_send_e_dev(:, :, :) = w_dev(:, dirps%n - 3:dirps%n, :)

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

!      call transeq_fused_thom_pp<<<derps%blocks, derps%threads>>>( &
!         du, u, conv, derps%ff_dev, derps%fs_dev, derps%fw_dev, derps%fp_dev, &
!         derps%sf_dev, derps%ss_dev, derps%sw_dev, derps%sp_dev, &
!         derps%n, self%nu, derps%alfai, derps%alsai, &
!         derps%afi, derps%bfi, derps%asi, derps%bsi, derps%csi, derps%dsi &
!      )
!      call transeq_fused_thom_pp<<<derps%blocks, derps%threads>>>( &
!         dv, v, conv, derps%ff_dev, derps%fs_dev, derps%fw_dev, derps%fp_dev, &
!         derps%sf_dev, derps%ss_dev, derps%sw_dev, derps%sp_dev, &
!         derps%n, self%nu, derps%alfai, derps%alsai, &
!         derps%afi, derps%bfi, derps%asi, derps%bsi, derps%csi, derps%dsi &
!      )
!      call transeq_fused_thom_pp<<<derps%blocks, derps%threads>>>( &
!         dw, w, conv, derps%ff_dev, derps%fs_dev, derps%fw_dev, derps%fp_dev, &
!         derps%sf_dev, derps%ss_dev, derps%sw_dev, derps%sp_dev, &
!         derps%n, self%nu, derps%alfai, derps%alsai, &
!         derps%afi, derps%bfi, derps%asi, derps%bsi, derps%csi, derps%dsi &
!      )

   end subroutine transeq_cuda_thom

   subroutine trans_x2y_cuda(self, u_y, v_y, w_y, u, v, w)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u_y, v_y, w_y
      class(field_t), intent(in) :: u, v, w

   end subroutine trans_x2y_cuda

   subroutine trans_x2z_cuda(self, u_z, v_z, w_z, u, v, w)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: u_z, v_z, w_z
      class(field_t), intent(in) :: u, v, w

   end subroutine trans_x2z_cuda

   subroutine sum_yzintox_cuda(self, du, dv, dw, &
                               du_y, dv_y, dw_y, du_z, dv_z, dw_z)
      implicit none

      class(cuda_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: du_y, dv_y, dw_y, du_z, dv_z, dw_z

   end subroutine sum_yzintox_cuda

end module m_cuda_backend

