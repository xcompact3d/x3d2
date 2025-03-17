module m_cuda_backend
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, move_data_loc, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      RDR_C2X, RDR_C2Y, RDR_C2Z, RDR_X2C, RDR_Y2C, RDR_Z2C, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, NULL_LOC, &
                      X_FACE, Y_FACE, Z_FACE
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t, tdsops_t

  use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_transeq_3fused, exec_dist_tds_compact
  use m_cuda_poisson_fft, only: cuda_poisson_fft_t
  use m_cuda_sendrecv, only: sendrecv_fields, sendrecv_3fields
  use m_cuda_tdsops, only: cuda_tdsops_t
  use m_cuda_kernels_dist, only: transeq_3fused_dist, transeq_3fused_subs
  use m_cuda_kernels_fieldops, only: axpby, buffer_copy, field_scale, &
                                     field_shift, scalar_product, &
                                     field_max_sum, field_set_y_face, &
                                     volume_integral
  use m_cuda_kernels_reorder, only: reorder_x2y, reorder_x2z, reorder_y2x, &
                                    reorder_y2z, reorder_z2x, reorder_z2y, &
                                    reorder_c2x, reorder_x2c, &
                                    sum_yintox, sum_zintox

  implicit none

  private :: transeq_halo_exchange, transeq_dist_component

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
    procedure :: scalar_product => scalar_product_cuda
    procedure :: field_max_mean => field_max_mean_cuda
    procedure :: field_scale => field_scale_cuda
    procedure :: field_shift => field_shift_cuda
    procedure :: field_set_face => field_set_face_cuda
    procedure :: field_volume_integral => field_volume_integral_cuda
    procedure :: copy_data_to_f => copy_data_to_f_cuda
    procedure :: copy_f_to_data => copy_f_to_data_cuda
    procedure :: init_poisson_fft => init_cuda_poisson_fft
    procedure :: transeq_cuda_dist
    procedure :: transeq_cuda_thom
    procedure :: tds_solve_dist
  end type cuda_backend_t

  interface cuda_backend_t
    module procedure init
  end interface cuda_backend_t

contains

  function init(mesh, allocator) result(backend)
    implicit none

    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator
    type(cuda_backend_t) :: backend

    type(cuda_poisson_fft_t) :: cuda_poisson_fft
    integer :: n_halo, n_groups

    call backend%base_init()

    select type (allocator)
    type is (cuda_allocator_t)
      ! class level access to the allocator
      backend%allocator => allocator
    end select
    backend%mesh => mesh

    backend%xthreads = dim3(SZ, 1, 1)
    backend%xblocks = dim3(backend%mesh%get_n_groups(DIR_X), 1, 1)
    backend%ythreads = dim3(SZ, 1, 1)
    backend%yblocks = dim3(backend%mesh%get_n_groups(DIR_Y), 1, 1)
    backend%zthreads = dim3(SZ, 1, 1)
    backend%zblocks = dim3(backend%mesh%get_n_groups(DIR_Z), 1, 1)

    n_halo = 4
    ! Buffer size should be big enough for the largest MPI exchange.
    n_groups = maxval([backend%mesh%get_n_groups(DIR_X), &
                       backend%mesh%get_n_groups(DIR_Y), &
                       backend%mesh%get_n_groups(DIR_Z)])

    allocate (backend%u_send_s_dev(SZ, n_halo, n_groups))
    allocate (backend%u_send_e_dev(SZ, n_halo, n_groups))
    allocate (backend%u_recv_s_dev(SZ, n_halo, n_groups))
    allocate (backend%u_recv_e_dev(SZ, n_halo, n_groups))
    allocate (backend%v_send_s_dev(SZ, n_halo, n_groups))
    allocate (backend%v_send_e_dev(SZ, n_halo, n_groups))
    allocate (backend%v_recv_s_dev(SZ, n_halo, n_groups))
    allocate (backend%v_recv_e_dev(SZ, n_halo, n_groups))
    allocate (backend%w_send_s_dev(SZ, n_halo, n_groups))
    allocate (backend%w_send_e_dev(SZ, n_halo, n_groups))
    allocate (backend%w_recv_s_dev(SZ, n_halo, n_groups))
    allocate (backend%w_recv_e_dev(SZ, n_halo, n_groups))

    allocate (backend%du_send_s_dev(SZ, 1, n_groups))
    allocate (backend%du_send_e_dev(SZ, 1, n_groups))
    allocate (backend%du_recv_s_dev(SZ, 1, n_groups))
    allocate (backend%du_recv_e_dev(SZ, 1, n_groups))
    allocate (backend%dud_send_s_dev(SZ, 1, n_groups))
    allocate (backend%dud_send_e_dev(SZ, 1, n_groups))
    allocate (backend%dud_recv_s_dev(SZ, 1, n_groups))
    allocate (backend%dud_recv_e_dev(SZ, 1, n_groups))
    allocate (backend%d2u_send_s_dev(SZ, 1, n_groups))
    allocate (backend%d2u_send_e_dev(SZ, 1, n_groups))
    allocate (backend%d2u_recv_s_dev(SZ, 1, n_groups))
    allocate (backend%d2u_recv_e_dev(SZ, 1, n_groups))

  end function init

  subroutine alloc_cuda_tdsops( &
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    )
    implicit none

    class(cuda_backend_t) :: self
    class(tdsops_t), allocatable, intent(inout) :: tdsops
    integer, intent(in) :: n_tds
    real(dp), intent(in) :: delta
    character(*), intent(in) :: operation, scheme
    integer, intent(in) :: bc_start, bc_end
    real(dp), optional, intent(in) :: stretch(:), stretch_correct(:)
    integer, optional, intent(in) :: n_halo
    character(*), optional, intent(in) :: from_to
    logical, optional, intent(in) :: sym
    real(dp), optional, intent(in) :: c_nu, nu0_nu

    allocate (cuda_tdsops_t :: tdsops)

    select type (tdsops)
    type is (cuda_tdsops_t)
      tdsops = cuda_tdsops_t(n_tds, delta, operation, scheme, bc_start, &
                             bc_end, stretch, stretch_correct, n_halo, &
                             from_to, sym, c_nu, nu0_nu)
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

    real(dp), device, pointer, dimension(:, :, :) :: u_dev, v_dev, w_dev, &
                                                     du_dev, dv_dev, dw_dev

    type(cuda_tdsops_t), pointer :: der1st, der1st_sym, der2nd, der2nd_sym

    call resolve_field_t(u_dev, u)
    call resolve_field_t(v_dev, v)
    call resolve_field_t(w_dev, w)

    call resolve_field_t(du_dev, du)
    call resolve_field_t(dv_dev, dv)
    call resolve_field_t(dw_dev, dw)

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

    call transeq_halo_exchange(self, u_dev, v_dev, w_dev, dirps%dir)

    call transeq_dist_component(self, du_dev, u_dev, u_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st, der1st_sym, der2nd, dirps%dir, &
                                blocks, threads)
    call transeq_dist_component(self, dv_dev, v_dev, u_dev, &
                                self%v_recv_s_dev, self%v_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st_sym, der1st, der2nd_sym, dirps%dir, &
                                blocks, threads)
    call transeq_dist_component(self, dw_dev, w_dev, u_dev, &
                                self%w_recv_s_dev, self%w_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st_sym, der1st, der2nd_sym, dirps%dir, &
                                blocks, threads)

    call du%set_data_loc(u%data_loc)
    call dv%set_data_loc(v%data_loc)
    call dw%set_data_loc(w%data_loc)

  end subroutine transeq_cuda_dist

  subroutine transeq_halo_exchange(self, u_dev, v_dev, w_dev, dir)
    class(cuda_backend_t) :: self
    real(dp), device, dimension(:, :, :), intent(in) :: u_dev, v_dev, w_dev
    integer, intent(in) :: dir
    integer :: n_halo, n, nproc_dir, pprev, pnext
    integer :: n_groups

    ! TODO: don't hardcode n_halo
    n_halo = 4
    n_groups = self%mesh%get_n_groups(dir)
    n = self%mesh%get_n(dir, VERT)
    nproc_dir = self%mesh%par%nproc_dir(dir)
    pprev = self%mesh%par%pprev(dir)
    pnext = self%mesh%par%pnext(dir)

    ! Copy halo data into buffer arrays
    call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, n)
    call copy_into_buffers(self%v_send_s_dev, self%v_send_e_dev, v_dev, n)
    call copy_into_buffers(self%w_send_s_dev, self%w_send_e_dev, w_dev, n)

    ! halo exchange
    call sendrecv_3fields( &
      self%u_recv_s_dev, self%u_recv_e_dev, &
      self%v_recv_s_dev, self%v_recv_e_dev, &
      self%w_recv_s_dev, self%w_recv_e_dev, &
      self%u_send_s_dev, self%u_send_e_dev, &
      self%v_send_s_dev, self%v_send_e_dev, &
      self%w_send_s_dev, self%w_send_e_dev, &
      SZ*n_halo*n_groups, nproc_dir, pprev, pnext)

  end subroutine transeq_halo_exchange

  subroutine transeq_dist_component(self, rhs_du_dev, u_dev, conv_dev, &
                                    u_recv_s_dev, u_recv_e_dev, &
                                    conv_recv_s_dev, conv_recv_e_dev, &
                                    tdsops_du, tdsops_dud, tdsops_d2u, &
                                    dir, blocks, threads)
    !! Computes RHS_x^u following:
    !!
    !! rhs_x^u = -0.5*(conv*du/dx + d(u*conv)/dx) + nu*d2u/dx2
    class(cuda_backend_t) :: self
    !> The result field, it is also used as temporary storage
    real(dp), device, dimension(:, :, :), intent(out) :: rhs_du_dev
    real(dp), device, dimension(:, :, :), intent(in) :: u_dev, conv_dev
    real(dp), device, dimension(:, :, :), intent(in) :: &
      u_recv_s_dev, u_recv_e_dev, &
      conv_recv_s_dev, conv_recv_e_dev
    class(cuda_tdsops_t), intent(in) :: tdsops_du, tdsops_dud, tdsops_d2u
    integer, intent(in) :: dir
    type(dim3), intent(in) :: blocks, threads

    class(field_t), pointer :: dud, d2u

    real(dp), device, pointer, dimension(:, :, :) :: dud_dev, d2u_dev

    ! Get some fields for storing the intermediate results
    dud => self%allocator%get_block(dir)
    d2u => self%allocator%get_block(dir)

    call resolve_field_t(dud_dev, dud)
    call resolve_field_t(d2u_dev, d2u)

    call exec_dist_transeq_3fused( &
      rhs_du_dev, &
      u_dev, u_recv_s_dev, u_recv_e_dev, &
      conv_dev, conv_recv_s_dev, conv_recv_e_dev, &
      dud_dev, d2u_dev, &
      self%du_send_s_dev, self%du_send_e_dev, &
      self%du_recv_s_dev, self%du_recv_e_dev, &
      self%dud_send_s_dev, self%dud_send_e_dev, &
      self%dud_recv_s_dev, self%dud_recv_e_dev, &
      self%d2u_send_s_dev, self%d2u_send_e_dev, &
      self%d2u_recv_s_dev, self%d2u_recv_e_dev, &
      tdsops_du, tdsops_d2u, self%nu, &
      self%mesh%par%nproc_dir(dir), self%mesh%par%pprev(dir), &
      self%mesh%par%pnext(dir), blocks, threads &
      )

    ! Release temporary blocks
    call self%allocator%release_block(dud)
    call self%allocator%release_block(d2u)

  end subroutine transeq_dist_component

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

  subroutine tds_solve_cuda(self, du, u, tdsops)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    type(dim3) :: blocks, threads

    ! Check if direction matches for both in/out fields and dirps
    if (u%dir /= du%dir) then
      error stop 'DIR mismatch between fields in tds_solve.'
    end if

    blocks = dim3(self%mesh%get_n_groups(u), 1, 1); threads = dim3(SZ, 1, 1)

    if (u%data_loc /= NULL_LOC) then
      call du%set_data_loc(move_data_loc(u%data_loc, u%dir, tdsops%move))
    end if

    call tds_solve_dist(self, du, u, tdsops, blocks, threads)

  end subroutine tds_solve_cuda

  subroutine tds_solve_dist(self, du, u, tdsops, blocks, threads)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops
    type(dim3), intent(in) :: blocks, threads

    real(dp), device, pointer, dimension(:, :, :) :: du_dev, u_dev

    type(cuda_tdsops_t), pointer :: tdsops_dev

    integer :: n_halo, n_groups, dir

    ! TODO: don't hardcode n_halo
    n_halo = 4
    dir = u%dir
    n_groups = self%mesh%get_n_groups(u)

    call resolve_field_t(du_dev, du)
    call resolve_field_t(u_dev, u)

    select type (tdsops)
    type is (cuda_tdsops_t); tdsops_dev => tdsops
    end select

    call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, &
                           tdsops_dev%n_tds)

    call sendrecv_fields(self%u_recv_s_dev, self%u_recv_e_dev, &
                         self%u_send_s_dev, self%u_send_e_dev, &
                         SZ*n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dir), &
                         self%mesh%par%pprev(dir), &
                         self%mesh%par%pnext(dir))

    ! call exec_dist
    call exec_dist_tds_compact( &
      du_dev, u_dev, &
      self%u_recv_s_dev, self%u_recv_e_dev, &
      self%du_send_s_dev, self%du_send_e_dev, &
      self%du_recv_s_dev, self%du_recv_e_dev, &
      tdsops_dev, self%mesh%par%nproc_dir(dir), &
      self%mesh%par%pprev(dir), self%mesh%par%pnext(dir), &
      blocks, threads &
      )

  end subroutine tds_solve_dist

  subroutine reorder_cuda(self, u_o, u_i, direction)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: u_o
    class(field_t), intent(in) :: u_i
    integer, intent(in) :: direction

    real(dp), device, pointer, dimension(:, :, :) :: u_o_d, u_i_d, u_temp_d
    class(field_t), pointer :: u_temp
    type(dim3) :: blocks, threads
    integer :: nx_padded, ny_padded, nz_padded
    integer, dimension(3) :: dims_padded

    call resolve_field_t(u_o_d, u_o)
    call resolve_field_t(u_i_d, u_i)

    dims_padded = self%mesh%get_padded_dims(DIR_C)
    nx_padded = dims_padded(1)
    ny_padded = dims_padded(2)
    nz_padded = dims_padded(3)

    select case (direction)
    case (RDR_X2Y)
      blocks = dim3(nx_padded/SZ, nz_padded, ny_padded/SZ)
      threads = dim3(SZ, SZ, 1)
      call reorder_x2y<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_X2Z)
      blocks = dim3(nx_padded, ny_padded/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_x2z<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_Y2X)
      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_y2x<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_Y2Z)
      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_y2z<<<blocks, threads>>>(u_o_d, u_i_d, & !&
                                            nx_padded, nz_padded)
    case (RDR_Z2X)
      blocks = dim3(nx_padded, ny_padded/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_z2x<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_Z2Y)
      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_z2y<<<blocks, threads>>>(u_o_d, u_i_d, & !&
                                            nx_padded, nz_padded)
    case (RDR_C2X)
      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_c2x<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_C2Y)
      ! First reorder from C to X, then from X to Y
      u_temp => self%allocator%get_block(DIR_X)
      call resolve_field_t(u_temp_d, u_temp)

      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_c2x<<<blocks, threads>>>(u_temp_d, u_i_d, nz_padded) !&

      blocks = dim3(nx_padded/SZ, nz_padded, ny_padded/SZ)
      threads = dim3(SZ, SZ, 1)
      call reorder_x2y<<<blocks, threads>>>(u_o_d, u_temp_d, nz_padded) !&

      call self%allocator%release_block(u_temp)
    case (RDR_C2Z)
      ! First reorder from C to X, then from X to Z
      u_temp => self%allocator%get_block(DIR_X)
      call resolve_field_t(u_temp_d, u_temp)

      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_c2x<<<blocks, threads>>>(u_temp_d, u_i_d, nz_padded) !&

      blocks = dim3(nx_padded, ny_padded/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_x2z<<<blocks, threads>>>(u_o_d, u_temp_d, nz_padded) !&

      call self%allocator%release_block(u_temp)
    case (RDR_X2C)
      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_x2c<<<blocks, threads>>>(u_o_d, u_i_d, nz_padded) !&
    case (RDR_Y2C)
      ! First reorder from Y to X, then from X to C
      u_temp => self%allocator%get_block(DIR_X)
      call resolve_field_t(u_temp_d, u_temp)

      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_y2x<<<blocks, threads>>>(u_temp_d, u_i_d, nz_padded) !&

      call reorder_x2c<<<blocks, threads>>>(u_o_d, u_temp_d, nz_padded) !&

      call self%allocator%release_block(u_temp)
    case (RDR_Z2C)
      ! First reorder from Z to X, then from X to C
      u_temp => self%allocator%get_block(DIR_X)
      call resolve_field_t(u_temp_d, u_temp)

      blocks = dim3(nx_padded, ny_padded/SZ, 1)
      threads = dim3(SZ, 1, 1)
      call reorder_z2x<<<blocks, threads>>>(u_temp_d, u_i_d, nz_padded) !&

      blocks = dim3(nx_padded/SZ, ny_padded/SZ, nz_padded)
      threads = dim3(SZ, SZ, 1)
      call reorder_x2c<<<blocks, threads>>>(u_o_d, u_temp_d, nz_padded) !&

      call self%allocator%release_block(u_temp)
    case default
      error stop 'Reorder direction is undefined.'
    end select

    ! reorder keeps the data_loc the same
    call u_o%set_data_loc(u_i%data_loc)

  end subroutine reorder_cuda

  subroutine sum_yintox_cuda(self, u, u_y)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_y

    real(dp), device, pointer, dimension(:, :, :) :: u_d, u_y_d
    type(dim3) :: blocks, threads
    integer, dimension(3) :: dims_padded

    call resolve_field_t(u_d, u)
    call resolve_field_t(u_y_d, u_y)

    dims_padded = self%mesh%get_padded_dims(DIR_C)

    blocks = dim3(dims_padded(1)/SZ, dims_padded(2)/SZ, dims_padded(3))
    threads = dim3(SZ, SZ, 1)
    call sum_yintox<<<blocks, threads>>>(u_d, u_y_d, dims_padded(3)) !&

  end subroutine sum_yintox_cuda

  subroutine sum_zintox_cuda(self, u, u_z)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_z

    real(dp), device, pointer, dimension(:, :, :) :: u_d, u_z_d
    type(dim3) :: blocks, threads
    integer, dimension(3) :: dims_padded

    call resolve_field_t(u_d, u)
    call resolve_field_t(u_z_d, u_z)

    dims_padded = self%mesh%get_padded_dims(DIR_C)

    blocks = dim3(dims_padded(1), dims_padded(2)/SZ, 1)
    threads = dim3(SZ, 1, 1)
    call sum_zintox<<<blocks, threads>>>(u_d, u_z_d, dims_padded(3)) !&

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

    call resolve_field_t(x_d, x)
    call resolve_field_t(y_d, y)

    nx = size(x_d, dim=2)
    blocks = dim3(size(x_d, dim=3), 1, 1)
    threads = dim3(SZ, 1, 1)
    call axpby<<<blocks, threads>>>(nx, a, x_d, b, y_d) !&

  end subroutine vecadd_cuda

  real(dp) function scalar_product_cuda(self, x, y) result(s)
    !! [[m_base_backend(module):scalar_product(interface)]]
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(in) :: x, y

    real(dp), device, pointer, dimension(:, :, :) :: x_d, y_d
    real(dp), device, allocatable :: sum_d
    integer :: dims(3), dims_padded(3), n, n_i, n_i_pad, n_j, ierr
    type(dim3) :: blocks, threads

    if ((x%data_loc == NULL_LOC) .or. (y%data_loc == NULL_LOC)) then
      error stop "You must set the data_loc before calling scalar product"
    end if
    if ((x%data_loc /= y%data_loc) .or. (x%dir /= y%dir)) then
      error stop "Called scalar product with incompatible fields"
    end if

    call resolve_field_t(x_d, x)
    call resolve_field_t(y_d, y)

    allocate (sum_d)
    sum_d = 0._dp

    dims = self%mesh%get_dims(x%data_loc)
    dims_padded = self%mesh%get_padded_dims(DIR_C)

    if (x%dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (x%dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (x%dir == DIR_Z) then
      n = dims(3); n_j = dims(2); n_i = dims(1); n_i_pad = dims_padded(1)
    else
      error stop 'scalar_product_cuda does not support DIR_C fields!'
    end if

    blocks = dim3(n_i, (n_j - 1)/SZ + 1, 1)
    threads = dim3(SZ, 1, 1)
    call scalar_product<<<blocks, threads>>>(sum_d, x_d, y_d, & !&
                                             n, n_i_pad, n_j)

    s = sum_d

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function scalar_product_cuda

  subroutine copy_into_buffers(u_send_s_dev, u_send_e_dev, u_dev, n)
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: u_send_s_dev, &
                                                         u_send_e_dev
    real(dp), device, dimension(:, :, :), intent(in) :: u_dev
    integer, intent(in) :: n

    type(dim3) :: blocks, threads
    integer :: n_halo = 4

    blocks = dim3(size(u_dev, dim=3), 1, 1)
    threads = dim3(SZ, 1, 1)
    call buffer_copy<<<blocks, threads>>>(u_send_s_dev, u_send_e_dev, & !&
                                          u_dev, n, n_halo)

  end subroutine copy_into_buffers

  subroutine field_max_mean_cuda(self, max_val, mean_val, f, enforced_data_loc)
    !! [[m_base_backend(module):field_max_mean(interface)]]
    implicit none

    class(cuda_backend_t) :: self
    real(dp), intent(out) :: max_val, mean_val
    class(field_t), intent(in) :: f
    integer, optional, intent(in) :: enforced_data_loc

    real(dp), device, pointer, dimension(:, :, :) :: f_d
    real(dp), device, allocatable :: max_d, sum_d
    integer :: data_loc, dims(3), dims_padded(3), n, n_i, n_i_pad, n_j, ierr
    type(dim3) :: blocks, threads

    if (f%data_loc == NULL_LOC .and. (.not. present(enforced_data_loc))) then
      error stop 'The input field to cuda::field_max_mean does not have a &
                  &valid f%data_loc. You may enforce a data_loc of your &
                  &choice as last argument to carry on at your own risk!'
    end if

    if (present(enforced_data_loc)) then
      data_loc = enforced_data_loc
    else
      data_loc = f%data_loc
    end if

    dims = self%mesh%get_dims(data_loc)
    dims_padded = self%mesh%get_padded_dims(DIR_C)

    call resolve_field_t(f_d, f)

    allocate (max_d, sum_d)
    max_d = 0._dp; sum_d = 0._dp

    if (f%dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Z) then
      n = dims(3); n_j = dims(2); n_i = dims(1); n_i_pad = dims_padded(1)
    else
      error stop 'field_max_mean does not support DIR_C fields!'
    end if

    blocks = dim3(n_i, (n_j - 1)/SZ + 1, 1)
    threads = dim3(SZ, 1, 1)
    call field_max_sum<<<blocks, threads>>>(max_d, sum_d, f_d, & !&
                                            n, n_i_pad, n_j)

    ! rank-local values, copy them first from device to host
    max_val = max_d
    mean_val = sum_d
    mean_val = mean_val/product(self%mesh%get_global_dims(data_loc))

    ! make sure all ranks have final values
    call MPI_Allreduce(MPI_IN_PLACE, max_val, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mean_val, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine field_max_mean_cuda

  subroutine field_scale_cuda(self, f, a)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    real(dp), device, pointer, dimension(:, :, :) :: f_d
    type(dim3) :: blocks, threads
    integer :: n

    call resolve_field_t(f_d, f)

    n = size(f_d, dim=2)
    blocks = dim3(size(f_d, dim=3), 1, 1)
    threads = dim3(SZ, 1, 1)
    call field_scale<<<blocks, threads>>>(f_d, a, n) !&

  end subroutine field_scale_cuda

  subroutine field_shift_cuda(self, f, a)
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    real(dp), device, pointer, dimension(:, :, :) :: f_d
    type(dim3) :: blocks, threads
    integer :: n

    call resolve_field_t(f_d, f)

    n = size(f_d, dim=2)
    blocks = dim3(size(f_d, dim=3), 1, 1)
    threads = dim3(SZ, 1, 1)
    call field_shift<<<blocks, threads>>>(f_d, a, n) !&

  end subroutine field_shift_cuda

  subroutine field_set_face_cuda(self, f, c_start, c_end, face)
    !! [[m_base_backend(module):field_set_face(subroutine)]]
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c_start, c_end
    integer, intent(in) :: face

    real(dp), device, pointer, dimension(:, :, :) :: f_d
    type(dim3) :: blocks, threads
    integer :: dims(3), nx, ny, nz

    if (f%dir /= DIR_X) then
      error stop 'Setting a field face is only supported for DIR_X fields.'
    end if

    if (f%data_loc == NULL_LOC) then
      error stop 'field_set_face require a valid data_loc.'
    end if

    call resolve_field_t(f_d, f)

    dims = self%mesh%get_dims(f%data_loc)

    select case (face)
    case (X_FACE)
      error stop 'Setting X_FACE is not yet supported.'
    case (Y_FACE)
      blocks = dim3((dims(1) - 1)/64 + 1, dims(3), 1)
      threads = dim3(64, 1, 1)
      call field_set_y_face<<<blocks, threads>>>(f_d, c_start, c_end, & !&
                                                 dims(1), dims(2), dims(3))
    case (Z_FACE)
      error stop 'Setting Z_FACE is not yet supported.'
    case default
      error stop 'face is undefined.'
    end select

  end subroutine field_set_face_cuda

  real(dp) function field_volume_integral_cuda(self, f) result(s)
    !! volume integral of a field
    implicit none

    class(cuda_backend_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer, dimension(:, :, :) :: f_d
    real(dp), device, allocatable :: integral_d
    integer :: dims(3), dims_padded(3), ierr
    type(dim3) :: blocks, threads

    if (f%data_loc == NULL_LOC) then
      error stop 'You must set the data_loc before calling volume integral.'
    end if
    if (f%dir /= DIR_X) then
      error stop 'Volume integral can only be called on DIR_X fields.'
    end if

    call resolve_field_t(f_d, f)

    allocate (integral_d)
    integral_d = 0._dp

    dims = self%mesh%get_dims(f%data_loc)
    dims_padded = self%mesh%get_padded_dims(DIR_C)

    blocks = dim3(dims(3), (dims(2) - 1)/SZ + 1, 1)
    threads = dim3(SZ, 1, 1)
    call volume_integral<<<blocks, threads>>>(integral_d, f_d, & !&
                                              dims(1), dims_padded(3), dims(2))

    s = integral_d

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function field_volume_integral_cuda

  subroutine copy_data_to_f_cuda(self, f, data)
    class(cuda_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(inout) :: data

    select type (f); type is (cuda_field_t); f%data_d = data; end select
  end subroutine copy_data_to_f_cuda

  subroutine copy_f_to_data_cuda(self, data, f)
    class(cuda_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    select type (f); type is (cuda_field_t); data = f%data_d; end select
  end subroutine copy_f_to_data_cuda

  subroutine init_cuda_poisson_fft(self, mesh, xdirps, ydirps, zdirps)
    implicit none

    class(cuda_backend_t) :: self
    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    allocate (cuda_poisson_fft_t :: self%poisson_fft)

    select type (poisson_fft => self%poisson_fft)
    type is (cuda_poisson_fft_t)
      poisson_fft = cuda_poisson_fft_t(mesh, xdirps, ydirps, zdirps)
    end select

  end subroutine init_cuda_poisson_fft

  subroutine resolve_field_t(u_dev, u)
    real(dp), device, pointer, dimension(:, :, :), intent(out) :: u_dev
    class(field_t), intent(in) :: u

    select type (u)
    type is (cuda_field_t)
      u_dev => u%data_d
    end select

  end subroutine resolve_field_t

end module m_cuda_backend

