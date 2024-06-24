module m_omp_backend
  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_ordering, only: get_index_reordering
  use m_common, only: dp, globs_t, get_dirs_from_rdr, VERT, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y
  use m_tdsops, only: dirps_t, tdsops_t, get_tds_n
  use m_omp_exec_dist, only: exec_dist_tds_compact, exec_dist_transeq_compact
  use m_omp_sendrecv, only: sendrecv_fields

  use m_omp_common, only: SZ
  use m_omp_poisson_fft, only: omp_poisson_fft_t
  use m_mesh, only: mesh_t

  implicit none

  private :: transeq_halo_exchange, transeq_dist_component

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
    procedure :: scalar_product => scalar_product_omp
    procedure :: copy_data_to_f => copy_data_to_f_omp
    procedure :: copy_f_to_data => copy_f_to_data_omp
    procedure :: init_poisson_fft => init_omp_poisson_fft
    procedure :: transeq_omp_dist
  end type omp_backend_t

  interface omp_backend_t
    module procedure init
  end interface omp_backend_t

contains

  function init(mesh, allocator) result(backend)
    implicit none

    class(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator
    type(omp_backend_t) :: backend

    integer :: n_halo, n_groups

    call backend%base_init()

    select type (allocator)
    type is (allocator_t)
      ! class level access to the allocator
      backend%allocator => allocator
    end select

    n_halo = 4
    backend%mesh => mesh
    n_groups = maxval([backend%mesh%get_n_groups(DIR_X), &
                       backend%mesh%get_n_groups(DIR_Y), &
                       backend%mesh%get_n_groups(DIR_Z)])

    allocate (backend%u_send_s(SZ, n_halo, n_groups))
    allocate (backend%u_send_e(SZ, n_halo, n_groups))
    allocate (backend%u_recv_s(SZ, n_halo, n_groups))
    allocate (backend%u_recv_e(SZ, n_halo, n_groups))
    allocate (backend%v_send_s(SZ, n_halo, n_groups))
    allocate (backend%v_send_e(SZ, n_halo, n_groups))
    allocate (backend%v_recv_s(SZ, n_halo, n_groups))
    allocate (backend%v_recv_e(SZ, n_halo, n_groups))
    allocate (backend%w_send_s(SZ, n_halo, n_groups))
    allocate (backend%w_send_e(SZ, n_halo, n_groups))
    allocate (backend%w_recv_s(SZ, n_halo, n_groups))
    allocate (backend%w_recv_e(SZ, n_halo, n_groups))

    allocate (backend%du_send_s(SZ, 1, n_groups))
    allocate (backend%du_send_e(SZ, 1, n_groups))
    allocate (backend%du_recv_s(SZ, 1, n_groups))
    allocate (backend%du_recv_e(SZ, 1, n_groups))
    allocate (backend%dud_send_s(SZ, 1, n_groups))
    allocate (backend%dud_send_e(SZ, 1, n_groups))
    allocate (backend%dud_recv_s(SZ, 1, n_groups))
    allocate (backend%dud_recv_e(SZ, 1, n_groups))
    allocate (backend%d2u_send_s(SZ, 1, n_groups))
    allocate (backend%d2u_send_e(SZ, 1, n_groups))
    allocate (backend%d2u_recv_s(SZ, 1, n_groups))
    allocate (backend%d2u_recv_e(SZ, 1, n_groups))

  end function init

  subroutine alloc_omp_tdsops( &
    self, tdsops, dir, operation, scheme, &
    n_halo, from_to, bc_start, bc_end, sym, c_nu, nu0_nu &
    )
    implicit none

    class(omp_backend_t) :: self
    class(tdsops_t), allocatable, intent(inout) :: tdsops
    integer, intent(in) :: dir
    character(*), intent(in) :: operation, scheme
    integer, optional, intent(in) :: n_halo
    character(*), optional, intent(in) :: from_to, bc_start, bc_end
    logical, optional, intent(in) :: sym
    real(dp), optional, intent(in) :: c_nu, nu0_nu
    integer :: tds_n
    real(dp) :: delta

    allocate (tdsops_t :: tdsops)

    select type (tdsops)
    type is (tdsops_t)
      tds_n = get_tds_n(self%mesh, dir, from_to)
      delta = self%mesh%geo%d(dir)
      tdsops = tdsops_t(tds_n, delta, operation, scheme, n_halo, from_to, &
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

    call transeq_halo_exchange(self, u, v, w, dirps%dir)

    call transeq_dist_component(self, du, u, u, &
                                self%u_recv_s, self%u_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st, dirps%der1st_sym, &
                                dirps%der2nd, dirps%dir)
    call transeq_dist_component(self, dv, v, u, &
                                self%v_recv_s, self%v_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st_sym, dirps%der1st, &
                                dirps%der2nd_sym, dirps%dir)
    call transeq_dist_component(self, dw, w, u, &
                                self%w_recv_s, self%w_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st_sym, dirps%der1st, &
                                dirps%der2nd_sym, dirps%dir)

  end subroutine transeq_omp_dist

  subroutine transeq_halo_exchange(self, u, v, w, dir)
    class(omp_backend_t) :: self
    class(field_t), intent(in) :: u, v, w
    integer, intent(in) :: dir
    integer :: n_halo, n, nproc_dir, pprev, pnext
    integer :: n_groups

    ! TODO: don't hardcode n_halo
    n_halo = 4
    n_groups = self%mesh%get_n_groups(dir)
    n = self%mesh%get_n(u)
    nproc_dir = self%mesh%par%nproc_dir(dir)
    pprev = self%mesh%par%pprev(dir)
    pnext = self%mesh%par%pnext(dir)

    call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, &
                           n, n_groups)
    call copy_into_buffers(self%v_send_s, self%v_send_e, v%data, &
                           n, n_groups)
    call copy_into_buffers(self%w_send_s, self%w_send_e, w%data, &
                           n, n_groups)

    call sendrecv_fields(self%u_recv_s, self%u_recv_e, &
                         self%u_send_s, self%u_send_e, &
                         SZ*n_halo*n_groups, &
                         nproc_dir, pprev, pnext)
    call sendrecv_fields(self%v_recv_s, self%v_recv_e, &
                         self%v_send_s, self%v_send_e, &
                         SZ*n_halo*n_groups, &
                         nproc_dir, pprev, pnext)
    call sendrecv_fields(self%w_recv_s, self%w_recv_e, &
                         self%w_send_s, self%w_send_e, &
                         SZ*n_halo*n_groups, &
                         nproc_dir, pprev, pnext)

  end subroutine transeq_halo_exchange

  subroutine transeq_dist_component(self, rhs, u, conv, &
                                    u_recv_s, u_recv_e, &
                                    conv_recv_s, conv_recv_e, &
                                    tdsops_du, tdsops_dud, tdsops_d2u, dir)
      !! Computes RHS_x^u following:
      !!
      !! rhs_x^u = -0.5*(conv*du/dx + d(u*conv)/dx) + nu*d2u/dx2
    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: rhs
    class(field_t), intent(in) :: u, conv
    real(dp), dimension(:, :, :), intent(in) :: u_recv_s, u_recv_e, &
                                                conv_recv_s, conv_recv_e
    class(tdsops_t), intent(in) :: tdsops_du
    class(tdsops_t), intent(in) :: tdsops_dud
    class(tdsops_t), intent(in) :: tdsops_d2u
    integer, intent(in) :: dir
    class(field_t), pointer :: du, d2u, dud

    du => self%allocator%get_block(dir, VERT)
    dud => self%allocator%get_block(dir, VERT)
    d2u => self%allocator%get_block(dir, VERT)

    call exec_dist_transeq_compact( &
      rhs%data, du%data, dud%data, d2u%data, &
      self%du_send_s, self%du_send_e, self%du_recv_s, self%du_recv_e, &
      self%dud_send_s, self%dud_send_e, self%dud_recv_s, self%dud_recv_e, &
      self%d2u_send_s, self%d2u_send_e, self%d2u_recv_s, self%d2u_recv_e, &
      u%data, u_recv_s, u_recv_e, &
      conv%data, conv_recv_s, conv_recv_e, &
      tdsops_du, tdsops_dud, tdsops_d2u, self%nu, &
      self%mesh%par%nproc_dir(dir), self%mesh%par%pprev(dir), &
      self%mesh%par%pnext(dir), self%mesh%get_n_groups(dir))

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

    ! Check if direction matches for both in/out fields and dirps
    if (dirps%dir /= du%dir .or. u%dir /= du%dir) then
      error stop 'DIR mismatch between fields and dirps in tds_solve.'
    end if

    call tds_solve_dist(self, du, u, dirps, tdsops)

  end subroutine tds_solve_omp

  subroutine tds_solve_dist(self, du, u, dirps, tdsops)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    type(dirps_t), intent(in) :: dirps
    class(tdsops_t), intent(in) :: tdsops
    integer :: n_halo, n_groups, dir

    ! TODO: don't hardcode n_halo
    n_halo = 4
    dir = u%dir
    n_groups = self%mesh%get_n_groups(u)

    call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, &
                           tdsops%tds_n, n_groups)

    ! halo exchange
    call sendrecv_fields(self%u_recv_s, self%u_recv_e, &
                         self%u_send_s, self%u_send_e, &
                         SZ*n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dir), &
                         self%mesh%par%pprev(dir), &
                         self%mesh%par%pnext(dir))

    call exec_dist_tds_compact( &
      du%data, u%data, self%u_recv_s, self%u_recv_e, &
      self%du_send_s, self%du_send_e, self%du_recv_s, self%du_recv_e, &
      tdsops, self%mesh%par%nproc_dir(dir), &
      self%mesh%par%pprev(dir), self%mesh%par%pnext(dir), &
      n_groups)

  end subroutine tds_solve_dist

  subroutine reorder_omp(self, u_, u, direction)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, intent(in) :: direction
    integer, dimension(3) :: dims
    integer :: i, j, k
    integer :: out_i, out_j, out_k
    integer :: dir_from, dir_to

    dims = self%mesh%get_padded_dims(u)
    call get_dirs_from_rdr(dir_from, dir_to, direction)

    !$omp parallel do private(out_i, out_j, out_k) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering( &
            out_i, out_j, out_k, i, j, k, dir_from, dir_to, self%mesh)
          u_%data(out_i, out_j, out_k) = u%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

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

    y%data = a*x%data + b*y%data ! fixme

  end subroutine vecadd_omp

  real(dp) function scalar_product_omp(self, x, y) result(s)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(in) :: x, y

    s = 0._dp

  end function scalar_product_omp

  subroutine copy_into_buffers(u_send_s, u_send_e, u, n, n_groups)
    implicit none

    real(dp), dimension(:, :, :), intent(out) :: u_send_s, u_send_e
    real(dp), dimension(:, :, :), intent(in) :: u
    integer, intent(in) :: n
    integer, intent(in) :: n_groups
    integer :: i, j, k
    integer :: n_halo = 4

    !$omp parallel do
    do k = 1, n_groups
      do j = 1, n_halo
        !$omp simd
        do i = 1, SZ
          u_send_s(i, j, k) = u(i, j, k)
          u_send_e(i, j, k) = u(i, n - n_halo + j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine copy_into_buffers

  subroutine copy_data_to_f_omp(self, f, data)
    class(omp_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    f%data = data
  end subroutine copy_data_to_f_omp

  subroutine copy_f_to_data_omp(self, data, f)
    class(omp_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    data = f%data
  end subroutine copy_f_to_data_omp

  subroutine init_omp_poisson_fft(self, mesh, xdirps, ydirps, zdirps)
    implicit none

    class(omp_backend_t) :: self
    class(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    allocate (omp_poisson_fft_t :: self%poisson_fft)

    select type (poisson_fft => self%poisson_fft)
    type is (omp_poisson_fft_t)
      poisson_fft = omp_poisson_fft_t(mesh, xdirps, ydirps, zdirps)
    end select

  end subroutine init_omp_poisson_fft

end module m_omp_backend

