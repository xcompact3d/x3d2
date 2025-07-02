module m_omp_backend
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, get_dirs_from_rdr, move_data_loc, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, NULL_LOC, &
                      X_FACE, Y_FACE, Z_FACE, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_ordering, only: get_index_reordering
  use m_tdsops, only: dirps_t, tdsops_t

  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_tds_compact, exec_dist_transeq_compact
  use m_omp_sendrecv, only: sendrecv_fields

  implicit none

  private :: transeq_halo_exchange, transeq_dist_component

  type, extends(base_backend_t) :: omp_backend_t
    !character(len=*), parameter :: name = 'omp'
    integer :: MPI_FP_PREC = dp
    real(dp), allocatable, dimension(:, :, :) :: &
      u_recv_s, u_recv_e, u_send_s, u_send_e, &
      v_recv_s, v_recv_e, v_send_s, v_send_e, &
      w_recv_s, w_recv_e, w_send_s, w_send_e, &
      spec_recv_s, spec_recv_e, spec_send_s, spec_send_e, &
      du_send_s, du_send_e, du_recv_s, du_recv_e, &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e
  contains
    procedure :: alloc_tdsops => alloc_omp_tdsops
    procedure :: transeq_x => transeq_x_omp
    procedure :: transeq_y => transeq_y_omp
    procedure :: transeq_z => transeq_z_omp
    procedure :: transeq_species => transeq_species_omp
    procedure :: tds_solve => tds_solve_omp
    procedure :: reorder => reorder_omp
    procedure :: sum_yintox => sum_yintox_omp
    procedure :: sum_zintox => sum_zintox_omp
    procedure :: veccopy => veccopy_omp
    procedure :: vecadd => vecadd_omp
    procedure :: vecmult => vecmult_omp
    procedure :: scalar_product => scalar_product_omp
    procedure :: field_max_mean => field_max_mean_omp
    procedure :: field_scale => field_scale_omp
    procedure :: field_shift => field_shift_omp
    procedure :: field_set_face => field_set_face_omp
    procedure :: field_volume_integral => field_volume_integral_omp
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

    type(mesh_t), target, intent(inout) :: mesh
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
    allocate (backend%spec_send_s(SZ, n_halo, n_groups))
    allocate (backend%spec_send_e(SZ, n_halo, n_groups))
    allocate (backend%spec_recv_s(SZ, n_halo, n_groups))
    allocate (backend%spec_recv_e(SZ, n_halo, n_groups))

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
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    )
    implicit none

    class(omp_backend_t) :: self
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

    allocate (tdsops_t :: tdsops)

    select type (tdsops)
    type is (tdsops_t)
      tdsops = tdsops_t(n_tds, delta, operation, scheme, bc_start, bc_end, &
                        stretch, stretch_correct, n_halo, from_to, sym, &
                        c_nu, nu0_nu)
    end select

  end subroutine alloc_omp_tdsops

  subroutine transeq_x_omp(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call self%transeq_omp_dist(du, dv, dw, u, v, w, nu, dirps)

  end subroutine transeq_x_omp

  subroutine transeq_y_omp(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    ! u, v, w is reordered so that we pass v, u, w
    call self%transeq_omp_dist(dv, du, dw, v, u, w, nu, dirps)

  end subroutine transeq_y_omp

  subroutine transeq_z_omp(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    ! u, v, w is reordered so that we pass w, u, v
    call self%transeq_omp_dist(dw, du, dv, w, u, v, nu, dirps)

  end subroutine transeq_z_omp

  subroutine transeq_species_omp(self, dspec, u, v, w, spec, nu, dirps, sync)
    !! Compute the convection and diffusion for the given field
    !! in the given direction.
    !! Halo exchange for the given field is necessary
    !! When sync is true, halo exchange of momentum is necessary
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: dspec
    class(field_t), intent(in) :: u, v, w, spec
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps
    logical, intent(in) :: sync

    integer :: n_halo, n_groups

    ! TODO: don't hardcode n_halo
    n_halo = 4
    n_groups = self%mesh%get_n_groups(dirps%dir)

    ! Halo exchange for momentum if needed
    if (sync) call transeq_halo_exchange(self, u, v, w, dirps%dir)

    ! Halo exchange for the given field
    call copy_into_buffers(self%spec_send_s, self%spec_send_e, spec%data, &
                           dirps%der1st%n_tds, n_groups)
    call sendrecv_fields(self%spec_recv_s, self%spec_recv_e, &
                         self%spec_send_s, self%spec_send_e, &
                         SZ*n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dirps%dir), &
                         self%mesh%par%pprev(dirps%dir), &
                         self%mesh%par%pnext(dirps%dir))

    ! combine convection and diffusion
    if (dirps%dir == DIR_X) then
      call transeq_dist_component(self, dspec, spec, u, nu, &
                                  self%spec_recv_s, self%spec_recv_e, &
                                  self%u_recv_s, self%u_recv_e, &
                                  dirps%der1st, dirps%der1st_sym, &
                                  dirps%der2nd, dirps%dir)
    else if (dirps%dir == DIR_Y) then
      call transeq_dist_component(self, dspec, spec, v, nu, &
                                  self%spec_recv_s, self%spec_recv_e, &
                                  self%v_recv_s, self%v_recv_e, &
                                  dirps%der1st, dirps%der1st_sym, &
                                  dirps%der2nd, dirps%dir)
    else
      call transeq_dist_component(self, dspec, spec, w, nu, &
                                  self%spec_recv_s, self%spec_recv_e, &
                                  self%w_recv_s, self%w_recv_e, &
                                  dirps%der1st, dirps%der1st_sym, &
                                  dirps%der2nd, dirps%dir)
    end if

  end subroutine transeq_species_omp

  subroutine transeq_omp_dist(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call transeq_halo_exchange(self, u, v, w, dirps%dir)

    call transeq_dist_component(self, du, u, u, nu, &
                                self%u_recv_s, self%u_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st, dirps%der1st_sym, &
                                dirps%der2nd, dirps%dir)
    call transeq_dist_component(self, dv, v, u, nu, &
                                self%v_recv_s, self%v_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st_sym, dirps%der1st, &
                                dirps%der2nd_sym, dirps%dir)
    call transeq_dist_component(self, dw, w, u, nu, &
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

  subroutine transeq_dist_component(self, rhs_du, u, conv, nu, &
                                    u_recv_s, u_recv_e, &
                                    conv_recv_s, conv_recv_e, &
                                    tdsops_du, tdsops_dud, tdsops_d2u, dir)
    !! Computes RHS_x^u following:
    !!
    !! rhs_x^u = -0.5*(conv*du/dx + d(u*conv)/dx) + nu*d2u/dx2
    class(omp_backend_t) :: self
    !> The result field, it is also used as temporary storage
    class(field_t), intent(inout) :: rhs_du
    class(field_t), intent(in) :: u, conv
    real(dp), intent(in) :: nu
    real(dp), dimension(:, :, :), intent(in) :: u_recv_s, u_recv_e, &
                                                conv_recv_s, conv_recv_e
    class(tdsops_t), intent(in) :: tdsops_du
    class(tdsops_t), intent(in) :: tdsops_dud
    class(tdsops_t), intent(in) :: tdsops_d2u
    integer, intent(in) :: dir
    class(field_t), pointer :: d2u, dud

    dud => self%allocator%get_block(dir)
    d2u => self%allocator%get_block(dir)

    call exec_dist_transeq_compact( &
      rhs_du%data, dud%data, d2u%data, &
      self%du_send_s, self%du_send_e, self%du_recv_s, self%du_recv_e, &
      self%dud_send_s, self%dud_send_e, self%dud_recv_s, self%dud_recv_e, &
      self%d2u_send_s, self%d2u_send_e, self%d2u_recv_s, self%d2u_recv_e, &
      u%data, u_recv_s, u_recv_e, &
      conv%data, conv_recv_s, conv_recv_e, &
      tdsops_du, tdsops_dud, tdsops_d2u, nu, &
      self%mesh%par%nproc_dir(dir), self%mesh%par%pprev(dir), &
      self%mesh%par%pnext(dir), self%mesh%get_n_groups(dir))

    call rhs_du%set_data_loc(u%data_loc)

    call self%allocator%release_block(dud)
    call self%allocator%release_block(d2u)

  end subroutine transeq_dist_component

  subroutine tds_solve_omp(self, du, u, tdsops)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    ! Check if direction matches for both in/out fields
    if (u%dir /= du%dir) then
      error stop 'DIR mismatch between fields in tds_solve.'
    end if

    if (u%data_loc /= NULL_LOC) then
      call du%set_data_loc(move_data_loc(u%data_loc, u%dir, tdsops%move))
    end if

    call tds_solve_dist(self, du, u, tdsops)

  end subroutine tds_solve_omp

  subroutine tds_solve_dist(self, du, u, tdsops)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops
    integer :: n_halo, n_groups, dir

    ! TODO: don't hardcode n_halo
    n_halo = 4
    dir = u%dir
    n_groups = self%mesh%get_n_groups(u)

    call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, &
                           tdsops%n_tds, n_groups)

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

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omp

  subroutine sum_yintox_omp(self, u, u_)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call sum_intox_omp(self, u, u_, DIR_Y)

  end subroutine sum_yintox_omp

  subroutine sum_zintox_omp(self, u, u_)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call sum_intox_omp(self, u, u_, DIR_Z)

  end subroutine sum_zintox_omp

  subroutine sum_intox_omp(self, u, u_, dir_to)

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_
    integer, intent(in) :: dir_to

    integer :: dir_from
    integer, dimension(3) :: dims
    integer :: i, j, k    ! Working indices
    integer :: ii, jj, kk ! Transpose indices

    dir_from = DIR_X

    dims = self%mesh%get_padded_dims(u)
    !$omp parallel do private(i, ii, jj, kk) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(ii, jj, kk, i, j, k, &
                                    dir_from, dir_to, self%mesh)
          u%data(i, j, k) = u%data(i, j, k) + u_%data(ii, jj, kk)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine sum_intox_omp

  subroutine veccopy_omp(self, dst, src)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: dst
    class(field_t), intent(in) :: src
    integer, dimension(3) :: dims
    integer :: i, j, k, ii

    integer :: nvec, remstart

    if (src%dir /= dst%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    dims = self%mesh%get_padded_dims(src)
    nvec = dims(1)/SZ
    remstart = nvec*SZ + 1

    !$omp parallel do private(i, ii) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        ! Execute inner vectorised loops
        do ii = 1, nvec
          !$omp simd
          do i = 1, SZ
            dst%data(i + (ii - 1)*SZ, j, k) = &
              src%data(i + (ii - 1)*SZ, j, k)
          end do
          !$omp end simd
        end do

        ! Remainder loop
        do i = remstart, dims(1)
          dst%data(i, j, k) = src%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine veccopy_omp

  subroutine vecadd_omp(self, a, x, b, y)
    implicit none

    class(omp_backend_t) :: self
    real(dp), intent(in) :: a
    class(field_t), intent(in) :: x
    real(dp), intent(in) :: b
    class(field_t), intent(inout) :: y
    integer, dimension(3) :: dims
    integer :: i, j, k, ii

    integer :: nvec, remstart

    if (x%dir /= y%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    dims = self%mesh%get_padded_dims(x)
    nvec = dims(1)/SZ
    remstart = nvec*SZ + 1

    !$omp parallel do private(i, ii) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        ! Execute inner vectorised loops
        do ii = 1, nvec
          !$omp simd
          do i = 1, SZ
            y%data(i + (ii - 1)*SZ, j, k) = &
              a*x%data(i + (ii - 1)*SZ, j, k) + &
              b*y%data(i + (ii - 1)*SZ, j, k)
          end do
          !$omp end simd
        end do

        ! Remainder loop
        do i = remstart, dims(1)
          y%data(i, j, k) = a*x%data(i, j, k) + b*y%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine vecadd_omp

  subroutine vecmult_omp(self, y, x)
    !! [[m_base_backend(module):vecmult(interface)]]
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: y
    class(field_t), intent(in) :: x
    integer :: i, j, k

    !$omp parallel do
    do k = 1, size(y%data, 3)
      do j = 1, size(y%data, 2)
        !$omp simd
        do i = 1, SZ
          y%data(i, j, k) = y%data(i, j, k)*x%data(i, j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine vecmult_omp

  real(dp) function scalar_product_omp(self, x, y) result(s)
    !! [[m_base_backend(module):scalar_product(interface)]]
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(in) :: x, y
    class(field_t), pointer :: x_, y_
    integer, dimension(3) :: dims
    integer :: i, j, k, ii
    integer :: nvec, remstart
    integer :: ierr

    if ((x%data_loc == NULL_LOC) .or. (y%data_loc == NULL_LOC)) then
      error stop "You must set the data_loc before calling scalar product"
    end if
    if (x%data_loc /= y%data_loc) then
      error stop "Called scalar product with incompatible fields"
    end if

    ! Reorient data into temporary DIR_C storage
    x_ => self%allocator%get_block(DIR_C, x%data_loc)
    call self%get_field_data(x_%data, x)
    y_ => self%allocator%get_block(DIR_C, y%data_loc)
    call self%get_field_data(y_%data, y)

    dims = self%mesh%get_field_dims(x_)

    nvec = dims(1)/SZ
    remstart = nvec*SZ + 1

    s = 0.0_dp
    !$omp parallel do reduction(+:s) private(i, ii) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        ! Execute inner vectorised loops
        do ii = 1, nvec
          !$omp simd reduction(+:s)
          do i = 1, SZ
            s = s + x_%data(i + (ii - 1)*SZ, j, k)* &
                y_%data(i + (ii - 1)*SZ, j, k)
          end do
          !$omp end simd
        end do

        ! Remainder loop
        do i = remstart, dims(1)
          s = s + x_%data(i, j, k)*y_%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

    ! Release temporary storage
    call self%allocator%release_block(x_)
    call self%allocator%release_block(y_)

    ! Reduce the result
    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, &
                       ierr)

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

  subroutine field_max_mean_omp(self, max_val, mean_val, f, enforced_data_loc)
    !! [[m_base_backend(module):field_max_mean(interface)]]
    implicit none

    class(omp_backend_t) :: self
    real(dp), intent(out) :: max_val, mean_val
    class(field_t), intent(in) :: f
    integer, optional, intent(in) :: enforced_data_loc

    real(dp) :: val, max_p, sum_p, max_pncl, sum_pncl
    integer :: data_loc, dims(3), dims_padded(3), n, n_i, n_i_pad, n_j
    integer :: i, j, k, k_i, k_j, ierr

    if (f%data_loc == NULL_LOC .and. (.not. present(enforced_data_loc))) then
      error stop 'The input field to omp::field_max_mean does not have a &
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

    if (f%dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Z) then
      n = dims(3); n_j = dims(1); n_i = dims(2); n_i_pad = dims_padded(2)
    else
      error stop 'field_max_mean does not support DIR_C fields!'
    end if

    sum_p = 0._dp
    max_p = 0._dp
    !$omp parallel do collapse(2) reduction(+:sum_p) reduction(max:max_p) &
    !$omp private(k, val, sum_pncl, max_pncl)
    do k_j = 1, (n_j - 1)/SZ + 1 ! loop over stacked groups
      do k_i = 1, n_i
        k = k_j + (k_i - 1)*((n_j - 1)/SZ + 1)
        sum_pncl = 0._dp
        max_pncl = 0._dp
        do j = 1, n
          ! loop over only non-padded entries in the present group
          do i = 1, min(SZ, n_j - (k_j - 1)*SZ)
            val = abs(f%data(i, j, k))
            sum_pncl = sum_pncl + val
            max_pncl = max(max_pncl, val)
          end do
        end do
        sum_p = sum_p + sum_pncl
        max_p = max(max_p, max_pncl)
      end do
    end do
    !$omp end parallel do

    ! rank-local values
    max_val = max_p
    mean_val = sum_p/product(self%mesh%get_global_dims(data_loc))

    ! make sure all ranks have final values
    call MPI_Allreduce(MPI_IN_PLACE, max_val, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mean_val, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine field_max_mean_omp

  subroutine field_scale_omp(self, f, a)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    f%data = a*f%data
  end subroutine field_scale_omp

  subroutine field_shift_omp(self, f, a)
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    f%data = f%data + a
  end subroutine field_shift_omp

  subroutine field_set_face_omp(self, f, c_start, c_end, face)
    !! [[m_base_backend(module):field_set_face(subroutine)]]
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c_start, c_end
    integer, intent(in) :: face

    integer :: dims(3), k, j, i_mod, k_end

    if (f%dir /= DIR_X) then
      error stop 'Setting a field face is only supported for DIR_X fields.'
    end if

    if (f%data_loc == NULL_LOC) then
      error stop 'field_set_face require a valid data_loc.'
    end if

    dims = self%mesh%get_dims(f%data_loc)

    select case (face)
    case (X_FACE)
      error stop 'Setting X_FACE is not yet supported.'
    case (Y_FACE)
      i_mod = mod(dims(2) - 1, SZ) + 1
      !$omp parallel do private(k_end)
      do k = 1, dims(3)
        k_end = k + (dims(2) - 1)/SZ*dims(3)
        do j = 1, dims(1)
          f%data(1, j, k) = c_start
          f%data(i_mod, j, k_end) = c_end
        end do
      end do
      !$omp end parallel do
    case (Z_FACE)
      error stop 'Setting Z_FACE is not yet supported.'
    case default
      error stop 'face is undefined.'
    end select

  end subroutine field_set_face_omp

  real(dp) function field_volume_integral_omp(self, f) result(s)
    !! volume integral of a field
    implicit none

    class(omp_backend_t) :: self
    class(field_t), intent(in) :: f

    real(dp) :: sum_p, sum_pncl
    integer :: dims(3), stacked, i, j, k, k_i, k_j, ierr

    if (f%data_loc == NULL_LOC) then
      error stop 'You must set the data_loc before calling volume integral.'
    end if
    if (f%dir /= DIR_X) then
      error stop 'Volume integral can only be called on DIR_X fields.'
    end if

    dims = self%mesh%get_dims(f%data_loc)
    stacked = (dims(2) - 1)/SZ + 1

    sum_p = 0._dp
    !$omp parallel do collapse(2) reduction(+:sum_p) private(k, sum_pncl)
    do k_j = 1, stacked ! loop over stacked groups
      do k_i = 1, dims(3)
        k = k_j + (k_i - 1)*stacked
        sum_pncl = 0._dp
        do j = 1, dims(1)
          ! loop over only non-padded entries in the present group
          do i = 1, min(SZ, dims(2) - (k_j - 1)*SZ)
            sum_pncl = sum_pncl + f%data(i, j, k)
          end do
        end do
        sum_p = sum_p + sum_pncl
      end do
    end do
    !$omp end parallel do

    ! rank-local values
    s = sum_p

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function field_volume_integral_omp

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
#ifdef WITH_2DECOMPFFT
    use m_omp_poisson_fft, only: omp_poisson_fft_t
#endif

    implicit none

    class(omp_backend_t) :: self
    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

#ifdef WITH_2DECOMPFFT
    allocate (omp_poisson_fft_t :: self%poisson_fft)

    select type (poisson_fft => self%poisson_fft)
    type is (omp_poisson_fft_t)
      poisson_fft = omp_poisson_fft_t(mesh, xdirps, ydirps, zdirps)
    end select
#else
    error stop 'This build does not support FFT based Poisson solver &
                &on the OpenMP backend!'
#endif

  end subroutine init_omp_poisson_fft

end module m_omp_backend

