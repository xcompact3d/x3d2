!!! src/backend/omptgt/backend.f90
!!
!! OpenMP target offload backend implementation, following the CUDA backend.
!!
!! Fields are device-resident: `omptgt_allocator_t` hands out
!! `omptgt_field_t`, whose storage is a raw `omp_target_alloc` pointer with no
!! host array behind it. `resolve_field_t` views that storage as a contiguous
!! Fortran pointer, which is then passed to a kernel taking it as an
!! explicit-shape dummy named in `OMP_DEVICE_ADDR` (see m_omptgt_common), so a
!! launch moves no field data between host and device.
!!
!! Not implemented, and so inherited from `base_backend_t` (whose defaults
!! stop the run naming the operation), are the operations the CUDA backend
!! builds on external libraries: the FFT Poisson solver (cuFFT/cuFFTMp) and
!! the device field export used by GPU-aware ADIOS2.

module m_omptgt_backend

  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_MAX, MPI_SUM, &
                   MPI_Allreduce

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, MPI_X3D2_DP, move_data_loc, get_dirs_from_rdr, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, NULL_LOC, &
                      X_FACE, Y_FACE, Z_FACE, BC_DIRICHLET
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t, tdsops_t

  use omp_lib, only: omp_get_default_device

  use m_omptgt_device_info, only: query_device_memory
  use m_omptgt_common, only: SZ, dev_alloc, dev_free, dev_upload, &
                             dev_download, device_sync
  use m_omptgt_exec_dist, only: exec_dist_transeq_3fused, &
                                exec_dist_tds_compact
  use m_omptgt_exec_thom, only: exec_thom_tds_compact
  use m_omptgt_sendrecv, only: sendrecv_fields, sendrecv_3fields
  use m_omptgt_tdsops, only: omptgt_tdsops_t, omptgt_tdsops_init
  use m_omptgt_kernels_fieldops, only: copy, axpby, pwmul, field_scale, &
                                       field_shift, buffer_copy, &
                                       vorticity_from_gradients, &
                                       qcriterion_from_gradients, &
                                       smagorinsky_from_gradients, &
                                       sgs_stress_from_gradients, &
                                       scalar_product, scalar_product_cart, &
                                       vector_norm_squared, field_max_sum, &
                                       volume_integral, slice_max_sum, &
                                       field_set_x_face, field_set_y_face, &
                                       field_set_y_plane, &
                                       field_set_abl_wall_stress, &
                                       field_set_x_face_from_field, &
                                       field_set_y_face_from_field, &
                                       field_add_x_face_from_field, &
                                       field_add_y_face_from_field, &
                                       field_plane_partial_sums, &
                                       field_plane_sums_reduce, &
                                       labelled_sums, labelled_axpy, &
                                       reorder, sum_intox

  implicit none

  type, extends(base_backend_t) :: omptgt_backend_t
    !> Device halo buffers, as in cuda_backend_t
    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      u_recv_s_dev => null(), u_recv_e_dev => null(), &
      u_send_s_dev => null(), u_send_e_dev => null(), &
      v_recv_s_dev => null(), v_recv_e_dev => null(), &
      v_send_s_dev => null(), v_send_e_dev => null(), &
      w_recv_s_dev => null(), w_recv_e_dev => null(), &
      w_send_s_dev => null(), w_send_e_dev => null(), &
      du_send_s_dev => null(), du_send_e_dev => null(), &
      du_recv_s_dev => null(), du_recv_e_dev => null(), &
      dud_send_s_dev => null(), dud_send_e_dev => null(), &
      dud_recv_s_dev => null(), dud_recv_e_dev => null(), &
      d2u_send_s_dev => null(), d2u_send_e_dev => null(), &
      d2u_recv_s_dev => null(), d2u_recv_e_dev => null()
  contains
    procedure :: alloc_tdsops => alloc_omptgt_tdsops
    procedure :: transeq_x => transeq_x_omptgt
    procedure :: transeq_y => transeq_y_omptgt
    procedure :: transeq_z => transeq_z_omptgt
    procedure :: transeq_species => transeq_species_omptgt
    procedure :: tds_solve => tds_solve_omptgt
    procedure :: thom_solve => thom_solve_omptgt
    procedure :: reorder => reorder_omptgt
    procedure :: sum_yintox => sum_yintox_omptgt
    procedure :: sum_zintox => sum_zintox_omptgt
    procedure :: veccopy => veccopy_omptgt
    procedure :: vecadd => vecadd_omptgt
    procedure :: vecmult => vecmult_omptgt
    procedure :: scalar_product => scalar_product_omptgt
    procedure :: vector_norm_squared => vector_norm_squared_omptgt
    procedure :: field_labelled_sums => field_labelled_sums_omptgt
    procedure :: field_labelled_axpy => field_labelled_axpy_omptgt
    procedure :: field_max_mean => field_max_mean_omptgt
    procedure :: slice_max_sum => slice_max_sum_omptgt
    procedure :: field_scale => field_scale_omptgt
    procedure :: field_shift => field_shift_omptgt
    procedure :: field_set_face => field_set_face_omptgt
    procedure :: field_set_y_plane => field_set_y_plane_omptgt
    procedure :: field_plane_sums => field_plane_sums_omptgt
    procedure :: field_set_abl_wall_stress => field_set_abl_wall_stress_omptgt
    procedure :: field_set_face_from_field => field_set_face_from_field_omptgt
    procedure :: field_add_face_from_field => field_add_face_from_field_omptgt
    procedure :: compute_vorticity => compute_vorticity_omptgt
    procedure :: compute_qcriterion => compute_qcriterion_omptgt
    procedure :: compute_smagorinsky_nut => compute_smagorinsky_nut_omptgt
    procedure :: compute_sgs_stress => compute_sgs_stress_omptgt
    procedure :: field_volume_integral => field_volume_integral_omptgt
    procedure :: copy_data_to_f => copy_data_to_f_omptgt
    procedure :: copy_f_to_data => copy_f_to_data_omptgt
    procedure :: sync => sync_omptgt
    procedure :: get_device_bw_info => get_device_bw_info_omptgt
    procedure :: transeq_omptgt_dist
    procedure :: tds_solve_dist
  end type omptgt_backend_t

  interface omptgt_backend_t
    module procedure omptgt_backend_init
  end interface

  private
  public :: omptgt_backend_t, resolve_field_t

contains

  type(omptgt_backend_t) function omptgt_backend_init(mesh, allocator) &
    result(backend)
    !! Constructs the backend over a device-resident allocator and allocates
    !! its halo buffers, sized for the largest exchange of any direction.

    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator

    integer :: n_groups, n_halo

    call backend%base_init()

    backend%allocator => allocator
    backend%mesh => mesh

    n_groups = maxval([backend%allocator%get_n_groups(DIR_X), &
                       backend%allocator%get_n_groups(DIR_Y), &
                       backend%allocator%get_n_groups(DIR_Z)])
    n_halo = backend%n_halo

    call dev_alloc(backend%u_send_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%u_send_e_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%u_recv_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%u_recv_e_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%v_send_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%v_send_e_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%v_recv_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%v_recv_e_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%w_send_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%w_send_e_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%w_recv_s_dev, SZ, n_halo, n_groups)
    call dev_alloc(backend%w_recv_e_dev, SZ, n_halo, n_groups)

    call dev_alloc(backend%du_send_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%du_send_e_dev, SZ, 1, n_groups)
    call dev_alloc(backend%du_recv_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%du_recv_e_dev, SZ, 1, n_groups)
    call dev_alloc(backend%dud_send_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%dud_send_e_dev, SZ, 1, n_groups)
    call dev_alloc(backend%dud_recv_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%dud_recv_e_dev, SZ, 1, n_groups)
    call dev_alloc(backend%d2u_send_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%d2u_send_e_dev, SZ, 1, n_groups)
    call dev_alloc(backend%d2u_recv_s_dev, SZ, 1, n_groups)
    call dev_alloc(backend%d2u_recv_e_dev, SZ, 1, n_groups)

  end function omptgt_backend_init

  subroutine alloc_omptgt_tdsops( &
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu, &
    filter_alpha &
    )
    implicit none

    class(omptgt_backend_t) :: self
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
    real(dp), optional, intent(in) :: filter_alpha

    allocate (omptgt_tdsops_t :: tdsops)

    ! Built in place: see m_omptgt_tdsops for why not by assignment
    select type (tdsops)
    type is (omptgt_tdsops_t)
      call omptgt_tdsops_init(tdsops, n_tds, delta, operation, scheme, &
                              bc_start, bc_end, stretch, stretch_correct, &
                              n_halo, from_to, sym, c_nu, nu0_nu, &
                              filter_alpha)
    end select

  end subroutine alloc_omptgt_tdsops

  subroutine sync_omptgt(self)
    !! Waits for outstanding device work (see device_sync).

    class(omptgt_backend_t) :: self

    call device_sync()

  end subroutine sync_omptgt

  subroutine get_device_bw_info_omptgt(self, mem_clock_rt, mem_bus_width, &
                                       available)
    !! Reports the memory clock rate (kHz) and bus width (bits) of the
    !! default device. OpenMP has no query for them, so they come from the
    !! vendor runtime (see m_omptgt_device_info); available is false when
    !! that runtime cannot be found.

    class(omptgt_backend_t) :: self
    integer, intent(out) :: mem_clock_rt
    integer, intent(out) :: mem_bus_width
    logical, intent(out) :: available

    call query_device_memory(omp_get_default_device(), mem_clock_rt, &
                             mem_bus_width, available)

  end subroutine get_device_bw_info_omptgt

  ! --------------------------------------------------------------------------
  ! Transport equation
  ! --------------------------------------------------------------------------

  subroutine transeq_x_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call self%transeq_omptgt_dist(du, dv, dw, u, v, w, nu, dirps)

  end subroutine transeq_x_omptgt

  subroutine transeq_y_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    ! u, v, w is reordered so that we pass v, u, w
    call self%transeq_omptgt_dist(dv, du, dw, v, u, w, nu, dirps)

  end subroutine transeq_y_omptgt

  subroutine transeq_z_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    ! u, v, w is reordered so that we pass w, u, v
    call self%transeq_omptgt_dist(dw, du, dv, w, u, v, nu, dirps)

  end subroutine transeq_z_omptgt

  subroutine transeq_species_omptgt(self, dspec, uvw, spec, nu, dirps, sync)
    !! Computes the convection and diffusion of a scalar species in the
    !! direction of dirps. The species is always halo-exchanged; the
    !! transporting velocity only when sync is true, otherwise the u buffers
    !! still hold its halos from the preceding momentum transport.
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: dspec
    class(field_t), intent(in) :: uvw, spec
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps
    logical, intent(in) :: sync

    integer :: n_groups, n
    type(omptgt_tdsops_t), pointer :: der1st, der1st_sym, der2nd, der2nd_sym
    real(dp), pointer, contiguous, dimension(:, :, :) :: u_dev, spec_dev, &
                                                         dspec_dev

    call resolve_field_t(u_dev, uvw)
    call resolve_field_t(spec_dev, spec)
    call resolve_field_t(dspec_dev, dspec)

    call resolve_dirps(dirps, der1st, der1st_sym, der2nd, der2nd_sym)

    n_groups = self%allocator%get_n_groups(dirps%dir)
    n = self%mesh%get_n(dirps%dir, VERT)

    if (sync) then
      call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, n, &
                             n_groups)
      call sendrecv_fields(self%u_recv_s_dev, self%u_recv_e_dev, &
                           self%u_send_s_dev, self%u_send_e_dev, &
                           SZ*self%n_halo*n_groups, &
                           self%mesh%par%nproc_dir(dirps%dir), &
                           self%mesh%par%pprev(dirps%dir), &
                           self%mesh%par%pnext(dirps%dir))
    end if

    call copy_into_buffers(self%v_send_s_dev, self%v_send_e_dev, spec_dev, &
                           n, n_groups)
    call sendrecv_fields(self%v_recv_s_dev, self%v_recv_e_dev, &
                         self%v_send_s_dev, self%v_send_e_dev, &
                         SZ*self%n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dirps%dir), &
                         self%mesh%par%pprev(dirps%dir), &
                         self%mesh%par%pnext(dirps%dir))

    call transeq_dist_component(self, dspec_dev, spec_dev, u_dev, nu, &
                                self%v_recv_s_dev, self%v_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st, der1st_sym, der2nd, dirps%dir)

    call dspec%set_data_loc(spec%data_loc)

  end subroutine transeq_species_omptgt

  subroutine transeq_omptgt_dist(self, du, dv, dw, u, v, w, nu, dirps)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      u_dev, v_dev, w_dev, du_dev, dv_dev, dw_dev
    type(omptgt_tdsops_t), pointer :: der1st, der1st_sym, der2nd, der2nd_sym

    call resolve_field_t(u_dev, u)
    call resolve_field_t(v_dev, v)
    call resolve_field_t(w_dev, w)

    call resolve_field_t(du_dev, du)
    call resolve_field_t(dv_dev, dv)
    call resolve_field_t(dw_dev, dw)

    call resolve_dirps(dirps, der1st, der1st_sym, der2nd, der2nd_sym)

    call transeq_halo_exchange(self, u_dev, v_dev, w_dev, dirps%dir)

    call transeq_dist_component(self, du_dev, u_dev, u_dev, nu, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st, der1st_sym, der2nd, dirps%dir)
    call transeq_dist_component(self, dv_dev, v_dev, u_dev, nu, &
                                self%v_recv_s_dev, self%v_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st_sym, der1st, der2nd_sym, dirps%dir)
    call transeq_dist_component(self, dw_dev, w_dev, u_dev, nu, &
                                self%w_recv_s_dev, self%w_recv_e_dev, &
                                self%u_recv_s_dev, self%u_recv_e_dev, &
                                der1st_sym, der1st, der2nd_sym, dirps%dir)

    call du%set_data_loc(u%data_loc)
    call dv%set_data_loc(v%data_loc)
    call dw%set_data_loc(w%data_loc)

  end subroutine transeq_omptgt_dist

  subroutine transeq_halo_exchange(self, u_dev, v_dev, w_dev, dir)
    class(omptgt_backend_t) :: self
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      u_dev, v_dev, w_dev
    integer, intent(in) :: dir

    integer :: n, n_groups

    n_groups = self%allocator%get_n_groups(dir)
    n = self%mesh%get_n(dir, VERT)

    call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, n, &
                           n_groups)
    call copy_into_buffers(self%v_send_s_dev, self%v_send_e_dev, v_dev, n, &
                           n_groups)
    call copy_into_buffers(self%w_send_s_dev, self%w_send_e_dev, w_dev, n, &
                           n_groups)

    call sendrecv_3fields( &
      self%u_recv_s_dev, self%u_recv_e_dev, &
      self%v_recv_s_dev, self%v_recv_e_dev, &
      self%w_recv_s_dev, self%w_recv_e_dev, &
      self%u_send_s_dev, self%u_send_e_dev, &
      self%v_send_s_dev, self%v_send_e_dev, &
      self%w_send_s_dev, self%w_send_e_dev, &
      SZ*self%n_halo*n_groups, self%mesh%par%nproc_dir(dir), &
      self%mesh%par%pprev(dir), self%mesh%par%pnext(dir))

  end subroutine transeq_halo_exchange

  subroutine transeq_dist_component(self, rhs_du_dev, u_dev, conv_dev, nu, &
                                    u_recv_s_dev, u_recv_e_dev, &
                                    conv_recv_s_dev, conv_recv_e_dev, &
                                    tdsops_du, tdsops_dud, tdsops_d2u, dir)
    !! Computes RHS_x^u following:
    !!
    !! rhs_x^u = -0.5*(conv*du/dx + d(u*conv)/dx) + nu*d2u/dx2
    class(omptgt_backend_t) :: self
    !> The result field, it is also used as temporary storage
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: rhs_du_dev
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      u_dev, conv_dev
    real(dp), intent(in) :: nu
    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      u_recv_s_dev, u_recv_e_dev, conv_recv_s_dev, conv_recv_e_dev
    type(omptgt_tdsops_t), intent(in) :: tdsops_du, tdsops_dud, tdsops_d2u
    integer, intent(in) :: dir

    class(field_t), pointer :: dud, d2u
    real(dp), pointer, contiguous, dimension(:, :, :) :: dud_dev, d2u_dev

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
      tdsops_du, tdsops_dud, tdsops_d2u, nu, &
      self%mesh%par%nproc_dir(dir), self%mesh%par%pprev(dir), &
      self%mesh%par%pnext(dir), self%allocator%get_n_groups(dir) &
      )

    call self%allocator%release_block(dud)
    call self%allocator%release_block(d2u)

  end subroutine transeq_dist_component

  ! --------------------------------------------------------------------------
  ! Tridiagonal solves
  ! --------------------------------------------------------------------------

  subroutine tds_solve_omptgt(self, du, u, tdsops)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    ! Thomas is exact along an undecomposed line (see tdsops%prefer_thomas)
    if (tdsops%prefer_thomas) then
      call self%thom_solve(du, u, tdsops)
      return
    end if

    if (u%dir /= du%dir) then
      error stop 'DIR mismatch between fields in tds_solve.'
    end if

    if (u%data_loc /= NULL_LOC) then
      call du%set_data_loc(move_data_loc(u%data_loc, u%dir, tdsops%move))
    end if

    call self%tds_solve_dist(du, u, tdsops)

  end subroutine tds_solve_omptgt

  subroutine thom_solve_omptgt(self, du, u, tdsops)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    real(dp), pointer, contiguous, dimension(:, :, :) :: du_dev, u_dev

    if (u%dir /= du%dir) then
      error stop 'DIR mismatch between fields in thom_solve.'
    end if

    if (u%data_loc /= NULL_LOC) then
      call du%set_data_loc(move_data_loc(u%data_loc, u%dir, tdsops%move))
    end if

    call resolve_field_t(du_dev, du)
    call resolve_field_t(u_dev, u)

    select type (tdsops)
    type is (omptgt_tdsops_t)
      call exec_thom_tds_compact(du_dev, u_dev, tdsops, &
                                 self%allocator%get_n_groups(u%dir))
    class default
      error stop 'Expected omptgt_tdsops_t in thom_solve_omptgt.'
    end select

  end subroutine thom_solve_omptgt

  subroutine tds_solve_dist(self, du, u, tdsops)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    real(dp), pointer, contiguous, dimension(:, :, :) :: du_dev, u_dev
    integer :: n_groups, dir

    dir = u%dir
    n_groups = self%allocator%get_n_groups(dir)

    call resolve_field_t(du_dev, du)
    call resolve_field_t(u_dev, u)

    select type (tdsops)
    type is (omptgt_tdsops_t)
      call copy_into_buffers(self%u_send_s_dev, self%u_send_e_dev, u_dev, &
                             tdsops%n_tds, n_groups)

      call sendrecv_fields(self%u_recv_s_dev, self%u_recv_e_dev, &
                           self%u_send_s_dev, self%u_send_e_dev, &
                           SZ*self%n_halo*n_groups, &
                           self%mesh%par%nproc_dir(dir), &
                           self%mesh%par%pprev(dir), &
                           self%mesh%par%pnext(dir))

      call exec_dist_tds_compact( &
        du_dev, u_dev, &
        self%u_recv_s_dev, self%u_recv_e_dev, &
        self%du_send_s_dev, self%du_send_e_dev, &
        self%du_recv_s_dev, self%du_recv_e_dev, &
        tdsops, self%mesh%par%nproc_dir(dir), &
        self%mesh%par%pprev(dir), self%mesh%par%pnext(dir), n_groups &
        )
    class default
      error stop 'Expected omptgt_tdsops_t in tds_solve_dist.'
    end select

  end subroutine tds_solve_dist

  ! --------------------------------------------------------------------------
  ! Reordering
  ! --------------------------------------------------------------------------

  subroutine reorder_omptgt(self, u_, u, direction)
    !! Reorders `u` into `u_` between the two data layouts encoded in
    !! `direction`, in a single pass for every pair of layouts (the CUDA
    !! backend goes through DIR_X for the pairs involving DIR_C).
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, intent(in) :: direction

    real(dp), pointer, contiguous, dimension(:, :, :) :: u_o_d, u_i_d
    integer :: dir_from, dir_to

    call get_dirs_from_rdr(dir_from, dir_to, direction)

    call resolve_field_t(u_o_d, u_)
    call resolve_field_t(u_i_d, u)

    call reorder(u_o_d, u_i_d, dir_from, dir_to, &
                 self%allocator%get_padded_dims(DIR_C))

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omptgt

  subroutine sum_yintox_omptgt(self, u, u_)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call sum_intox_omptgt(self, u, u_, DIR_Y)

  end subroutine sum_yintox_omptgt

  subroutine sum_zintox_omptgt(self, u, u_)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call sum_intox_omptgt(self, u, u_, DIR_Z)

  end subroutine sum_zintox_omptgt

  subroutine sum_intox_omptgt(self, u, u_, dir_from)
    !! u += u_, u a DIR_X field and u_ the same quantity in layout dir_from.
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_
    integer, intent(in) :: dir_from

    real(dp), pointer, contiguous, dimension(:, :, :) :: u_d, u_d_

    call resolve_field_t(u_d, u)
    call resolve_field_t(u_d_, u_)

    call sum_intox(u_d, u_d_, dir_from, self%allocator%get_padded_dims(DIR_C))

  end subroutine sum_intox_omptgt

  ! --------------------------------------------------------------------------
  ! Pointwise operations
  ! --------------------------------------------------------------------------

  subroutine veccopy_omptgt(self, dst, src)
    !! Copies `src` into `dst`. Both fields must share the same direction.
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: dst
    class(field_t), intent(in) :: src

    real(dp), pointer, contiguous, dimension(:, :, :) :: dst_d, src_d

    if (src%dir /= dst%dir) then
      error stop "Called vector copy with incompatible fields"
    end if

    call resolve_field_t(dst_d, dst)
    call resolve_field_t(src_d, src)
    if (size(dst_d) < size(src_d)) then
      error stop "Destination field is smaller than the source"
    end if

    call copy(size(src_d), dst_d, src_d)

  end subroutine veccopy_omptgt

  subroutine vecadd_omptgt(self, a, x, b, y)
    !! Computes y = a*x + b*y.
    class(omptgt_backend_t) :: self
    real(dp), intent(in) :: a
    class(field_t), intent(in) :: x
    real(dp), intent(in) :: b
    class(field_t), intent(inout) :: y

    real(dp), pointer, contiguous, dimension(:, :, :) :: x_d, y_d

    if (x%dir /= y%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    call resolve_field_t(x_d, x)
    call resolve_field_t(y_d, y)

    call axpby(size(x_d), a, x_d, b, y_d)

  end subroutine vecadd_omptgt

  subroutine vecmult_omptgt(self, y, x)
    !! [[m_base_backend(module):vecmult_base(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: y
    class(field_t), intent(in) :: x

    real(dp), pointer, contiguous, dimension(:, :, :) :: x_d, y_d

    call resolve_field_t(x_d, x)
    call resolve_field_t(y_d, y)

    call pwmul(y_d, x_d, size(y_d))

  end subroutine vecmult_omptgt

  subroutine field_scale_omptgt(self, f, a)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d

    call resolve_field_t(f_d, f)
    call field_scale(f_d, a, size(f_d))

  end subroutine field_scale_omptgt

  subroutine field_shift_omptgt(self, f, a)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d

    call resolve_field_t(f_d, f)
    call field_shift(f_d, a, size(f_d))

  end subroutine field_shift_omptgt

  subroutine compute_vorticity_omptgt( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      vort_d, dudy_d, dudz_d, dvdx_d, dvdz_d, dwdx_d, dwdy_d

    call resolve_field_t(vort_d, field_out)
    call resolve_field_t(dudy_d, dudy)
    call resolve_field_t(dudz_d, dudz)
    call resolve_field_t(dvdx_d, dvdx)
    call resolve_field_t(dvdz_d, dvdz)
    call resolve_field_t(dwdx_d, dwdx)
    call resolve_field_t(dwdy_d, dwdy)

    call vorticity_from_gradients(vort_d, dudy_d, dudz_d, dvdx_d, dvdz_d, &
                                  dwdx_d, dwdy_d, size(vort_d))

  end subroutine compute_vorticity_omptgt

  subroutine compute_qcriterion_omptgt( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      qcrit_d, dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, dwdx_d, &
      dwdy_d, dwdz_d

    call resolve_field_t(qcrit_d, field_out)
    call resolve_field_t(dudx_d, dudx)
    call resolve_field_t(dudy_d, dudy)
    call resolve_field_t(dudz_d, dudz)
    call resolve_field_t(dvdx_d, dvdx)
    call resolve_field_t(dvdy_d, dvdy)
    call resolve_field_t(dvdz_d, dvdz)
    call resolve_field_t(dwdx_d, dwdx)
    call resolve_field_t(dwdy_d, dwdy)
    call resolve_field_t(dwdz_d, dwdz)

    call qcriterion_from_gradients(qcrit_d, dudx_d, dudy_d, dudz_d, dvdx_d, &
                                   dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d, &
                                   size(qcrit_d))

  end subroutine compute_qcriterion_omptgt

  subroutine compute_smagorinsky_nut_omptgt( &
    self, nut, mixing_length_sq, dudx, dudy, dudz, dvdx, dvdy, dvdz, &
    dwdx, dwdy, dwdz)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: nut
    class(field_t), intent(in) :: mixing_length_sq
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      nut_d, mixing_length_sq_d, dudx_d, dudy_d, dudz_d, &
      dvdx_d, dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d

    call resolve_field_t(nut_d, nut)
    call resolve_field_t(mixing_length_sq_d, mixing_length_sq)
    call resolve_field_t(dudx_d, dudx)
    call resolve_field_t(dudy_d, dudy)
    call resolve_field_t(dudz_d, dudz)
    call resolve_field_t(dvdx_d, dvdx)
    call resolve_field_t(dvdy_d, dvdy)
    call resolve_field_t(dvdz_d, dvdz)
    call resolve_field_t(dwdx_d, dwdx)
    call resolve_field_t(dwdy_d, dwdy)
    call resolve_field_t(dwdz_d, dwdz)

    call smagorinsky_from_gradients(nut_d, mixing_length_sq_d, dudx_d, &
                                    dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, &
                                    dwdx_d, dwdy_d, dwdz_d, size(nut_d))

  end subroutine compute_smagorinsky_nut_omptgt

  subroutine compute_sgs_stress_omptgt( &
    self, stress, nut, gradient_a, gradient_b, scale_a, scale_b)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: stress
    class(field_t), intent(in) :: nut, gradient_a, gradient_b
    real(dp), intent(in) :: scale_a, scale_b

    real(dp), pointer, contiguous, dimension(:, :, :) :: &
      stress_d, nut_d, gradient_a_d, gradient_b_d

    call resolve_field_t(stress_d, stress)
    call resolve_field_t(nut_d, nut)
    call resolve_field_t(gradient_a_d, gradient_a)
    call resolve_field_t(gradient_b_d, gradient_b)

    call sgs_stress_from_gradients(stress_d, nut_d, gradient_a_d, &
                                   gradient_b_d, scale_a, scale_b, &
                                   size(stress_d))

  end subroutine compute_sgs_stress_omptgt

  ! --------------------------------------------------------------------------
  ! Reductions
  ! --------------------------------------------------------------------------

  subroutine pencil_extents(self, dir, data_loc, n, n_i, n_groups_pad, n_j)
    !! Line length `n`, stacked extent `n_j` and repeated extent `n_i` of a
    !! pencil field in layout `dir`, with the number of groups the padded
    !! stacked extent makes (see get_index_dir in m_ordering).
    class(omptgt_backend_t) :: self
    integer, intent(in) :: dir, data_loc
    integer, intent(out) :: n, n_i, n_groups_pad, n_j

    integer :: dims(3), dims_padded(3)

    dims = self%mesh%get_dims(data_loc)
    dims_padded = self%allocator%get_padded_dims(DIR_C)

    if (dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3)
      n_groups_pad = dims_padded(2)/SZ
    else if (dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3)
      n_groups_pad = dims_padded(1)/SZ
    else if (dir == DIR_Z) then
      n = dims(3); n_j = dims(1); n_i = dims(2)
      n_groups_pad = dims_padded(1)/SZ
    else
      error stop 'pencil_extents: DIR_C has no pencil layout.'
    end if

  end subroutine pencil_extents

  real(dp) function scalar_product_omptgt(self, x, y) result(s)
    !! [[m_base_backend(module):scalar_product_base(function)]]
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: x, y

    real(dp), pointer, contiguous, dimension(:, :, :) :: x_d, y_d
    integer :: n, n_i, n_groups_pad, n_j, ierr

    if ((x%data_loc == NULL_LOC) .or. (y%data_loc == NULL_LOC)) then
      error stop "You must set the data_loc before calling scalar product"
    end if
    if ((x%data_loc /= y%data_loc) .or. (x%dir /= y%dir)) then
      error stop "Called scalar product with incompatible fields"
    end if

    call resolve_field_t(x_d, x)
    call resolve_field_t(y_d, y)

    if (x%dir == DIR_C) then
      s = scalar_product_cart(x_d, y_d, shape(x_d), &
                              self%mesh%get_dims(x%data_loc))
    else
      call pencil_extents(self, x%dir, x%data_loc, n, n_i, n_groups_pad, n_j)
      s = scalar_product(x_d, y_d, size(x_d, 2), n, n_i, n_groups_pad, n_j)
    end if

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function scalar_product_omptgt

  real(dp) function vector_norm_squared_omptgt(self, a, b, c) &
    result(norm_squared)
    !! Global sum of a**2 + b**2 + c**2, with one MPI reduction.
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: a, b, c

    real(dp), pointer, contiguous, dimension(:, :, :) :: a_d, b_d, c_d
    real(dp) :: local_sum
    integer :: n, n_i, n_groups_pad, n_j, ierr

    if (a%data_loc == NULL_LOC .or. b%data_loc == NULL_LOC .or. &
        c%data_loc == NULL_LOC) then
      error stop 'You must set data_loc before computing a vector norm.'
    end if
    if (a%data_loc /= b%data_loc .or. a%data_loc /= c%data_loc) then
      error stop 'Vector-norm fields must use the same data location.'
    end if
    if (a%dir /= DIR_X .or. b%dir /= DIR_X .or. c%dir /= DIR_X) then
      error stop 'Vector-norm fields must use DIR_X layout.'
    end if

    call resolve_field_t(a_d, a)
    call resolve_field_t(b_d, b)
    call resolve_field_t(c_d, c)

    call pencil_extents(self, DIR_X, a%data_loc, n, n_i, n_groups_pad, n_j)
    local_sum = vector_norm_squared(a_d, b_d, c_d, size(a_d, 2), n, n_i, &
                                    n_groups_pad, n_j)

    call MPI_Allreduce(local_sum, norm_squared, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function vector_norm_squared_omptgt

  subroutine field_labelled_sums_omptgt(self, sums, g, label, a, b, c)
    !! [[m_base_backend(module):field_labelled_sums_base(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: sums(:, :)
    class(field_t), intent(in) :: g, label, a, b, c

    real(dp), pointer, contiguous, dimension(:, :, :) :: g_d, label_d, a_d, &
                                                         b_d, c_d
    real(dp), pointer, contiguous :: sums_d(:)
    real(dp), allocatable :: sums_flat(:)
    integer :: n, n_i, n_groups_pad, n_j

    call check_labelled_fields(g, label, a, b, c)

    call resolve_field_t(g_d, g)
    call resolve_field_t(label_d, label)
    call resolve_field_t(a_d, a)
    call resolve_field_t(b_d, b)
    call resolve_field_t(c_d, c)

    ! sums is passed flat: the device helpers have no 2D download
    allocate (sums_flat(size(sums)))
    sums_flat = 0._dp
    call dev_alloc(sums_d, size(sums))
    call dev_upload(sums_d, sums_flat)

    call pencil_extents(self, DIR_X, g%data_loc, n, n_i, n_groups_pad, n_j)
    call labelled_sums(sums_d, g_d, label_d, a_d, b_d, c_d, size(sums, 2), &
                       size(g_d, 2), n, n_i, n_groups_pad, n_j)

    call dev_download(sums_flat, sums_d)
    call dev_free(sums_d)
    sums = reshape(sums_flat, shape(sums))

  end subroutine field_labelled_sums_omptgt

  subroutine field_labelled_axpy_omptgt(self, coeffs, g, label, a, b, c)
    !! [[m_base_backend(module):field_labelled_axpy_base(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    real(dp), intent(in) :: coeffs(:, :)
    class(field_t), intent(in) :: g, label
    class(field_t), intent(inout) :: a, b, c

    real(dp), pointer, contiguous, dimension(:, :, :) :: g_d, label_d, a_d, &
                                                         b_d, c_d
    real(dp), pointer, contiguous :: coeffs_d(:, :)
    integer :: n, n_i, n_groups_pad, n_j

    call check_labelled_fields(g, label, a, b, c)

    call resolve_field_t(g_d, g)
    call resolve_field_t(label_d, label)
    call resolve_field_t(a_d, a)
    call resolve_field_t(b_d, b)
    call resolve_field_t(c_d, c)

    call dev_alloc(coeffs_d, size(coeffs, 1), size(coeffs, 2))
    call dev_upload(coeffs_d, coeffs)

    call pencil_extents(self, DIR_X, g%data_loc, n, n_i, n_groups_pad, n_j)
    call labelled_axpy(coeffs_d, g_d, label_d, a_d, b_d, c_d, &
                       size(coeffs, 2), size(g_d, 2), n, n_i, n_groups_pad, &
                       n_j)

    ! The kernel may still be reading coeffs_d (see device_sync)
    call device_sync()
    call dev_free(coeffs_d)

  end subroutine field_labelled_axpy_omptgt

  subroutine check_labelled_fields(g, label, a, b, c)
    !! Checks the fields passed to the labelled field operations share a
    !! valid data location and the DIR_X layout.
    implicit none

    class(field_t), intent(in) :: g, label, a, b, c

    if (g%data_loc == NULL_LOC .or. label%data_loc == NULL_LOC .or. &
        a%data_loc == NULL_LOC .or. b%data_loc == NULL_LOC .or. &
        c%data_loc == NULL_LOC) then
      error stop 'You must set data_loc before labelled field operations.'
    end if
    if (label%data_loc /= g%data_loc .or. a%data_loc /= g%data_loc .or. &
        b%data_loc /= g%data_loc .or. c%data_loc /= g%data_loc) then
      error stop 'Labelled field operations need one data location.'
    end if
    if (g%dir /= DIR_X .or. label%dir /= DIR_X .or. a%dir /= DIR_X .or. &
        b%dir /= DIR_X .or. c%dir /= DIR_X) then
      error stop 'Labelled field operations need DIR_X layout.'
    end if
  end subroutine check_labelled_fields

  subroutine field_max_mean_omptgt(self, max_val, mean_val, f, &
                                   enforced_data_loc)
    !! [[m_base_backend(module):field_max_mean_base(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: max_val, mean_val
    class(field_t), intent(in) :: f
    integer, optional, intent(in) :: enforced_data_loc

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    integer :: data_loc, n, n_i, n_groups_pad, n_j, ierr

    if (f%data_loc == NULL_LOC .and. (.not. present(enforced_data_loc))) then
      error stop 'The input field to omptgt::field_max_mean does not have a &
                  &valid f%data_loc. You may enforce a data_loc of your &
                  &choice as last argument to carry on at your own risk!'
    end if

    if (present(enforced_data_loc)) then
      data_loc = enforced_data_loc
    else
      data_loc = f%data_loc
    end if

    if (f%dir == DIR_C) then
      error stop 'field_max_mean does not support DIR_C fields!'
    end if
    call pencil_extents(self, f%dir, data_loc, n, n_i, n_groups_pad, n_j)

    call resolve_field_t(f_d, f)
    call field_max_sum(max_val, mean_val, f_d, size(f_d, 2), n, n_i, &
                       n_groups_pad, n_j)
    mean_val = mean_val/product(self%mesh%get_global_dims(data_loc))

    call MPI_Allreduce(MPI_IN_PLACE, max_val, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mean_val, 1, MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine field_max_mean_omptgt

  subroutine slice_max_sum_omptgt(self, max_val, sum_val, f, i_slice, &
                                  enforced_data_loc)
    !! [[m_base_backend(module):slice_max_sum_base(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: max_val, sum_val
    class(field_t), intent(in) :: f
    integer, intent(in) :: i_slice
    integer, optional, intent(in) :: enforced_data_loc

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    integer :: data_loc, n, n_i, n_groups_pad, n_j

    if (f%data_loc == NULL_LOC .and. (.not. present(enforced_data_loc))) then
      error stop 'The input field to omptgt::slice_max_sum does not have a &
                  &valid f%data_loc. You may enforce a data_loc of your &
                  &choice as last argument to carry on at your own risk!'
    end if

    if (present(enforced_data_loc)) then
      data_loc = enforced_data_loc
    else
      data_loc = f%data_loc
    end if

    if (f%dir == DIR_C) then
      error stop 'slice_max_sum does not support DIR_C fields!'
    end if
    call pencil_extents(self, f%dir, data_loc, n, n_i, n_groups_pad, n_j)

    if (i_slice < 1 .or. i_slice > n) then
      error stop 'slice_max_sum: i_slice out of range'
    end if

    ! Rank-local values; caller is responsible for MPI_Allreduce.
    call resolve_field_t(f_d, f)
    call slice_max_sum(max_val, sum_val, f_d, size(f_d, 2), i_slice, n_i, &
                       n_groups_pad, n_j)

  end subroutine slice_max_sum_omptgt

  real(dp) function field_volume_integral_omptgt(self, f) result(s)
    !! volume integral of a field
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    integer :: n, n_i, n_groups_pad, n_j, ierr

    if (f%data_loc == NULL_LOC) then
      error stop 'You must set the data_loc before calling volume integral.'
    end if
    if (f%dir /= DIR_X) then
      error stop 'Volume integral can only be called on DIR_X fields.'
    end if

    call resolve_field_t(f_d, f)
    call pencil_extents(self, DIR_X, f%data_loc, n, n_i, n_groups_pad, n_j)
    s = volume_integral(f_d, size(f_d, 2), n, n_i, n_groups_pad, n_j)

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function field_volume_integral_omptgt

  subroutine field_plane_sums_omptgt(self, sums, f)
    !! [[m_base_backend(module):field_plane_sums(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: sums(:)
    class(field_t), intent(in) :: f

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    real(dp), pointer, contiguous :: partial_d(:, :), sums_d(:)
    integer :: dims(3)

    if (f%dir /= DIR_X) &
      error stop 'field_plane_sums is only supported for DIR_X fields.'
    if (f%data_loc == NULL_LOC) &
      error stop 'field_plane_sums requires a valid data_loc.'

    dims = self%mesh%get_dims(f%data_loc)
    if (size(sums) < dims(2)) &
      error stop 'field_plane_sums: sums is smaller than the y extent.'

    call resolve_field_t(f_d, f)
    call dev_alloc(partial_d, dims(2), dims(3))
    call dev_alloc(sums_d, dims(2))

    ! Sum along x for every (y, z) line, then over z for every y, each in a
    ! fixed order so the result is reproducible.
    call field_plane_partial_sums(partial_d, f_d, size(f_d, 2), dims(1), &
                                  dims(2), dims(3), n_yb_pad(self))
    call field_plane_sums_reduce(sums_d, partial_d, dims(2), dims(3))

    call dev_download(sums(1:dims(2)), sums_d)
    call dev_free(partial_d)
    call dev_free(sums_d)

  end subroutine field_plane_sums_omptgt

  ! --------------------------------------------------------------------------
  ! Faces and planes
  ! --------------------------------------------------------------------------

  subroutine field_set_face_omptgt(self, f, c_start, c_end, face, &
                                   bc_start, bc_end, flow_rate_diff)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c_start, c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start, bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    integer :: dims(3)
    integer :: bc_s, bc_e
    real(dp) :: flow_rate_diff_val

    if (f%dir /= DIR_X) then
      error stop 'Setting a field face is only supported for DIR_X fields.'
    end if
    if (f%data_loc == NULL_LOC) then
      error stop 'field_set_face requires a valid data_loc.'
    end if

    bc_s = BC_DIRICHLET
    bc_e = BC_DIRICHLET
    flow_rate_diff_val = 0._dp
    if (present(bc_start)) bc_s = bc_start
    if (present(bc_end)) bc_e = bc_end
    if (present(flow_rate_diff)) flow_rate_diff_val = flow_rate_diff

    call resolve_field_t(f_d, f)
    dims = self%mesh%get_dims(f%data_loc)

    select case (face)
    case (X_FACE)
      call field_set_x_face(f_d, c_start, c_end, bc_s, bc_e, &
                            flow_rate_diff_val, size(f_d, 2), &
                            dims(1), dims(2), dims(3), n_yb_pad(self))
    case (Y_FACE)
      ! Only Dirichlet faces carry a prescribed value; a non-Dirichlet face
      ! is left untouched so the two y-faces can differ (ABL: no-slip floor,
      ! free-slip lid).
      call field_set_y_face(f_d, c_start, c_end, bc_s == BC_DIRICHLET, &
                            bc_e == BC_DIRICHLET, size(f_d, 2), &
                            dims(1), dims(2), dims(3), n_yb_pad(self))
    case (Z_FACE)
      error stop 'Setting Z_FACE is not yet supported.'
    case default
      error stop 'face is undefined.'
    end select

  end subroutine field_set_face_omptgt

  subroutine field_set_y_plane_omptgt(self, f, c, plane)
    !! [[m_base_backend(module):field_set_y_plane(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c
    integer, intent(in) :: plane

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d
    integer :: dims(3)

    if (f%dir /= DIR_X) &
      error stop 'field_set_y_plane is only supported for DIR_X fields.'
    if (f%data_loc == NULL_LOC) &
      error stop 'field_set_y_plane requires a valid data_loc.'

    call resolve_field_t(f_d, f)
    dims = self%mesh%get_dims(f%data_loc)
    if (plane < 1 .or. plane > dims(2)) &
      error stop 'field_set_y_plane: plane outside the domain.'

    call field_set_y_plane(f_d, c, plane, size(f_d, 2), dims(1), dims(3), &
                           n_yb_pad(self))

  end subroutine field_set_y_plane_omptgt

  subroutine field_set_abl_wall_stress_omptgt( &
    self, stress, u, w, sample_plane, stress_plane, drag_coeff, component)
    !! [[m_base_backend(module):field_set_abl_wall_stress(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: stress
    class(field_t), intent(in) :: u, w
    integer, intent(in) :: sample_plane, stress_plane, component
    real(dp), intent(in) :: drag_coeff

    real(dp), pointer, contiguous, dimension(:, :, :) :: stress_d, u_d, w_d
    integer :: dims(3)

    if (stress%dir /= DIR_X .or. u%dir /= DIR_X .or. w%dir /= DIR_X) &
      error stop 'ABL wall stress requires DIR_X fields.'
    if (stress%data_loc /= VERT .or. u%data_loc /= VERT .or. &
        w%data_loc /= VERT) &
      error stop 'ABL wall stress requires vertex fields.'
    dims = self%mesh%get_dims(VERT)
    if (min(sample_plane, stress_plane) < 1 .or. &
        max(sample_plane, stress_plane) > dims(2)) &
      error stop 'ABL wall stress plane is outside the domain.'
    if (component /= 1 .and. component /= 3) &
      error stop 'Invalid ABL wall stress component.'

    call resolve_field_t(stress_d, stress)
    call resolve_field_t(u_d, u)
    call resolve_field_t(w_d, w)

    call field_set_abl_wall_stress(stress_d, u_d, w_d, drag_coeff, &
                                   sample_plane, stress_plane, &
                                   size(stress_d, 2), dims(1), dims(3), &
                                   n_yb_pad(self), component)

  end subroutine field_set_abl_wall_stress_omptgt

  subroutine field_set_face_from_field_omptgt(self, f, f_start, c_end, face, &
                                              bc_start, bc_end, &
                                              flow_rate_diff)
    !! Set a face of `f` using values supplied by another field `f_start`.
    !! Both fields must be DIR_X VERT pencil-layout. X_FACE uses the inlet
    !! plane (pencil index 1) of f_start plus a convective outflow update
    !! at the right face. Y_FACE copies the bottom and top y-pencil planes
    !! of f_start onto the corresponding planes of f (Dirichlet only).
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: f_start
    real(dp), intent(in) :: c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start, bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d, f_start_d
    integer :: dims(3)
    integer :: bc_s, bc_e
    real(dp) :: flow_rate_diff_val

    if (f%dir /= DIR_X) &
      error stop 'field_set_face_from_field: only supported for DIR_X fields.'
    if (f_start%dir /= DIR_X) &
      error stop 'field_set_face_from_field: f_start must be DIR_X.'

    bc_s = BC_DIRICHLET
    bc_e = BC_DIRICHLET
    flow_rate_diff_val = 0._dp
    if (present(bc_start)) bc_s = bc_start
    if (present(bc_end)) bc_e = bc_end
    if (present(flow_rate_diff)) flow_rate_diff_val = flow_rate_diff

    call resolve_field_t(f_d, f)
    call resolve_field_t(f_start_d, f_start)

    ! The BC fields are always built as VERT in the case setup
    dims = self%mesh%get_dims(VERT)

    select case (face)
    case (X_FACE)
      call field_set_x_face_from_field(f_d, f_start_d, c_end, bc_s, bc_e, &
                                       flow_rate_diff_val, size(f_d, 2), &
                                       dims(1), dims(2), dims(3), &
                                       n_yb_pad(self))
    case (Y_FACE)
      if (bc_s /= BC_DIRICHLET .or. bc_e /= BC_DIRICHLET) then
        error stop &
          'field_set_face_from_field: Y_FACE only supports BC_DIRICHLET.'
      end if
      call field_set_y_face_from_field(f_d, f_start_d, size(f_d, 2), &
                                       dims(1), dims(2), dims(3), &
                                       n_yb_pad(self))
    case (Z_FACE)
      error stop 'field_set_face_from_field: Z_FACE is not yet supported.'
    case default
      error stop 'field_set_face_from_field: face is undefined.'
    end select

  end subroutine field_set_face_from_field_omptgt

  subroutine field_add_face_from_field_omptgt(self, f, g, face, bc_start, &
                                              bc_end)
    !! [[m_base_backend(module):field_add_face_from_field(subroutine)]]
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: g
    integer, intent(in) :: face, bc_start, bc_end

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d, g_d
    integer :: dims(3)
    logical :: set_start, set_end

    if (f%dir /= DIR_X .or. g%dir /= DIR_X) &
      error stop 'field_add_face_from_field: only supported for DIR_X fields.'
    if (f%data_loc == NULL_LOC) &
      error stop 'field_add_face_from_field: requires a valid data_loc.'

    set_start = (bc_start == BC_DIRICHLET)
    set_end = (bc_end == BC_DIRICHLET)

    call resolve_field_t(f_d, f)
    call resolve_field_t(g_d, g)
    dims = self%mesh%get_dims(f%data_loc)

    select case (face)
    case (X_FACE)
      call field_add_x_face_from_field(f_d, g_d, set_start, set_end, &
                                       size(f_d, 2), dims(1), dims(2), &
                                       dims(3), n_yb_pad(self))
    case (Y_FACE)
      call field_add_y_face_from_field(f_d, g_d, set_start, set_end, &
                                       size(f_d, 2), dims(1), dims(2), &
                                       dims(3), n_yb_pad(self))
    case default
      error stop 'field_add_face_from_field: only X_FACE and Y_FACE supported.'
    end select

  end subroutine field_add_face_from_field_omptgt

  ! --------------------------------------------------------------------------
  ! Host transfers
  ! --------------------------------------------------------------------------

  subroutine copy_data_to_f_omptgt(self, f, data)
    !! Copies a host array of the field's shape into a device-resident field.
    class(omptgt_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d

    call resolve_field_t(f_d, f)
    call dev_upload(f_d, data)

  end subroutine copy_data_to_f_omptgt

  subroutine copy_f_to_data_omptgt(self, data, f)
    !! Copies a device-resident field into a host array of its shape.
    class(omptgt_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    real(dp), pointer, contiguous, dimension(:, :, :) :: f_d

    call resolve_field_t(f_d, f)
    call dev_download(data, f_d)

  end subroutine copy_f_to_data_omptgt

  ! --------------------------------------------------------------------------
  ! Helpers
  ! --------------------------------------------------------------------------

  integer function n_yb_pad(self)
    !! Number of y-blocks per z station of a DIR_X field: the stride between
    !! z stations in its group index (see get_index_dir in m_ordering).
    class(omptgt_backend_t) :: self

    integer :: dims_padded(3)

    dims_padded = self%allocator%get_padded_dims(DIR_C)
    n_yb_pad = dims_padded(2)/SZ

  end function n_yb_pad

  subroutine copy_into_buffers(u_send_s_dev, u_send_e_dev, u_dev, n, &
                               n_groups)
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      u_send_s_dev, u_send_e_dev, u_dev
    integer, intent(in) :: n, n_groups

    call buffer_copy(u_send_s_dev, u_send_e_dev, u_dev, size(u_dev, 2), n, &
                     n_groups)

  end subroutine copy_into_buffers

  subroutine resolve_dirps(dirps, der1st, der1st_sym, der2nd, der2nd_sym)
    !! Points at the omptgt operators the transport equation needs.
    type(dirps_t), target, intent(in) :: dirps
    type(omptgt_tdsops_t), pointer, intent(out) :: der1st, der1st_sym, &
                                                   der2nd, der2nd_sym

    select type (tdsops => dirps%der1st)
    type is (omptgt_tdsops_t); der1st => tdsops
    class default; error stop 'Expected omptgt_tdsops_t for der1st.'
    end select
    select type (tdsops => dirps%der1st_sym)
    type is (omptgt_tdsops_t); der1st_sym => tdsops
    class default; error stop 'Expected omptgt_tdsops_t for der1st_sym.'
    end select
    select type (tdsops => dirps%der2nd)
    type is (omptgt_tdsops_t); der2nd => tdsops
    class default; error stop 'Expected omptgt_tdsops_t for der2nd.'
    end select
    select type (tdsops => dirps%der2nd_sym)
    type is (omptgt_tdsops_t); der2nd_sym => tdsops
    class default; error stop 'Expected omptgt_tdsops_t for der2nd_sym.'
    end select

  end subroutine resolve_dirps

  subroutine resolve_field_t(u_dev, u)
    !! Views the device storage of `u` as an array of its current shape. The
    !! pointer holds a device address: pass it to a kernel, never index it on
    !! the host.
    !!
    !! The storage is reached through the type-bound get_dev_ptr, which
    !! omptgt_field_t overrides, rather than a select type: CCE was seen to
    !! take the class default branch of a select type here for a field that
    !! was an omptgt_field_t (perf_transeq, whose fields are host associated
    !! into an internal procedure), while dispatching on the same field
    !! worked.
    use iso_c_binding, only: c_ptr, c_associated, c_f_pointer

    real(dp), pointer, contiguous, dimension(:, :, :), intent(out) :: u_dev
    class(field_t), intent(in) :: u

    type(c_ptr) :: ptr

    ptr = u%get_dev_ptr()
    if (.not. c_associated(ptr)) then
      print '(a,i0,a,i0)', 'omptgt backend: field id ', u%id, ', dir ', &
        u%dir
      error stop 'omptgt backend: field has no device storage.'
    end if
    call c_f_pointer(ptr, u_dev, u%get_shape())

  end subroutine resolve_field_t

end module m_omptgt_backend
