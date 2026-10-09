module m_base_backend
  use iso_c_binding, only: c_ptr

  use m_allocator, only: allocator_t
  use m_common, only: dp, DIR_C, get_rdr_from_dirs
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: tdsops_t, dirps_t

  implicit none

  type, abstract :: base_backend_t
      !! base_backend class defines all the operations that the solver
      !! class requires.
      !!
      !! For example, transport equation in solver class evaluates the
      !! derivatives in x, y, and z directions, and reorders the input
      !! fields as required. Then finally, combines all the directional
      !! derivatives to obtain the divergence of U*.
      !!
      !! All these high level operations solver class executes are
      !! defined here. Every backend implementation extends the present
      !! backend class to define the specifics of these operations based
      !! on the target architecture.
      !!
      !! None of the operations is deferred: each one defaults to
      !! `not_implemented`, which stops the run naming the operation. A
      !! backend therefore overrides only what it actually implements, and
      !! an operation added for one backend does not force a stub into
      !! every other backend. A backend reaching an operation it does not
      !! override fails loudly at the call site instead of silently doing
      !! nothing.
      !!
      !! `supports_device_field_export` has a genuinely correct answer for
      !! a backend that does not override it (`.false.`), so it defaults to
      !! that rather than to an error.

    !> DistD2 implementation is hardcoded for 4 halo layers for all backends
    integer :: n_halo = 4
    type(mesh_t), pointer :: mesh
    class(allocator_t), pointer :: allocator
    class(poisson_fft_t), pointer :: poisson_fft
  contains
    procedure :: transeq_x => transeq_x_base
    procedure :: transeq_y => transeq_y_base
    procedure :: transeq_z => transeq_z_base
    procedure :: transeq_species => transeq_species_base
    procedure :: tds_solve => tds_solve_base
    procedure :: thom_solve => thom_solve_base
    procedure :: reorder => reorder_base
    procedure :: sum_yintox => sum_yintox_base
    procedure :: sum_zintox => sum_zintox_base
    procedure :: veccopy => veccopy_base
    procedure :: vecadd => vecadd_base
    procedure :: vecmult => vecmult_base
    procedure :: scalar_product => scalar_product_base
    procedure :: vector_norm_squared => vector_norm_squared_base
    procedure :: field_max_mean => field_max_mean_base
    procedure :: slice_max_sum => slice_max_sum_base
    procedure :: field_plane_sums => field_plane_sums_base
    procedure :: field_scale => field_scale_base
    procedure :: field_shift => field_shift_base
    procedure :: field_volume_integral => field_volume_integral_base
    procedure :: field_set_face => field_set_face_base
    procedure :: field_set_y_plane => field_set_y_plane_base
    procedure :: field_set_abl_wall_stress => field_set_abl_wall_stress_base
    procedure :: field_set_face_from_field => field_set_face_from_field_base
    procedure :: field_add_face_from_field => field_add_face_from_field_base
    procedure :: compute_vorticity => compute_vorticity_base
    procedure :: compute_qcriterion => compute_qcriterion_base
    procedure :: compute_smagorinsky_nut => compute_smagorinsky_nut_base
    procedure :: compute_sgs_stress => compute_sgs_stress_base
    procedure :: copy_data_to_f => copy_data_to_f_base
    procedure :: copy_f_to_data => copy_f_to_data_base
    procedure :: alloc_tdsops => alloc_tdsops_base
    procedure :: init_poisson_fft => init_poisson_fft_base
    procedure :: sync => sync_base
    procedure :: get_device_bw_info => get_device_bw_info_base
    procedure :: base_init
    procedure :: get_field_data
    procedure :: set_field_data
    procedure :: supports_device_field_export
    procedure :: export_field_to_device
  end type base_backend_t

  private :: not_implemented

contains

  subroutine not_implemented(operation)
    !! Stops the run, naming the backend operation that has no
    !! implementation. The name is printed rather than passed as the stop
    !! code because not every compiler here accepts a variable stop code.
    implicit none

    character(*), intent(in) :: operation

    print *, trim(operation)//': not implemented by the active backend'
    error stop 'Backend operation not implemented'

  end subroutine not_implemented

  subroutine base_init(self)
    implicit none

    class(base_backend_t) :: self

  end subroutine base_init

  subroutine get_field_data(self, data, f, dir)
   !! Extract data from field `f` optionally reordering into `dir` orientation.
   !! To output in same orientation as `f`, use `call ...%get_field_data(data, f, f%dir)`
    implicit none

    class(base_backend_t) :: self
    real(dp), dimension(:, :, :), intent(out) :: data !! Output array
    class(field_t), intent(in) :: f !! Field
    integer, optional, intent(in) :: dir !! Desired orientation of output array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(f%dir, direction)

    ! Carry out a reorder if we need, and copy from field to data array
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction)
      call self%reorder(f_temp, f, rdr_dir)
      call self%copy_f_to_data(data, f_temp)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_f_to_data(data, f)
    end if

  end subroutine get_field_data

  subroutine set_field_data(self, f, data, dir)
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f !! Field
    real(dp), dimension(:, :, :), intent(in) :: data !! Input array
    integer, optional, intent(in) :: dir !! Orientation of input array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(direction, f%dir)

    ! Carry out a reorder if we need, and copy from data array to field
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction, f%data_loc)
      call self%copy_data_to_f(f_temp, data)
      call self%reorder(f, f_temp, rdr_dir)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_data_to_f(f, data)
    end if

  end subroutine set_field_data

  logical function supports_device_field_export(self, f)
    !! Whether `export_field_to_device` can pack field `f` without staging
    !! it through host memory. Backends without device memory keep this
    !! default and report no support.
    implicit none

    class(base_backend_t), intent(in) :: self
    class(field_t), intent(in) :: f

    supports_device_field_export = .false.
  end function supports_device_field_export

  subroutine export_field_to_device(self, buffer, f, dims, to_sp)
    !! Pack the first `dims` Cartesian points of field `f` into a
    !! contiguous device buffer, in (i, j, k) order without padding.
    !!
    !! `buffer` is the device address of product(dims) elements, of kind
    !! `sp` when `to_sp` is true and `dp` otherwise. The buffer holds the
    !! complete result when this returns, so a library that reads device
    !! memory outside the backend's stream can consume it directly.
    !! Only called when `supports_device_field_export(f)` is true.
    implicit none

    class(base_backend_t) :: self
    type(c_ptr), intent(in) :: buffer
    class(field_t), intent(in) :: f
    integer, intent(in) :: dims(3)
    logical, intent(in) :: to_sp

    call not_implemented('export_field_to_device')
  end subroutine export_field_to_device

  subroutine sync_base(self)
    !! Waits for outstanding device work.
    implicit none

    class(base_backend_t) :: self

    call not_implemented('sync')
  end subroutine sync_base

  subroutine get_device_bw_info_base(self, mem_clock_rt, mem_bus_width, &
                                     available)
    !! Reports the memory clock rate and bus width of the target device,
    !! used to work out a theoretical peak bandwidth.
    implicit none

    class(base_backend_t) :: self
    integer, intent(out) :: mem_clock_rt
    integer, intent(out) :: mem_bus_width
    logical, intent(out) :: available

    call not_implemented('get_device_bw_info')
  end subroutine get_device_bw_info_base

  subroutine transeq_x_base(self, du, dv, dw, u, v, w, nu, dirps)
       !! transeq equation obtains the derivatives direction by
       !! direction, and the exact algorithm used to obtain these
       !! derivatives are decided at runtime. Backend implementations
       !! are responsible from directing calls to transeq_x into
       !! the correct algorithm.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_x')

  end subroutine transeq_x_base

  subroutine transeq_y_base(self, du, dv, dw, u, v, w, nu, dirps)
       !! As transeq_x, for the y direction.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_y')

  end subroutine transeq_y_base

  subroutine transeq_z_base(self, du, dv, dw, u, v, w, nu, dirps)
       !! As transeq_x, for the z direction.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_z')

  end subroutine transeq_z_base

  subroutine transeq_species_base(self, dspec, uvw, spec, nu, dirps, sync)
       !! As transeq_x, for a scalar species transported by uvw.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: dspec
    class(field_t), intent(in) :: uvw, spec
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps
    logical, intent(in) :: sync

    call not_implemented('transeq_species')

  end subroutine transeq_species_base

  subroutine tds_solve_base(self, du, u, tdsops)
    !! Applies a tridiagonal operator to u, the exact algorithm used
    !! being decided at runtime. Backend implementations are responsible
    !! from directing calls to tds_solve to the correct algorithm.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    call not_implemented('tds_solve')

  end subroutine tds_solve_base

  subroutine thom_solve_base(self, du, u, tdsops)
    !! As tds_solve, solving the tridiagonal system with the Thomas
    !! algorithm.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    call not_implemented('thom_solve')

  end subroutine thom_solve_base

  subroutine reorder_base(self, u_, u, direction)
       !! reorder subroutines are straightforward, they rearrange
       !! data into our specialist data structure so that regardless
       !! of the direction tridiagonal systems are solved efficiently
       !! and fast.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, intent(in) :: direction

    call not_implemented('reorder')

  end subroutine reorder_base

  subroutine sum_yintox_base(self, u, u_)
       !! sum_yintox combines a y directional field into the
       !! corresponding x directional field.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call not_implemented('sum_yintox')

  end subroutine sum_yintox_base

  subroutine sum_zintox_base(self, u, u_)
       !! As sum_yintox, for a z directional field.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call not_implemented('sum_zintox')

  end subroutine sum_zintox_base

  subroutine veccopy_base(self, dst, src)
       !! copy vectors: y = x
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: dst
    class(field_t), intent(in) :: src

    call not_implemented('veccopy')

  end subroutine veccopy_base

  subroutine vecadd_base(self, a, x, b, y)
       !! adds two vectors together: y = a*x + b*y
    implicit none

    class(base_backend_t) :: self
    real(dp), intent(in) :: a
    class(field_t), intent(in) :: x
    real(dp), intent(in) :: b
    class(field_t), intent(inout) :: y

    call not_implemented('vecadd')

  end subroutine vecadd_base

  subroutine vecmult_base(self, y, x)
      !! pointwise multiplication between two vectors: y(:) = y(:) * x(:)
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: y
    class(field_t), intent(in) :: x

    call not_implemented('vecmult')

  end subroutine vecmult_base

  real(dp) function scalar_product_base(self, x, y) result(s)
       !! Calculates the scalar product of two input fields
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(in) :: x, y

    s = 0._dp
    call not_implemented('scalar_product')

  end function scalar_product_base

  real(dp) function vector_norm_squared_base(self, a, b, c) &
    result(norm_squared)
    !! Computes the global discrete squared norm of a three-component field.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(in) :: a, b, c

    norm_squared = 0._dp
    call not_implemented('vector_norm_squared')

  end function vector_norm_squared_base

  subroutine field_max_mean_base(self, max_val, mean_val, f, &
                                 enforced_data_loc)
    !! Obtains maximum and mean values in a field
    implicit none

    class(base_backend_t) :: self
    real(dp), intent(out) :: max_val, mean_val
    class(field_t), intent(in) :: f
    integer, optional, intent(in) :: enforced_data_loc

    max_val = 0._dp
    mean_val = 0._dp
    call not_implemented('field_max_mean')

  end subroutine field_max_mean_base

  subroutine slice_max_sum_base(self, max_val, sum_val, f, &
                                i_slice, enforced_data_loc)
  !! Reduces a single slice of f at index i_slice along f's DIR axis.
  !! Returns signed max (not abs) and signed sum. No division by count.
  !! Caller is responsible for MPI_Allreduce across ranks.
    implicit none

    class(base_backend_t) :: self
    real(dp), intent(out) :: max_val, sum_val
    class(field_t), intent(in) :: f
    integer, intent(in) :: i_slice
    integer, optional, intent(in) :: enforced_data_loc

    max_val = 0._dp
    sum_val = 0._dp
    call not_implemented('slice_max_sum')

  end subroutine slice_max_sum_base

  subroutine field_plane_sums_base(self, sums, f)
    !! Sum a DIR_X field over x and z on every local y-plane:
    !! sums(j) is the sum over the plane of y index j. The order of
    !! summation is fixed, so the result is reproducible.
    implicit none

    class(base_backend_t) :: self
    real(dp), intent(out) :: sums(:)
    class(field_t), intent(in) :: f

    sums = 0._dp
    call not_implemented('field_plane_sums')

  end subroutine field_plane_sums_base

  subroutine field_scale_base(self, f, a)
    !! Scales a field by a
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    call not_implemented('field_scale')

  end subroutine field_scale_base

  subroutine field_shift_base(self, f, a)
    !! Shifts a field by a
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    call not_implemented('field_shift')

  end subroutine field_shift_base

  real(dp) function field_volume_integral_base(self, f) result(s)
    !! Reduces field to a scalar by integrating over the volume
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(in) :: f

    s = 0._dp
    call not_implemented('field_volume_integral')

  end function field_volume_integral_base

  subroutine field_set_face_base(self, f, c_start, c_end, face, &
                                 bc_start, bc_end, flow_rate_diff)
    !! A field is a subdomain with a rectangular cuboid shape.
    !! It has 6 faces, and these faces are either a subdomain boundary
    !! or a global domain boundary based on the location of the subdomain.
    !! This subroutine allows us to set any of these faces to a value,
    !! 'c_start' and 'c_end' for faces at opposite sides.
    !! 'face' is one of X_FACE, Y_FACE, Z_FACE from common.f90
    !! Optionally, bc_start/bc_end select the BC type (default BC_DIRICHLET).
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c_start, c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start
    integer, optional, intent(in) :: bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    call not_implemented('field_set_face')

  end subroutine field_set_face_base

  subroutine field_set_y_plane_base(self, f, c, plane)
    !! Set one interior y-plane of a DIR_X field to a constant.
    !!
    !! field_set_face reaches only the two boundary planes. The neutral ABL
    !! wall model needs the first plane above a no-slip floor, where the
    !! resolved gradient would otherwise contribute a second, spurious
    !! stress on top of the modelled one.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c
    integer, intent(in) :: plane !! 1-based y-vertex index

    call not_implemented('field_set_y_plane')

  end subroutine field_set_y_plane_base

  subroutine field_set_abl_wall_stress_base(self, stress, u, w, &
                                            sample_plane, stress_plane, &
                                            drag_coeff, component)
    !! Replace one SGS stress plane using the local sampled velocity.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: stress
    class(field_t), intent(in) :: u, w
    integer, intent(in) :: sample_plane, stress_plane, component
    real(dp), intent(in) :: drag_coeff

    call not_implemented('field_set_abl_wall_stress')

  end subroutine field_set_abl_wall_stress_base

  subroutine field_set_face_from_field_base(self, f, f_start, c_end, face, &
                                            bc_start, bc_end, flow_rate_diff)
    !! As field_set_face but with a spatially-varying inlet face field
    !! instead of a scalar c_start.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: f_start
    real(dp), intent(in) :: c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start
    integer, optional, intent(in) :: bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    call not_implemented('field_set_face_from_field')

  end subroutine field_set_face_from_field_base

  subroutine field_add_face_from_field_base(self, f, g, face, bc_start, bc_end)
    !! Add the face planes of `g` onto those of `f`, on Dirichlet faces
    !! only. Both fields are DIR_X; X_FACE and Y_FACE are supported.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: g
    integer, intent(in) :: face, bc_start, bc_end

    call not_implemented('field_add_face_from_field')

  end subroutine field_add_face_from_field_base

  subroutine compute_vorticity_base( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    !! Computes the vorticity magnitude from velocity gradients
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_vorticity')

  end subroutine compute_vorticity_base

  subroutine compute_qcriterion_base( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    !! Computes the Q-criterion from velocity gradients
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_qcriterion')

  end subroutine compute_qcriterion_base

  subroutine compute_smagorinsky_nut_base( &
    self, nut, mixing_length_sq, dudx, dudy, dudz, dvdx, dvdy, dvdz, &
    dwdx, dwdy, dwdz)
    !! Computes nut=l_s^2*sqrt(2*S_ij*S_ij) from velocity gradients.
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: nut
    class(field_t), intent(in) :: mixing_length_sq
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_smagorinsky_nut')

  end subroutine compute_smagorinsky_nut_base

  subroutine compute_sgs_stress_base( &
    self, stress, nut, gradient_a, gradient_b, scale_a, scale_b)
    !! Forms nut*(scale_a*gradient_a + scale_b*gradient_b).
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: stress
    class(field_t), intent(in) :: nut, gradient_a, gradient_b
    real(dp), intent(in) :: scale_a, scale_b

    call not_implemented('compute_sgs_stress')

  end subroutine compute_sgs_stress_base

  subroutine copy_data_to_f_base(self, f, data)
       !! Copy a regular 3D array in host memory into the specialist
       !! data structure field that lives on device or host
    implicit none

    class(base_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    call not_implemented('copy_data_to_f')

  end subroutine copy_data_to_f_base

  subroutine copy_f_to_data_base(self, data, f)
       !! Copy the specialist data structure from device or host back
       !! to a regular 3D data array in host memory.
    implicit none

    class(base_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    data = 0._dp
    call not_implemented('copy_f_to_data')

  end subroutine copy_f_to_data_base

  subroutine alloc_tdsops_base( &
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu, &
    filter_alpha &
    )
    !! Allocates the backend's own tdsops type and fills in the
    !! coefficients of the requested operation.
    implicit none

    class(base_backend_t) :: self
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

    call not_implemented('alloc_tdsops')

  end subroutine alloc_tdsops_base

  subroutine init_poisson_fft_base(self, mesh, xdirps, ydirps, zdirps, lowmem)
    !! Allocates and initialises the backend's FFT based Poisson solver.
    implicit none

    class(base_backend_t) :: self
    type(mesh_t), target, intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

    call not_implemented('init_poisson_fft')

  end subroutine init_poisson_fft_base

end module m_base_backend
