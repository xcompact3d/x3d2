!!! src/backend/omptgt/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Fields are device-resident: `omptgt_allocator_t` hands out
!! `omptgt_field_t`, whose storage is a raw `omp_target_alloc` pointer with no
!! host array behind it. Each operation therefore takes the field's
!! `get_dev_ptr()`, casts it to an array pointer with `c_f_pointer`, and drives
!! its loop nest inside a target region that names the pointer in an
!! `is_device_ptr` clause.
!!
!! Operations with no offloaded implementation yet call `not_implemented`,
!! which stops the run naming the operation.

module m_omptgt_backend

  use iso_c_binding, only: c_ptr, c_f_pointer

  use m_mpi, only: MPI_COMM_WORLD, MPI_SUM, MPI_Allreduce

  use m_common, only: dp, DIR_C, DIR_X, NULL_LOC, MPI_X3D2_DP, &
                      get_dirs_from_rdr

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t
  use m_field, only: field_t
  use m_ordering, only: get_index_reordering
  use m_tdsops, only: tdsops_t, dirps_t

  use m_omptgt_common, only: SZ
  use m_omptgt_allocator, only: omptgt_field_t

  implicit none

  type, extends(base_backend_t) :: omptgt_backend_t
  contains
    ! Offloaded operations
    procedure :: copy_f_to_data => copy_f_to_data_omptgt
    procedure :: copy_data_to_f => copy_data_to_f_omptgt
    procedure :: reorder => reorder_omptgt
    procedure :: vecadd => vecadd_omptgt
    procedure :: veccopy => veccopy_omptgt
    procedure :: vector_norm_squared => vector_norm_squared_omptgt
    procedure :: sync => sync_omptgt
    procedure :: get_device_bw_info => get_device_bw_info_omptgt
    ! Not offloaded yet
    procedure :: alloc_tdsops => alloc_tdsops_omptgt
    procedure :: transeq_x => transeq_x_omptgt
    procedure :: transeq_y => transeq_y_omptgt
    procedure :: transeq_z => transeq_z_omptgt
    procedure :: transeq_species => transeq_species_omptgt
    procedure :: tds_solve => tds_solve_omptgt
    procedure :: thom_solve => thom_solve_omptgt
    procedure :: sum_yintox => sum_yintox_omptgt
    procedure :: sum_zintox => sum_zintox_omptgt
    procedure :: vecmult => vecmult_omptgt
    procedure :: scalar_product => scalar_product_omptgt
    procedure :: field_max_mean => field_max_mean_omptgt
    procedure :: slice_max_sum => slice_max_sum_omptgt
    procedure :: field_scale => field_scale_omptgt
    procedure :: field_shift => field_shift_omptgt
    procedure :: field_volume_integral => field_volume_integral_omptgt
    procedure :: field_set_face => field_set_face_omptgt
    procedure :: field_set_face_from_field => &
      field_set_face_from_field_omptgt
    procedure :: compute_vorticity => compute_vorticity_omptgt
    procedure :: compute_qcriterion => compute_qcriterion_omptgt
    procedure :: compute_smagorinsky_nut => compute_smagorinsky_nut_omptgt
    procedure :: compute_sgs_stress => compute_sgs_stress_omptgt
    procedure :: init_poisson_fft => init_poisson_fft_omptgt
  end type

  interface omptgt_backend_t
    module procedure omptgt_backend_init
  end interface

  private
  public :: omptgt_backend_t

contains

  type(omptgt_backend_t) function omptgt_backend_init(mesh, allocator) &
    result(backend)
    !! Constructs the backend over a device-resident allocator.

    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator

    call backend%base_init()

    backend%allocator => allocator
    backend%mesh => mesh
  end function

  subroutine sync_omptgt(self)
    !! Waits for outstanding device work.
    !!
    !! Every target region in this backend is synchronous: none carries a
    !! `nowait` clause, so the host has already waited for the device by the
    !! time the region's enclosing call returns. That makes this a no-op
    !! rather than something unimplemented. It stops being one the moment an
    !! offloaded region here is made asynchronous.

    class(omptgt_backend_t) :: self

  end subroutine

  subroutine get_device_bw_info_omptgt(self, mem_clock_rt, mem_bus_width, &
                                       available)
    !! Reports no memory bandwidth figures: OpenMP offers no portable query
    !! for the memory clock or bus width of the target device, and this
    !! backend deliberately does not reach past OpenMP to a vendor runtime.

    class(omptgt_backend_t) :: self
    integer, intent(out) :: mem_clock_rt
    integer, intent(out) :: mem_bus_width
    logical, intent(out) :: available

    mem_clock_rt = 0
    mem_bus_width = 0
    available = .false.

  end subroutine

  subroutine veccopy_omptgt(self, dst, src)
    !! Copies `src` into `dst`. Both fields must be device-resident and share
    !! the same direction.

    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: dst
    class(field_t), intent(in) :: src

    if (src%dir /= dst%dir) then
      error stop "Called vector copy with incompatible fields"
    end if

    select type (dst)
    type is (omptgt_field_t)
      select type (src)
      type is (omptgt_field_t)
        call veccopy_offload_(dst%get_dev_ptr(), dst%get_shape(), &
                              src%get_dev_ptr(), src%get_shape())
      class default
        error stop "Called omptgt vector copy with unsupported source vector"
      end select
    class default
      error stop &
        "Called omptgt vector copy with unsupported destination vector"
    end select
  end subroutine

  subroutine veccopy_offload_(dst_ptr, n_dst, src_ptr, n_src)
    !! Offloaded kernel behind `veccopy_omptgt`, copying over the destination's
    !! extents.

    type(c_ptr), intent(in) :: dst_ptr, src_ptr
    integer, dimension(3), intent(in) :: n_dst, n_src

    real(dp), dimension(:, :, :), pointer :: dst, src
    integer :: i, j, k

    if (any(n_src < n_dst)) then
      error stop "Source field is smaller than the destination"
    end if

    call c_f_pointer(dst_ptr, dst, shape=n_dst)
    call c_f_pointer(src_ptr, src, shape=n_src)
    !$omp target is_device_ptr(dst_ptr, src_ptr)
    !$omp teams loop collapse(3)
    do k = 1, n_dst(3)
      do j = 1, n_dst(2)
        do i = 1, n_dst(1)
          dst(i, j, k) = src(i, j, k)
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  subroutine vecadd_omptgt(self, a, x, b, y)
    !! Computes y = a*x + b*y. Both fields must be device-resident.

    class(omptgt_backend_t) :: self
    real(dp), intent(in) :: a
    class(field_t), intent(in) :: x
    real(dp), intent(in) :: b
    class(field_t), intent(inout) :: y

    if (x%dir /= y%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    select type (x)
    type is (omptgt_field_t)
      select type (y)
      type is (omptgt_field_t)
        call vecadd_offload(self, a, x, b, y)
      class default
        error stop "Called omptgt vector add with unsupported result vector"
      end select
    class default
      error stop "Called omptgt vector add with unsupported source vector"
    end select

  end subroutine

  subroutine vecadd_offload(self, a, x, b, y)
    !! Device implementation of `vecadd`, which looks up the padded extents of
    !! the layout the fields are in.

    class(omptgt_backend_t) :: self
    real(dp), intent(in) :: a
    type(omptgt_field_t), intent(in) :: x
    real(dp), intent(in) :: b
    type(omptgt_field_t), intent(inout) :: y

    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(x%dir)

    call vecadd_offload_(dims, a, x%get_dev_ptr(), x%get_shape(), &
                         b, y%get_dev_ptr(), y%get_shape())

  end subroutine

  subroutine vecadd_offload_(dims, a, x_ptr, n_x, b, y_ptr, n_y)
    !! Offloaded kernel evaluating y = a*x + b*y over the padded domain.
    integer, dimension(3), intent(in) :: dims
    real(dp), intent(in) :: a
    type(c_ptr), intent(in) :: x_ptr
    integer, dimension(3), intent(in) :: n_x
    real(dp), intent(in) :: b
    type(c_ptr), intent(in) :: y_ptr
    integer, dimension(3), intent(in) :: n_y

    real(dp), dimension(:, :, :), pointer :: x, y
    integer :: i, j, k

    call c_f_pointer(x_ptr, x, shape=n_x)
    call c_f_pointer(y_ptr, y, shape=n_y)
    !$omp target is_device_ptr(x_ptr, y_ptr)
    !$omp teams loop collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          y(i, j, k) = a*x(i, j, k) + b*y(i, j, k)
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target
  end subroutine

  real(dp) function vector_norm_squared_omptgt(self, a, b, c) &
    result(norm_squared)
    !! Global sum of a**2 + b**2 + c**2, with one MPI reduction.

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: a, b, c

    real(dp) :: local_sum
    integer :: dims(3), ierr

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

    dims = self%mesh%get_dims(a%data_loc)

    select type (a)
    type is (omptgt_field_t)
      select type (b)
      type is (omptgt_field_t)
        select type (c)
        type is (omptgt_field_t)
          call vector_norm_squared_offload_( &
            local_sum, a%get_dev_ptr(), b%get_dev_ptr(), c%get_dev_ptr(), &
            a%get_shape(), dims)
        class default
          error stop "Called omptgt vector norm with unsupported vector"
        end select
      class default
        error stop "Called omptgt vector norm with unsupported vector"
      end select
    class default
      error stop "Called omptgt vector norm with unsupported vector"
    end select

    call MPI_Allreduce(local_sum, norm_squared, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function vector_norm_squared_omptgt

  subroutine vector_norm_squared_offload_(local_sum, a_dev, b_dev, c_dev, &
                                          n_fld, dims)
    !! Rank-local sum of a**2 + b**2 + c**2 over the physical points only.
    real(dp), intent(out) :: local_sum
    ! Named *_dev rather than *_ptr so the third one does not shadow the
    ! c_ptr type imported from iso_c_binding.
    type(c_ptr), intent(in) :: a_dev, b_dev, c_dev
    integer, dimension(3), intent(in) :: n_fld
    integer, dimension(3), intent(in) :: dims

    real(dp), dimension(:, :, :), pointer :: a, b, c
    integer :: i, j, k, k_i, k_j, n_i, stacked

    ! Pencils are stacked SZ points at a time along y, and a pencil group
    ! index runs fastest within a given z station.
    stacked = (dims(2) - 1)/SZ + 1

    call c_f_pointer(a_dev, a, shape=n_fld)
    call c_f_pointer(b_dev, b, shape=n_fld)
    call c_f_pointer(c_dev, c, shape=n_fld)

    local_sum = 0._dp
    !$omp target map(tofrom:local_sum) is_device_ptr(a_dev, b_dev, c_dev)
    !$omp teams loop collapse(3) reduction(+:local_sum) private(i, k, n_i)
    do k_j = 1, stacked
      do k_i = 1, dims(3)
        do j = 1, dims(1)
          k = k_j + (k_i - 1)*stacked
          ! The last group along y is partially filled with padding.
          n_i = min(SZ, dims(2) - (k_j - 1)*SZ)
          do i = 1, n_i
            local_sum = local_sum + a(i, j, k)**2 + b(i, j, k)**2 &
                        + c(i, j, k)**2
          end do
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  subroutine copy_data_to_f_omptgt(self, f, data)
    !! Copies a host array into a device-resident field.
    class(omptgt_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(f%dir)

    ! XXX: This could be improved following cuda/backend.f90:resolve_field_t()
    select type (f)
    type is (omptgt_field_t)
      call copy_data_to_f_omptgt_(f%get_dev_ptr(), f%get_shape(), data, dims)
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_data_to_f_omptgt

  subroutine copy_data_to_f_omptgt_(f_ptr, n_f, d, dims)
    !! Offloaded kernel behind `copy_data_to_f`. The host array is mapped in
    !! for the duration of the region, so this transfers over the bus.
    type(c_ptr), intent(in) :: f_ptr
    integer, dimension(3), intent(in) :: n_f
    real(dp), dimension(:, :, :), intent(in) :: d
    integer, dimension(3), intent(in) :: dims

    real(dp), dimension(:, :, :), pointer :: f_arr
    integer :: i, j, k

    call c_f_pointer(f_ptr, f_arr, shape=n_f)
    ! XXX: This could be improved following cuda/backend.f90:resolve_field_t()
    !$omp target map(to:d) is_device_ptr(f_ptr)
    !$omp teams loop collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          f_arr(i, j, k) = d(i, j, k)
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  subroutine copy_f_to_data_omptgt(self, data, f)
    !! Copies a device-resident field back into a host array.
    class(omptgt_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(f%dir)

    select type (f)
    type is (omptgt_field_t)
      call copy_f_to_data_omptgt_(data, f%get_dev_ptr(), f%get_shape(), dims)
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_f_to_data_omptgt

  subroutine copy_f_to_data_omptgt_(data, f_ptr, n_f, dims)
    !! Offloaded kernel behind `copy_f_to_data`. The host array is mapped back
    !! out of the region, so this transfers over the bus.
    real(dp), dimension(:, :, :), intent(out) :: data
    type(c_ptr), intent(in) :: f_ptr
    integer, dimension(3), intent(in) :: n_f
    integer, dimension(3), intent(in) :: dims

    real(dp), dimension(:, :, :), pointer :: f_arr
    integer :: i, j, k

    call c_f_pointer(f_ptr, f_arr, shape=n_f)
    !$omp target map(from:data) is_device_ptr(f_ptr)
    !$omp teams loop collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          data(i, j, k) = f_arr(i, j, k)
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  subroutine reorder_omptgt(self, u_, u, direction)
    !! Reorders `u` into `u_` between the two data layouts encoded in
    !! `direction`. Both fields are device-resident.
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, intent(in) :: direction
    integer, dimension(3) :: dims, cart_padded
    integer :: dir_from, dir_to

    dims = self%allocator%get_padded_dims(u%dir)
    cart_padded = self%allocator%get_padded_dims(DIR_C)
    call get_dirs_from_rdr(dir_from, dir_to, direction)

    ! XXX: This could be improved following cuda/backend.f90:resolve_field_t()
    select type (u_)
    type is (omptgt_field_t)
      select type (u)
      type is (omptgt_field_t)
        call reorder_omptgt_dd(u_%get_dev_ptr(), u_%get_shape(), &
                               u%get_dev_ptr(), u%get_shape(), dims, &
                               dir_from, dir_to, cart_padded)
      class default
        error stop "Called omptgt reorder with unsupported source field"
      end select
    class default
      error stop "Unsupported"
    end select

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omptgt

  subroutine reorder_point(u_, u, i, j, k, dir_from, dir_to, cart_padded)
    !$omp declare target
    !! Moves one point into its reordered slot.
    !!
    !! The reordered indices are locals of this routine rather than scalars
    !! privatised on the loop: an automatic local is per-call scratch by
    !! definition, so no private clause is needed. That matters because NVHPC
    !! miscompiles the private-clause spelling, handing get_index_reordering
    !! an invalid address for its intent(out) arguments, so the indices come
    !! back as garbage and the store off them faults.
    real(dp), dimension(:, :, :), intent(inout) :: u_
    real(dp), dimension(:, :, :), intent(in) :: u
    integer, intent(in) :: i, j, k
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    integer :: out_i, out_j, out_k

    call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                              dir_from, dir_to, SZ, cart_padded)
    u_(out_i, out_j, out_k) = u(i, j, k)

  end subroutine reorder_point

  subroutine reorder_omptgt_dd(u_ptr, n_u_, u_in_ptr, n_u, dims, dir_from, &
                               dir_to, cart_padded)
    !! Offloaded reordering of one device-resident field into another.
    type(c_ptr), intent(in) :: u_ptr, u_in_ptr
    integer, dimension(3), intent(in) :: n_u_, n_u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    real(dp), dimension(:, :, :), pointer :: u_, u
    integer :: i, j, k

    call c_f_pointer(u_ptr, u_, shape=n_u_)
    call c_f_pointer(u_in_ptr, u, shape=n_u)
    !$omp target is_device_ptr(u_ptr, u_in_ptr)
    !$omp teams loop collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call reorder_point(u_, u, i, j, k, dir_from, dir_to, cart_padded)
        end do
      end do
    end do
    !$omp end teams loop
    !$omp end target

  end subroutine

  subroutine not_implemented(operation)
    !! Stops the run, naming the operation that has no offloaded
    !! implementation yet. The name is printed rather than passed as the stop
    !! code because not every compiler here accepts a variable stop code.
    character(*), intent(in) :: operation

    print *, trim(operation)//': not implemented in the OMP_TGT backend yet'
    error stop 'not implemented in the OMP_TGT backend yet'

  end subroutine

  subroutine alloc_tdsops_omptgt( &
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    )
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

    call not_implemented('alloc_tdsops')

  end subroutine

  subroutine transeq_x_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_x')

  end subroutine

  subroutine transeq_y_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_y')

  end subroutine

  subroutine transeq_z_omptgt(self, du, dv, dw, u, v, w, nu, dirps)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps

    call not_implemented('transeq_z')

  end subroutine

  subroutine transeq_species_omptgt(self, dspec, uvw, spec, nu, dirps, sync)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: dspec
    class(field_t), intent(in) :: uvw, spec
    real(dp), intent(in) :: nu
    type(dirps_t), intent(in) :: dirps
    logical, intent(in) :: sync

    call not_implemented('transeq_species')

  end subroutine

  subroutine tds_solve_omptgt(self, du, u, tdsops)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    call not_implemented('tds_solve')

  end subroutine

  subroutine thom_solve_omptgt(self, du, u, tdsops)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: du
    class(field_t), intent(in) :: u
    class(tdsops_t), intent(in) :: tdsops

    call not_implemented('thom_solve')

  end subroutine

  subroutine sum_yintox_omptgt(self, u, u_)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call not_implemented('sum_yintox')

  end subroutine

  subroutine sum_zintox_omptgt(self, u, u_)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: u
    class(field_t), intent(in) :: u_

    call not_implemented('sum_zintox')

  end subroutine

  subroutine vecmult_omptgt(self, y, x)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: y
    class(field_t), intent(in) :: x

    call not_implemented('vecmult')

  end subroutine

  real(dp) function scalar_product_omptgt(self, x, y) result(s)
    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: x, y

    s = 0._dp
    call not_implemented('scalar_product')

  end function

  subroutine field_max_mean_omptgt(self, max_val, mean_val, f, &
                                   enforced_data_loc)
    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: max_val, mean_val
    class(field_t), intent(in) :: f
    integer, optional, intent(in) :: enforced_data_loc

    max_val = 0._dp
    mean_val = 0._dp
    call not_implemented('field_max_mean')

  end subroutine

  subroutine slice_max_sum_omptgt(self, max_val, sum_val, f, i_slice, &
                                  enforced_data_loc)
    class(omptgt_backend_t) :: self
    real(dp), intent(out) :: max_val, sum_val
    class(field_t), intent(in) :: f
    integer, intent(in) :: i_slice
    integer, optional, intent(in) :: enforced_data_loc

    max_val = 0._dp
    sum_val = 0._dp
    call not_implemented('slice_max_sum')

  end subroutine

  subroutine field_scale_omptgt(self, f, a)
    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    call not_implemented('field_scale')

  end subroutine

  subroutine field_shift_omptgt(self, f, a)
    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f
    real(dp), intent(in) :: a

    call not_implemented('field_shift')

  end subroutine

  real(dp) function field_volume_integral_omptgt(self, f) result(s)
    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: f

    s = 0._dp
    call not_implemented('field_volume_integral')

  end function

  subroutine field_set_face_omptgt(self, f, c_start, c_end, face, &
                                   bc_start, bc_end, flow_rate_diff)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    real(dp), intent(in) :: c_start, c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start
    integer, optional, intent(in) :: bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    call not_implemented('field_set_face')

  end subroutine

  subroutine field_set_face_from_field_omptgt(self, f, f_start, c_end, face, &
                                              bc_start, bc_end, &
                                              flow_rate_diff)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: f_start
    real(dp), intent(in) :: c_end
    integer, intent(in) :: face
    integer, optional, intent(in) :: bc_start
    integer, optional, intent(in) :: bc_end
    real(dp), optional, intent(in) :: flow_rate_diff

    call not_implemented('field_set_face_from_field')

  end subroutine

  subroutine compute_vorticity_omptgt( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_vorticity')

  end subroutine

  subroutine compute_qcriterion_omptgt( &
    self, field_out, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: field_out
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_qcriterion')

  end subroutine

  subroutine compute_smagorinsky_nut_omptgt( &
    self, nut, mixing_length_sq, dudx, dudy, dudz, dvdx, dvdy, dvdz, &
    dwdx, dwdy, dwdz)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: nut
    class(field_t), intent(in) :: mixing_length_sq
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz

    call not_implemented('compute_smagorinsky_nut')

  end subroutine

  subroutine compute_sgs_stress_omptgt( &
    self, stress, nut, gradient_a, gradient_b, scale_a, scale_b)
    class(omptgt_backend_t) :: self
    class(field_t), intent(inout) :: stress
    class(field_t), intent(in) :: nut, gradient_a, gradient_b
    real(dp), intent(in) :: scale_a, scale_b

    call not_implemented('compute_sgs_stress')

  end subroutine

  subroutine init_poisson_fft_omptgt(self, mesh, xdirps, ydirps, zdirps, &
                                     lowmem)
    class(omptgt_backend_t) :: self
    type(mesh_t), target, intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

    call not_implemented('init_poisson_fft')

  end subroutine

end module
