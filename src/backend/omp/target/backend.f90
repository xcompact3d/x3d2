!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

  use iso_c_binding, only: c_ptr, c_f_pointer

  use mpi

  use m_common, only: dp, DIR_C, DIR_X, NULL_LOC, MPI_X3D2_DP, &
                      get_dirs_from_rdr

  use m_allocator, only: allocator_t
  use m_mesh, only: mesh_t
  use m_field, only: field_t
  use m_ordering, only: get_index_reordering

  use m_omp_common, only: SZ
  use m_omp_backend, only: omp_backend_t

  use m_omptgt_allocator, only: omptgt_field_t

  implicit none

  type, extends(omp_backend_t) :: omptgt_backend_t
  contains
    procedure :: copy_f_to_data => copy_f_to_data_omptgt
    procedure :: copy_data_to_f => copy_data_to_f_omptgt
    procedure :: reorder => reorder_omptgt
    procedure :: vecadd => vecadd_omptgt
    procedure :: veccopy => veccopy_omptgt
    procedure :: vector_norm_squared => vector_norm_squared_omptgt
  end type

  interface omptgt_backend_t
    module procedure omptgt_backend_init
  end interface

  private
  public :: omptgt_backend_t

contains

  type(omptgt_backend_t) function omptgt_backend_init(mesh, allocator) &
    result(backend)

    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator

    backend%omp_backend_t = omp_backend_t(mesh, allocator)
  end function

  subroutine veccopy_omptgt(self, dst, src)

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
        error stop "Device/host fallback not yet implemented"
      end select
    class default
      call self%omp_backend_t%vecadd(a, x, b, y)
    end select

  end subroutine

  subroutine vecadd_offload(self, a, x, b, y)

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
          error stop "Called omptgt vector copy with unsupported source vector"
        end select
      class default
        error stop "Called omptgt vector copy with unsupported source vector"
      end select
    class default
      error stop "Called omptgt vector copy with unsupported source vector"
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
        call reorder_omptgt_dh(u_%get_dev_ptr(), u_%get_shape(), u%data, &
                               dims, dir_from, dir_to, cart_padded)
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

  subroutine reorder_omptgt_dh(u_ptr, n_u_, u, dims, dir_from, dir_to, &
                               cart_padded)
    type(c_ptr), intent(in) :: u_ptr
    integer, dimension(3), intent(in) :: n_u_
    real(dp), dimension(:, :, :), pointer, intent(in) :: u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    real(dp), dimension(:, :, :), pointer :: u_
    integer :: i, j, k

    call c_f_pointer(u_ptr, u_, shape=n_u_)
    !$omp target map(to:u) is_device_ptr(u_ptr)
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

end module
