!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

  use mpi

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, NULL_LOC, MPI_X3D2_DP, &
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
    procedure :: scalar_product => scalar_product_omptgt
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
        call veccopy_offload_(dst%data_tgt, src%data_tgt)
      class default
        error stop "Called omptgt vector copy with unsupported source vector"
      end select
    class default
      error stop &
        "Called omptgt vector copy with unsupported destination vector"
    end select
  end subroutine

  subroutine veccopy_offload_(dst, src)

    real(dp), dimension(:, :, :), intent(inout) :: dst
    real(dp), dimension(:, :, :), intent(in) :: src

    integer, dimension(3) :: n
    integer :: i, j, k

    n = shape(dst)

    !$omp target teams loop collapse(3) has_device_addr(dst, src)
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          dst(i, j, k) = src(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop

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

    call vecadd_offload_(dims, a, x%data_tgt, b, y%data_tgt)

  end subroutine

  subroutine vecadd_offload_(dims, a, x, b, y)
    integer, dimension(3), intent(in) :: dims
    real(dp), intent(in) :: a
    real(dp), dimension(:, :, :), intent(in) :: x
    real(dp), intent(in) :: b
    real(dp), dimension(:, :, :), intent(inout) :: y

    integer :: i, j, k

    !$omp target teams loop collapse(3) has_device_addr(x, y)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          y(i, j, k) = a*x(i, j, k) + b*y(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop
  end subroutine

  real(dp) function scalar_product_omptgt(self, x, y) result(s)
    implicit none

    class(omptgt_backend_t) :: self
    class(field_t), intent(in) :: x, y

    real(dp) :: local_sum
    integer :: dims(3), dims_padded(3), n, n_i, n_i_pad, n_j, ierr

    if ((x%data_loc == NULL_LOC) .or. (y%data_loc == NULL_LOC)) then
      error stop "You must set the data_loc before calling scalar product"
    end if
    if ((x%data_loc /= y%data_loc) .or. (x%dir /= y%dir)) then
      error stop "Called scalar product with incompatible fields"
    end if

    dims = self%mesh%get_dims(x%data_loc)
    dims_padded = self%allocator%get_padded_dims(DIR_C)

    if (x%dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (x%dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (x%dir == DIR_Z) then
      n = dims(3); n_j = dims(2); n_i = dims(1); n_i_pad = dims_padded(1)
    else
      error stop 'scalar_product_cuda does not support DIR_C fields!'
    end if

    select type (x)
    type is (omptgt_field_t)
      select type (y)
      type is (omptgt_field_t)
        call scalar_product_offload_(local_sum, x%data_tgt, y%data_tgt, &
                n, n_j, n_i, n_i_pad)
      class default
        error stop "Called omptgt vector copy with unsupported source vector"
      end select
    class default
      error stop "Called omptgt vector copy with unsupported source vector"
    end select

    call MPI_Allreduce(local_sum, s, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function scalar_product_omptgt

  subroutine scalar_product_offload_(local_sum, x, y, n, n_i, n_i_pad, n_j)
    real(dp), intent(out) :: local_sum
    real(dp), dimension(:, :, :), intent(in) :: x, y
    integer, intent(in) :: n, n_i, n_i_pad, n_j

    integer :: i, j, k, k_i, k_j

    local_sum = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:local_sum) private(j, k) &
    !$omp has_device_addr(x, y)
    do k_j = 1, (n_j - 1)/SZ + 1
      do k_i = 1, n_i
        do i = 1, SZ
          k = k_j + (k_i - 1)*((n_j - 1)/SZ + 1)
          do j = 1, n
            local_sum = local_sum + x(i, j, k) * y(i, j, k)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
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
          call vector_norm_squared_offload_(local_sum, a%data_tgt, &
                                            b%data_tgt, c%data_tgt, dims)
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

  subroutine vector_norm_squared_offload_(local_sum, a, b, c, dims)
    !! Rank-local sum of a**2 + b**2 + c**2 over the physical points only.
    real(dp), intent(out) :: local_sum
    real(dp), dimension(:, :, :), intent(in) :: a, b, c
    integer, dimension(3), intent(in) :: dims

    integer :: i, j, k, k_i, k_j, n_i, stacked

    ! Pencils are stacked SZ points at a time along y, and a pencil group
    ! index runs fastest within a given z station.
    stacked = (dims(2) - 1)/SZ + 1

    local_sum = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:local_sum) private(i, k, n_i) &
    !$omp has_device_addr(a, b, c)
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
    !$omp end target teams distribute parallel do

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
      call copy_data_to_f_omptgt_(f%data_tgt, data, dims)
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_data_to_f_omptgt

  subroutine copy_data_to_f_omptgt_(f_arr, d, dims)
    real(dp), dimension(:, :, :), intent(inout) :: f_arr
    real(dp), dimension(:, :, :), intent(in) :: d
    integer, dimension(3), intent(in) :: dims

    integer :: i, j, k

    ! XXX: This could be improved following cuda/backend.f90:resolve_field_t()
    !$omp target teams loop collapse(3) map(to:d) has_device_addr(f_arr)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          f_arr(i, j, k) = d(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop

  end subroutine

  subroutine copy_f_to_data_omptgt(self, data, f)
    class(omptgt_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(f%dir)

    select type (f)
    type is (omptgt_field_t)
      call copy_f_to_data_omptgt_(data, f%data_tgt, dims)
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_f_to_data_omptgt

  subroutine copy_f_to_data_omptgt_(data, f_arr, dims)
    real(dp), dimension(:, :, :), intent(out) :: data
    real(dp), dimension(:, :, :), intent(in) :: f_arr
    integer, dimension(3), intent(in) :: dims

    integer :: i, j, k

    !$omp target teams loop collapse(3) map(from:data) has_device_addr(f_arr)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          data(i, j, k) = f_arr(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop

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
        call reorder_omptgt_dd(u_%data_tgt, u%data_tgt, dims, dir_from, &
                               dir_to, cart_padded)
      class default
        call reorder_omptgt_dh(u_%data_tgt, u%data, dims, dir_from, dir_to, &
                               cart_padded)
      end select
    class default
      error stop "Unsupported"
    end select

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omptgt

  subroutine reorder_omptgt_dd(u_, u, dims, dir_from, dir_to, cart_padded)
    real(dp), dimension(:, :, :), pointer :: u_
    real(dp), dimension(:, :, :), pointer, intent(in) :: u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    integer :: i, j, k
    integer :: out_i, out_j, out_k

    !$omp target teams loop collapse(3) private(out_i, out_j, out_k) has_device_addr(u_, u)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u_(out_i, out_j, out_k) = u(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop

  end subroutine

  subroutine reorder_omptgt_dh(u_, u, dims, dir_from, dir_to, cart_padded)
    real(dp), dimension(:, :, :), pointer :: u_
    real(dp), dimension(:, :, :), pointer, intent(in) :: u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    integer :: i, j, k
    integer :: out_i, out_j, out_k

    !$omp target teams loop collapse(3) private(out_i, out_j, out_k) map(to:u) has_device_addr(u_)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u_(out_i, out_j, out_k) = u(i, j, k)
        end do
      end do
    end do
    !$omp end target teams loop

  end subroutine

end module
