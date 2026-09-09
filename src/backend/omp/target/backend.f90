!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

  use m_common, only: dp, DIR_C, get_dirs_from_rdr, &
                      RDR_C2X, RDR_X2C, RDR_X2Y, RDR_Y2X, RDR_X2Z, RDR_Z2X

  use m_allocator, only: allocator_t
  use m_mesh, only: mesh_t
  use m_field, only: field_t
  use m_ordering, only: get_index_reordering

  use m_omp_common, only: SZ
  use m_omp_backend, only: omp_backend_t

  use m_omptgt_allocator, only: omptgt_field_t
  use m_omptgt_kernels_reorder, only: reorder_omptgt_c2x, reorder_omptgt_x2c, &
                                      reorder_omptgt_x2y, reorder_omptgt_y2x, &
                                      reorder_omptgt_x2z, reorder_omptgt_z2x

  implicit none

  type, extends(omp_backend_t) :: omptgt_backend_t
  contains
    procedure :: copy_f_to_data => copy_f_to_data_omptgt
    procedure :: copy_data_to_f => copy_data_to_f_omptgt
    procedure :: reorder => reorder_omptgt
    procedure :: vecadd => vecadd_omptgt
    procedure :: veccopy => veccopy_omptgt
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
        call reorder_offload(u_%data_tgt, u%data_tgt, direction, dims, &
                             dir_from, dir_to, cart_padded)
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

  subroutine reorder_offload(u_, u, direction, dims, dir_from, dir_to, &
                             cart_padded)
    !! Picks a kernel for a device to device reorder.
    !!
    !! Directional storages differ in which cartesian direction their leading
    !! dimension runs along: DIR_C, DIR_Y and DIR_Z run theirs along x, DIR_X
    !! runs its along y. A reorder between DIR_X and any other direction
    !! therefore transposes, and gets a kernel that stages the data through a
    !! tile in team-shared memory so that both the read and the write stay
    !! contiguous. The rest are contiguous on both sides already and go
    !! through the generic index map.
    real(dp), dimension(:, :, :), pointer :: u_
    real(dp), dimension(:, :, :), pointer, intent(in) :: u
    integer, intent(in) :: direction
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    integer :: nx, ny, nz

    nx = cart_padded(1); ny = cart_padded(2); nz = cart_padded(3)

    select case (direction)
    case (RDR_C2X)
      call reorder_omptgt_c2x(u_, u, nx, ny, nz)
    case (RDR_X2C)
      call reorder_omptgt_x2c(u_, u, nx, ny, nz)
    case (RDR_X2Y)
      call reorder_omptgt_x2y(u_, u, nx, ny, nz)
    case (RDR_Y2X)
      call reorder_omptgt_y2x(u_, u, nx, ny, nz)
    case (RDR_X2Z)
      call reorder_omptgt_x2z(u_, u, nx, ny, nz)
    case (RDR_Z2X)
      call reorder_omptgt_z2x(u_, u, nx, ny, nz)
    case default
      call reorder_omptgt_dd(u_, u, dims, dir_from, dir_to, cart_padded)
    end select

  end subroutine reorder_offload

  subroutine reorder_omptgt_dd(u_, u, dims, dir_from, dir_to, cart_padded)
    !! Generic device to device reorder: walks the input and maps each index
    !! individually. Used for the reorders that need no transpose.
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
