!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

  use iso_c_binding, only: c_ptr, c_f_pointer

  use m_common, only: dp, DIR_C, get_dirs_from_rdr

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
        call veccopy_offload_(dst%get_dev_ptr(), dst%get_shape(), src%get_dev_ptr(), src%get_shape())
      class default
        error stop "Called omptgt vector copy with unsupported source vector"
      end select
    class default
      error stop &
        "Called omptgt vector copy with unsupported destination vector"
    end select
  end subroutine

  subroutine veccopy_offload_(cp_dst, n_dst, cp_src, n_src)
    type(c_ptr), intent(inout) :: cp_dst
    type(c_ptr), intent(in) :: cp_src
    integer, dimension(3), intent(in) :: n_dst, n_src

    real(dp), dimension(:, :, :), pointer :: dst
    real(dp), dimension(:, :, :), pointer :: src

    integer :: i, j, k

    if (any(n_src < n_dst)) then
      error stop "SRC array is smaller than destination"
    end if

    !$omp target teams is_device_ptr(cp_dst, cp_src)
    call c_f_pointer(cp_dst, dst, shape=n_dst)
    call c_f_pointer(cp_src, src, shape=n_src)
    !$omp loop collapse(3)
    do k = 1, n_dst(3)
      do j = 1, n_dst(2)
        do i = 1, n_dst(1)
          dst(i, j, k) = src(i, j, k)
        end do
      end do
    end do
    !$omp end loop
    !$omp end target teams

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

    type(c_ptr) :: cp_x, cp_y
    real(dp), dimension(:,:,:), pointer :: p_x, p_y

    integer, dimension(3) :: dims
    integer :: i, j, k

    dims = self%allocator%get_padded_dims(x%dir)

    cp_x = x%get_dev_ptr()
    cp_y = y%get_dev_ptr()

    !$omp target teams is_device_ptr(cp_x, cp_y)
    call c_f_pointer(cp_x, p_x, shape=x%get_shape())
    call c_f_pointer(cp_y, p_y, shape=y%get_shape())
    !$omp loop collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          p_y(i, j, k) = a*p_x(i, j, k) + b*p_y(i, j, k)
        end do
      end do
    end do
    !$omp end loop
    !$omp end target teams

  end subroutine

  subroutine copy_data_to_f_omptgt(self, f, data)
    class(omptgt_backend_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    type(c_ptr) :: f_dev_ptr
    real(dp), dimension(:,:,:), pointer :: p_f
    integer, dimension(3) :: dims
    integer :: i, j, k

    dims = self%allocator%get_padded_dims(f%dir)

    ! XXX: This could be improved following cuda/backend.f90:resolve_field_t()
    select type (f)
    type is (omptgt_field_t)
      f_dev_ptr = f%get_dev_ptr()
      !$omp target teams map(to:data) is_device_ptr(f_dev_ptr)
      call c_f_pointer(f_dev_ptr, p_f, shape=f%get_shape())
      !$omp loop collapse(3)
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            p_f(i, j, k) = data(i, j, k)
          end do
        end do
      end do
      !$omp end loop
      !$omp end target teams
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_data_to_f_omptgt

  subroutine copy_f_to_data_omptgt(self, data, f)
    class(omptgt_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(field_t), intent(in) :: f

    type(c_ptr) :: f_dev_ptr
    real(dp), dimension(:,:,:), pointer :: p_f
    integer, dimension(3) :: dims
    integer :: i, j, k

    dims = self%allocator%get_padded_dims(f%dir)

    select type (f)
    type is (omptgt_field_t)
      f_dev_ptr = f%get_dev_ptr()
      !$omp target teams map(from:data) is_device_ptr(f_dev_ptr)
      call c_f_pointer(f_dev_ptr, p_f, shape=f%get_shape())
      !$omp loop collapse(3)
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            data(i, j, k) = p_f(i, j, k)
          end do
        end do
      end do
      !$omp end loop
      !$omp end target teams
    class default
      error stop "Unsupported"
    end select

  end subroutine copy_f_to_data_omptgt

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
        call reorder_omptgt_dd(u_, u, dims, dir_from, &
                               dir_to, cart_padded)
      class default
        call reorder_omptgt_dh(u_, u%data, dims, dir_from, dir_to, &
                               cart_padded)
      end select
    class default
      error stop "Unsupported"
    end select

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omptgt

  subroutine reorder_omptgt_dd(u_, u, dims, dir_from, dir_to, cart_padded)
    type(omptgt_field_t) :: u_
    type(omptgt_field_t), intent(in) :: u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    type(c_ptr) :: cp_u_, cp_u
    real(dp), dimension(:,:,:), pointer :: p_u_, p_u
    integer :: i, j, k
    integer :: out_i, out_j, out_k

    cp_u_ = u_%get_dev_ptr()
    cp_u = u%get_dev_ptr()
    !$omp target teams is_device_ptr(cp_u_, cp_u)
    call c_f_pointer(cp_u_, p_u_, shape=u_%get_shape())
    call c_f_pointer(cp_u, p_u, shape=u%get_shape())
    !$omp loop collapse(3) private(out_i, out_j, out_k)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          p_u_(out_i, out_j, out_k) = p_u(i, j, k)
        end do
      end do
    end do
    !$omp end loop
    !$omp end target teams

  end subroutine

  subroutine reorder_omptgt_dh(u_, u, dims, dir_from, dir_to, cart_padded)
    type(omptgt_field_t) :: u_
    real(dp), dimension(:, :, :), pointer, intent(in) :: u
    integer, dimension(3), intent(in) :: dims
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3), intent(in) :: cart_padded

    type(c_ptr) :: cp_u_
    real(dp), dimension(:,:,:), pointer :: p_u_
    integer :: i, j, k
    integer :: out_i, out_j, out_k

    cp_u_ = u_%get_dev_ptr()
    !$omp target teams map(to:u) is_device_ptr(cp_u_)
    call c_f_pointer(cp_u_, p_u_, shape=u_%get_shape())
    !$omp loop collapse(3) private(out_i, out_j, out_k)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          p_u_(out_i, out_j, out_k) = u(i, j, k)
        end do
      end do
    end do
    !$omp end loop
    !$omp end target teams

  end subroutine

end module
