!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

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
  end type

  interface omptgt_backend_t
    module procedure omptgt_backend_init
  end interface

  private
  public :: omptgt_backend_t

contains

  type(omptgt_backend_t) function omptgt_backend_init(mesh, allocator) result(backend)
    
    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), target, intent(inout) :: allocator

    backend%omp_backend_t = omp_backend_t(mesh, allocator)
  end function

  subroutine vecadd_omptgt(self, a, x, b, y)

    class(omptgt_backend_t) :: self
    real(dp), intent(in) :: a
    class(field_t), intent(in) :: x
    real(dp), intent(in) :: b
    class(field_t), intent(inout) :: y

    if (x%dir /= y%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    select type(x)
    type is(omptgt_field_t)
      select type(y)
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

    integer :: i, j, k
    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(x%dir)

    !$omp target teams distribute parallel do default(shared) private(i, j, k) collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          y%data_tgt(i, j, k) = a * x%data_tgt(i, j, k) + b * y%data_tgt(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine

  subroutine copy_data_to_f_omptgt(self, f, data)
    class(omptgt_backend_t), intent(inout) :: self
    class(omptgt_field_t), intent(inout) :: f
    real(dp), dimension(:, :, :), intent(in) :: data

    integer :: i, j, k
    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(f%dir)
    print *, "-->", dims

    !$omp target teams distribute parallel do default(shared) private(i, j, k) collapse(3) map(to:data)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          f%data_tgt(i, j, k) = data(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine copy_data_to_f_omptgt
    
  subroutine copy_f_to_data_omptgt(self, data, f)
    class(omptgt_backend_t), intent(inout) :: self
    real(dp), dimension(:, :, :), intent(out) :: data
    class(omptgt_field_t), intent(in) :: f

    integer :: i, j, k
    integer, dimension(3) :: dims

    dims = self%allocator%get_padded_dims(f%dir)
    print *, "<--", dims

    print *, minval(data), maxval(data)
    !$omp target teams distribute parallel do default(shared) private(i, j, k) collapse(3) map(from:data)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          data(i, j, k) = f%data_tgt(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    print *, minval(data), maxval(data)
    
  end subroutine copy_f_to_data_omptgt

  subroutine reorder_omptgt(self, u_, u, direction)
    class(omptgt_backend_t) :: self
    class(omptgt_field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, intent(in) :: direction
    integer, dimension(3) :: dims, cart_padded
    integer :: i, j, k
    integer :: out_i, out_j, out_k
    integer :: dir_from, dir_to

    dims = self%allocator%get_padded_dims(u%dir)
    cart_padded = self%allocator%get_padded_dims(DIR_C)
    call get_dirs_from_rdr(dir_from, dir_to, direction)

    select type(u)
    type is(omptgt_field_t)
      call reorder_omptgt_dd(u_, u, dims, cart_padded, dir_from, dir_to)
    class default
      call reorder_omptgt_hd(u_, u, dims, cart_padded, dir_from, dir_to)
    end select

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omptgt

  subroutine reorder_omptgt_dd(u_, u, dims, cart_padded, dir_from, dir_to)
    class(omptgt_field_t), intent(inout) :: u_
    class(omptgt_field_t), intent(in) :: u
    integer, dimension(3) :: dims, cart_padded
    integer :: dir_from, dir_to

    integer :: i, j, k
    integer :: out_i, out_j, out_k

    !$omp target teams distribute parallel do private(out_i, out_j, out_k) collapse(3)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u_%data_tgt(out_i, out_j, out_k) = u%data_tgt(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine

  subroutine reorder_omptgt_hd(u_, u, dims, cart_padded, dir_from, dir_to)
    class(omptgt_field_t), intent(inout) :: u_
    class(field_t), intent(in) :: u
    integer, dimension(3) :: dims, cart_padded
    integer :: dir_from, dir_to

    integer :: i, j, k
    integer :: out_i, out_j, out_k

    !$omp target teams distribute parallel do private(out_i, out_j, out_k) collapse(3) map(to:u%data)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u_%data_tgt(out_i, out_j, out_k) = u%data(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine

end module
