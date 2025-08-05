!!! src/backend/omp/target/backend.f90
!!
!! OpenMP target offload backend implementation.
!!
!! Note this extends the CPU (host) OpenMP backend with the intention of being able to use fallback implementations where necessary.

module m_omptgt_backend

  use m_common, only: dp

  use m_allocator, only: allocator_t
  use m_mesh, only: mesh_t
  use m_field, only: field_t

  use m_omp_backend, only: omp_backend_t

  use m_omptgt_allocator, only: omptgt_field_t

  implicit none

  type, extends(omp_backend_t) :: omptgt_backend_t
  contains
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

    !$omp target teams distribute parallel do collapse(3) map(to : x) map(tofrom : y)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          y%data_tgt(i, j, k) = a * x%data_tgt(i, j, k) + b * y%data_tgt(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine

end module
