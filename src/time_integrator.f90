module m_time_integrator
   use m_allocator, only: allocator_t, field_t, flist_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp

   implicit none

   type :: time_intg_t
      integer :: istep, nsteps, nsubsteps, order, nvars, nolds
      type(flist_t), allocatable :: olds(:,:)
      class(base_backend_t), pointer :: backend
      class(allocator_t), pointer :: allocator
   contains
      procedure :: step
   end type time_intg_t

   interface time_intg_t
      module procedure constructor
   end interface time_intg_t

contains

   function constructor(backend, allocator, nvars)
      implicit none

      type(time_intg_t) :: constructor
      class(base_backend_t), pointer :: backend
      class(allocator_t), pointer :: allocator
      integer, intent(in), optional :: nvars

      integer :: i, j

      constructor%backend => backend
      constructor%allocator => allocator

      if (present(nvars)) then
         constructor%nvars = nvars
      else
         constructor%nvars = 3
      end if

      constructor%nolds = 0

      allocate(constructor%olds(constructor%nvars, constructor%nolds))

      ! Request all the storage for old timesteps
      do i = 1, constructor%nvars
         do j = 1, constructor%nolds
            constructor%olds(i, j)%ptr => allocator%get_block()
         end do
      end do

   end function constructor

   subroutine step(self, u, v, w, du, dv, dw, dt)
      implicit none

      class(time_intg_t), intent(in) :: self
      class(field_t), intent(inout) :: u, v, w
      class(field_t), intent(in) :: du, dv, dw

      real(dp), intent(in) :: dt

      call self%backend%vecadd(dt, du, 1._dp, u)
      call self%backend%vecadd(dt, dv, 1._dp, v)
      call self%backend%vecadd(dt, dw, 1._dp, w)

   end subroutine step

   subroutine adams_bashford_1st(vels, olds, coeffs)
      type(flist_t) :: vels(:), olds(:)
      real :: coeffs(:)

      !call vec_add(vels, olds, coeffs)
   end subroutine adams_bashford_1st

end module m_time_integrator
