module m_time_integrator
   use m_allocator, only: allocator_t, field_t, flist_t
   use m_base_backend, only: base_backend_t
   implicit none

   private

   type :: time_intg_t
      integer :: istep, nsteps, nsubsteps, order, nvars, nolds
      type(flist_t), allocatable :: olds(:,:)
      class(base_backend_t), pointer :: backend
      class(allocator_t), pointer :: allocator
      class(field_t), pointer :: u, v, w, du, dv, dw
   contains
      procedure :: run
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

      constructor%u => allocator%get_block()
      constructor%v => allocator%get_block()
      constructor%w => allocator%get_block()

      constructor%du => allocator%get_block()
      constructor%dv => allocator%get_block()
      constructor%dw => allocator%get_block()

      if (present(nvars)) then
         constructor%nvars = nvars
      else
         constructor%nvars = 3
      end if

      constructor%nolds = 1

      allocate(constructor%olds(constructor%nvars, constructor%nolds))

      ! Request all the storage for old timesteps
      do i = 1, constructor%nvars
         do j = 1, constructor%nolds
            constructor%olds(i, j)%ptr => allocator%get_block()
         end do
      end do

   end function constructor

   subroutine run(self, n_iter)
      implicit none

      class(time_intg_t), intent(in) :: self
      integer, intent(in) :: n_iter

      print*, 'start run'
      call self%backend%transeq(self%du, self%dv, self%dw, &
                                self%u, self%v, self%w)

      ! time integration
      print*, 'run'

   end subroutine run

   subroutine adams_bashford_1st(vels, olds, coeffs)
      type(flist_t) :: vels(:), olds(:)
      real :: coeffs(:)

      !call vec_add(vels, olds, coeffs)
   end subroutine adams_bashford_1st

end module m_time_integrator
