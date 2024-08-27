module m_case
  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, DIR_X
  use m_solver, only: solver_t
  use m_time_integrator, only: time_intg_t

  implicit none

  type :: case_t
    !! Base case type
    !!
    !! Each specific case extends this base class and overwrites its procedures
    class(base_backend_t), pointer :: backend
    class(solver_t), pointer :: solver
    !class(time_intg_t) :: time_integrator
    type(allocator_t) :: host_allocator
  contains
    procedure :: set_ICs
    procedure :: impose_BCs
    procedure :: post_transeq
    procedure :: postprocess
    procedure :: run
  end type case_t

  interface case_t
    module procedure init
  end interface case_t 

contains

  function init(solver) result(case)
    implicit none

    class(solver_t), target, intent(inout) :: solver
    type(case_t) :: case

    case%solver => solver

    call case%set_ICs()

  end function init

  subroutine set_ICs(self)
    implicit none

    class(case_t) :: self

  end subroutine set_ICs

  subroutine impose_BCs(self)
    implicit none

    class(case_t) :: self

  end subroutine impose_BCs

  subroutine post_transeq(self, du, dv, dw)
    !! Use this for forcings; e.g. constant pressure gradient
    implicit none

    class(case_t) :: self
    class(field_t), intent(inout) :: du, dv, dw

  end subroutine post_transeq

  subroutine postprocess(self, t)
    implicit none

    class(case_t) :: self
    real(dp), intent(in) :: t

    ! for example call enstrophy from vector_calculus
  end subroutine postprocess

  subroutine run(self)
    implicit none

    class(case_t) :: self

    class(field_t), pointer :: du, dv, dw

    real(dp) :: t
    integer :: i, j

    do i = 1, self%solver%n_iters
      do j = 1, self%solver%time_integrator%nstage
        du => self%solver%backend%allocator%get_block(DIR_X)
        dv => self%solver%backend%allocator%get_block(DIR_X)
        dw => self%solver%backend%allocator%get_block(DIR_X)

        call self%solver%transeq(du, dv, dw, self%solver%u, self%solver%v, self%solver%w)

        call self%post_transeq(du, dv, dw)

        ! time integration
        call self%solver%time_integrator%step( &
          self%solver%u, self%solver%v, self%solver%w, du, dv, dw, self%solver%dt &
          )

        call self%solver%backend%allocator%release_block(du)
        call self%solver%backend%allocator%release_block(dv)
        call self%solver%backend%allocator%release_block(dw)

        ! impose boundary conditions
        call self%impose_BCs()

        ! pressure
        call self%solver%pressure_correction()
      end do

      if (mod(i, self%solver%n_output) == 0) then
        t = i*self%solver%dt
        call self%postprocess(t)
      end if
    end do

  end subroutine run

end module m_case
