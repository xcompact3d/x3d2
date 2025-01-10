module m_base_case
  !! Provides the base case for running a simulation. New cases are
  !! implemented by extending this to specify the initial and boundary
  !! conditions, forcing terms and case-specific postprocessing and analysis.
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Z, DIR_C, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: solver_t, init

  implicit none

  type, abstract :: base_case_t
    class(solver_t), allocatable :: solver
  contains
    procedure(boundary_conditions), deferred :: boundary_conditions
    procedure(initial_conditions), deferred :: initial_conditions
    procedure(post_transeq), deferred :: post_transeq
    procedure(postprocess), deferred :: postprocess
    procedure :: case_init
    procedure :: run
    procedure :: print_enstrophy
    procedure :: print_div_max_mean
  end type base_case_t

  abstract interface
    subroutine boundary_conditions(self)
      !! Applies case-specific boundary coinditions
      import :: base_case_t
      implicit none

      class(base_case_t) :: self
    end subroutine boundary_conditions

    subroutine initial_conditions(self)
      !! Sets case-specific initial conditions
      import :: base_case_t
      implicit none

      class(base_case_t) :: self
    end subroutine initial_conditions

    subroutine post_transeq(self, du, dv, dw)
      !! Applies case-specific or model realated forcings after transeq
      import :: base_case_t
      import :: field_t
      implicit none

      class(base_case_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
    end subroutine post_transeq

    subroutine postprocess(self, t)
      !! Triggers case-specific postprocessings at user specified intervals
      import :: base_case_t
      import :: dp
      implicit none

      class(base_case_t) :: self
      real(dp), intent(in) :: t
    end subroutine postprocess
  end interface

contains

  subroutine case_init(self, backend, mesh, host_allocator)
    implicit none

    class(base_case_t) :: self
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator

    self%solver = init(backend, mesh, host_allocator)

    call self%initial_conditions()

  end subroutine case_init

  subroutine print_enstrophy(self, u, v, w)
    !! Reports the enstrophy
    implicit none

    class(base_case_t), intent(in) :: self
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: du, dv, dw
    real(dp) :: enstrophy

    du => self%solver%backend%allocator%get_block(DIR_X, VERT)
    dv => self%solver%backend%allocator%get_block(DIR_X, VERT)
    dw => self%solver%backend%allocator%get_block(DIR_X, VERT)

    call self%solver%curl(du, dv, dw, u, v, w)
    enstrophy = 0.5_dp*(self%solver%backend%scalar_product(du, du) &
                        + self%solver%backend%scalar_product(dv, dv) &
                        + self%solver%backend%scalar_product(dw, dw)) &
                /self%solver%ngrid
    if (self%solver%mesh%par%is_root()) print *, 'enstrophy:', enstrophy

    call self%solver%backend%allocator%release_block(du)
    call self%solver%backend%allocator%release_block(dv)
    call self%solver%backend%allocator%release_block(dw)

  end subroutine print_enstrophy

  subroutine print_div_max_mean(self, u, v, w)
    !! Reports the div(u) at cell centres
    implicit none

    class(base_case_t), intent(in) :: self
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: div_u
    class(field_t), pointer :: u_out
    real(dp) :: div_u_max, div_u_mean
    integer :: ierr

    div_u => self%solver%backend%allocator%get_block(DIR_Z)

    call self%solver%divergence_v2p(div_u, u, v, w)

    u_out => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(u_out%data, div_u)

    call self%solver%backend%allocator%release_block(div_u)

    div_u_max = maxval(abs(u_out%data))
    div_u_mean = sum(abs(u_out%data))/self%solver%ngrid

    call self%solver%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%solver%mesh%par%is_root()) &
      print *, 'div u max mean:', div_u_max, div_u_mean

  end subroutine print_div_max_mean

  subroutine run(self)
    !! Runs the solver forwards in time from t=t_0 to t=T, performing
    !! postprocessing/IO and reporting diagnostics.
    implicit none

    class(base_case_t), intent(inout) :: self

    class(field_t), pointer :: du, dv, dw
    class(field_t), pointer :: u_out, v_out, w_out

    real(dp) :: t
    integer :: i, j

    if (self%solver%mesh%par%is_root()) print *, 'initial conditions'
    t = 0._dp
    call self%postprocess(t)

    if (self%solver%mesh%par%is_root()) print *, 'start run'

    do i = 1, self%solver%n_iters
      do j = 1, self%solver%time_integrator%nstage
        ! first apply case-specific BCs
        call self%boundary_conditions()

        du => self%solver%backend%allocator%get_block(DIR_X)
        dv => self%solver%backend%allocator%get_block(DIR_X)
        dw => self%solver%backend%allocator%get_block(DIR_X)

        call self%solver%transeq(du, dv, dw, &
                                 self%solver%u, self%solver%v, self%solver%w)

        ! models that introduce source terms handled here
        call self%post_transeq(du, dv, dw)

        ! time integration
        call self%solver%time_integrator%step( &
          self%solver%u, self%solver%v, self%solver%w, du, dv, dw, &
          self%solver%dt &
          )

        call self%solver%backend%allocator%release_block(du)
        call self%solver%backend%allocator%release_block(dv)
        call self%solver%backend%allocator%release_block(dw)

        call self%solver%pressure_correction(self%solver%u, self%solver%v, &
                                             self%solver%w)
      end do

      if (mod(i, self%solver%n_output) == 0) then
        t = i*self%solver%dt
        call self%postprocess(t)
      end if
    end do

    if (self%solver%mesh%par%is_root()) print *, 'run end'

    ! Below is for demonstrating purpuses only, to be removed when we have
    ! proper I/O in place.
    u_out => self%solver%host_allocator%get_block(DIR_C)
    v_out => self%solver%host_allocator%get_block(DIR_C)
    w_out => self%solver%host_allocator%get_block(DIR_C)

    call self%solver%backend%get_field_data(u_out%data, self%solver%u)
    call self%solver%backend%get_field_data(v_out%data, self%solver%v)
    call self%solver%backend%get_field_data(w_out%data, self%solver%w)

    if (self%solver%mesh%par%is_root()) then
      print *, 'norms', norm2(u_out%data), norm2(v_out%data), norm2(w_out%data)
    end if

    call self%solver%host_allocator%release_block(u_out)
    call self%solver%host_allocator%release_block(v_out)
    call self%solver%host_allocator%release_block(w_out)

  end subroutine run

end module m_base_case
