module m_solver
   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t
   use m_tdsops, only: tdsops_t, dirps_t
   use m_time_integrator, only: time_intg_t

   implicit none

   type :: solver_t
      !! solver class defines the Incompact3D algorithm at a very high level.
      !!
      !! Procedures defined here that are part of the Incompact3D algorithm
      !! are: transeq, divergence, poisson, and gradient.
      !!
      !! The operations these high level procedures require are provided by
      !! the relavant backend implementations.
      !!
      !! transeq procedure obtains the derivations in x, y, and z directions
      !! using the transeq_x, transeq_y, and transeq_z operations provided by
      !! the backend.
      !! There are two different algorithms available for this operation, a
      !! distributed algorithm and the Thomas algorithm. At the solver class
      !! level it isn't known which algorithm will be executed, that is decided
      !! at run time and therefore backend implementations are responsible for
      !! executing the right subroutines.
      !!
      !! Allocator is responsible from giving us a field sized array when
      !! requested. For example, when the derivations in x direction are
      !! completed and we are ready for the y directional derivatives, we need
      !! three fields to reorder and store the velocities in y direction. Also,
      !! we need three more fields for storing the results, and the get_block
      !! method of the allocator is used to arrange all these memory
      !! assignments. Later, when a field is no more required, release_block
      !! method of the allocator can be used to make this field available
      !! for later use.

      real(dp) :: dt, nu

      class(field_t), pointer :: u, v, w, du, dv, dw

      class(base_backend_t), pointer :: backend
      class(dirps_t), pointer :: xdirps, ydirps, zdirps
      class(time_intg_t), pointer :: time_integrator
   contains
      procedure :: transeq
      procedure :: run
   end type solver_t

   interface solver_t
      module procedure init
   end interface solver_t
contains

   function init(backend, time_integrator, xdirps, ydirps, zdirps, globs) &
      result(solver)
      implicit none

      class(base_backend_t), target, intent(inout) :: backend
      class(time_intg_t), target, intent(inout) :: time_integrator
      class(dirps_t), target, intent(inout) :: xdirps, ydirps, zdirps
      class(globs_t), intent(in) :: globs
      type(solver_t) :: solver

      real(dp), allocatable, dimension(:, :, :) :: u_init, v_init, w_init
      integer :: dims(3)

      real(dp) :: dx, dy, dz
      integer :: nx, ny, nz

      solver%backend => backend
      solver%time_integrator => time_integrator

      solver%xdirps => xdirps
      solver%ydirps => ydirps
      solver%zdirps => zdirps

      solver%u => solver%backend%allocator%get_block()
      solver%v => solver%backend%allocator%get_block()
      solver%w => solver%backend%allocator%get_block()

      ! Set initial conditions
      dims(:) = solver%backend%allocator%dims(:)
      allocate(u_init(dims(1), dims(2), dims(3)))
      allocate(v_init(dims(1), dims(2), dims(3)))
      allocate(w_init(dims(1), dims(2), dims(3)))

      u_init = 0
      v_init = 0
      w_init = 0

      call solver%backend%set_fields( &
         solver%u, solver%v, solver%w, u_init, v_init, w_init &
      )

      deallocate(u_init, v_init, w_init)
      print*, 'initial conditions are set'

      ! Allocate fields for storing the RHS
      solver%du => solver%backend%allocator%get_block()
      solver%dv => solver%backend%allocator%get_block()
      solver%dw => solver%backend%allocator%get_block()

      nx = globs%nx_loc; ny = globs%ny_loc; nz = globs%nz_loc
      dx = globs%dx; dy = globs%dy; dz = globs%dz

      ! Allocate and set the tdsops
      call solver%backend%alloc_tdsops(solver%xdirps%der1st, nx, dx, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%ydirps%der1st, ny, dy, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%zdirps%der1st, nz, dz, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%xdirps%der1st_sym, nx, dx, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%ydirps%der1st_sym, ny, dy, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%zdirps%der1st_sym, nz, dz, &
                                       'first-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%xdirps%der2nd, nx, dx, &
                                       'second-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%ydirps%der2nd, ny, dy, &
                                       'second-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%zdirps%der2nd, nz, dz, &
                                       'second-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%xdirps%der2nd_sym, nx, dx, &
                                       'second-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%ydirps%der2nd_sym, ny, dy, &
                                       'second-deriv', 'compact6')
      call solver%backend%alloc_tdsops(solver%zdirps%der2nd_sym, nz, dz, &
                                       'second-deriv', 'compact6')

   end function init

   subroutine transeq(self, du, dv, dw, u, v, w)
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w

      class(field_t), pointer :: u_y, v_y, w_y, u_z, v_z, w_z, &
                                 du_y, dv_y, dw_y, du_z, dv_z, dw_z

      ! -1/2(nabla u curl u + u nabla u) + nu nablasq u

      ! call derivatives in x direction. Based on the run time arguments this
      ! executes a distributed algorithm or the Thomas algorithm.
      call self%backend%transeq_x(du, dv, dw, u, v, w, self%xdirps)

      ! request fields from the allocator
      u_y => self%backend%allocator%get_block()
      v_y => self%backend%allocator%get_block()
      w_y => self%backend%allocator%get_block()
      du_y => self%backend%allocator%get_block()
      dv_y => self%backend%allocator%get_block()
      dw_y => self%backend%allocator%get_block()

      ! reorder data from x orientation to y orientation
      call self%backend%trans_x2y(u_y, v_y, w_y, u, v, w)
      ! similar to the x direction, obtain derivatives in y.
      call self%backend%transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, self%ydirps)

      ! we don't need the velocities in y orientation any more, so release
      ! them to open up space.
      ! It is important that this doesn't actually deallocate any memory,
      ! it just makes the corresponding memory space available for use.
      call self%backend%allocator%release_block(u_y)
      call self%backend%allocator%release_block(v_y)
      call self%backend%allocator%release_block(w_y)

      ! just like in y direction, get some fields for the z derivatives.
      u_z => self%backend%allocator%get_block()
      v_z => self%backend%allocator%get_block()
      w_z => self%backend%allocator%get_block()
      du_z => self%backend%allocator%get_block()
      dv_z => self%backend%allocator%get_block()
      dw_z => self%backend%allocator%get_block()

      ! reorder from x to z
      call self%backend%trans_x2z(u_z, v_z, w_z, u, v, w)
      ! get the derivatives in z
      call self%backend%transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, self%zdirps)

      ! there is no need to keep velocities in z orientation around, so release
      call self%backend%allocator%release_block(u_z)
      call self%backend%allocator%release_block(v_z)
      call self%backend%allocator%release_block(w_z)

      ! gather all the contributions into the x result array
      ! this function does the data reordering and summations at once.
      call self%backend%sum_yzintox(du, dv, dw, &
                                    du_y, dv_y, dw_y, du_z, dv_z, dw_z)

      ! release all the unnecessary blocks.
      call self%backend%allocator%release_block(du_y)
      call self%backend%allocator%release_block(dv_y)
      call self%backend%allocator%release_block(dw_y)
      call self%backend%allocator%release_block(du_z)
      call self%backend%allocator%release_block(dv_z)
      call self%backend%allocator%release_block(dw_z)

   end subroutine transeq

   subroutine run(self, n_iter, u_out, v_out, w_out)
      implicit none

      class(solver_t), intent(in) :: self
      integer, intent(in) :: n_iter
      real(dp), dimension(:, :, :), intent(inout) :: u_out, v_out, w_out

      integer :: i

      print*, 'start run'

      do i = 1, n_iter
         call self%transeq(self%du, self%dv, self%dw, self%u, self%v, self%w)

         ! time integration
         call self%time_integrator%step(self%u, self%v, self%w, &
                                        self%du, self%dv, self%dw, self%dt)

         !! pressure stuff
         !call self%divergence(udiv, u, v, w)
         !call self%poisson(p, udiv)
         !call self%gradient(px, py, pz, p)

         !! velocity correction
         !call self%backend%vec3add(u, v, w, px, py, pz)
      end do

      print*, 'run end'

      call self%backend%get_fields( &
         u_out, v_out, w_out, self%du, self%dv, self%dw &
      )

   end subroutine run

end module m_solver
