module m_solver
   use m_allocator, only: allocator_t, field_t
   use m_base_backend, only: base_backend_t
   use m_common, only: dp, globs_t, RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2Y
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

      class(field_t), pointer :: u, v, w

      class(base_backend_t), pointer :: backend
      class(dirps_t), pointer :: xdirps, ydirps, zdirps
      class(time_intg_t), pointer :: time_integrator
   contains
      procedure :: transeq
      procedure :: divergence_v2p
      procedure :: gradient_p2v
      procedure :: curl
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

      nx = globs%nx_loc; ny = globs%ny_loc; nz = globs%nz_loc
      dx = globs%dx; dy = globs%dy; dz = globs%dz

      ! Allocate and set the tdsops
      call allocate_tdsops(solver%xdirps, nx, dx, solver%backend)
      call allocate_tdsops(solver%ydirps, ny, dy, solver%backend)
      call allocate_tdsops(solver%zdirps, nz, dz, solver%backend)

   end function init

   subroutine allocate_tdsops(dirps, nx, dx, backend)
      class(dirps_t), intent(inout) :: dirps
      real(dp), intent(in) :: dx
      integer, intent(in) :: nx
      class(base_backend_t), intent(in) :: backend

      call backend%alloc_tdsops(dirps%der1st, nx, dx, &
                                'first-deriv', 'compact6')
      call backend%alloc_tdsops(dirps%der1st_sym, nx, dx, &
                                'first-deriv', 'compact6')
      call backend%alloc_tdsops(dirps%der2nd, nx, dx, &
                                'second-deriv', 'compact6')
      call backend%alloc_tdsops(dirps%der2nd_sym, nx, dx, &
                                'second-deriv', 'compact6')
      call backend%alloc_tdsops(dirps%interpl_v2p, nx, dx, &
                                'interpolate', 'classic', from_to='v2p')
      call backend%alloc_tdsops(dirps%interpl_p2v, nx, dx, &
                                'interpolate', 'classic', from_to='p2v')
      call backend%alloc_tdsops(dirps%stagder_v2p, nx, dx, &
                                'stag-deriv', 'compact6', from_to='v2p')
      call backend%alloc_tdsops(dirps%stagder_p2v, nx, dx, &
                                'stag-deriv', 'compact6', from_to='p2v')


   end subroutine

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
      call self%backend%reorder(u_y, u, RDR_X2Y)
      call self%backend%reorder(v_y, v, RDR_X2Y)
      call self%backend%reorder(w_y, w, RDR_X2Y)

      ! similar to the x direction, obtain derivatives in y.
      call self%backend%transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, self%ydirps)

      ! we don't need the velocities in y orientation any more, so release
      ! them to open up space.
      ! It is important that this doesn't actually deallocate any memory,
      ! it just makes the corresponding memory space available for use.
      call self%backend%allocator%release_block(u_y)
      call self%backend%allocator%release_block(v_y)
      call self%backend%allocator%release_block(w_y)

      call self%backend%sum_yintox(du, du_y)
      call self%backend%sum_yintox(dv, dv_y)
      call self%backend%sum_yintox(dw, dw_y)

      call self%backend%allocator%release_block(du_y)
      call self%backend%allocator%release_block(dv_y)
      call self%backend%allocator%release_block(dw_y)

      ! just like in y direction, get some fields for the z derivatives.
      u_z => self%backend%allocator%get_block()
      v_z => self%backend%allocator%get_block()
      w_z => self%backend%allocator%get_block()
      du_z => self%backend%allocator%get_block()
      dv_z => self%backend%allocator%get_block()
      dw_z => self%backend%allocator%get_block()

      ! reorder from x to z
      call self%backend%reorder(u_z, u, RDR_X2Z)
      call self%backend%reorder(v_z, v, RDR_X2Z)
      call self%backend%reorder(w_z, w, RDR_X2Z)

      ! get the derivatives in z
      call self%backend%transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, self%zdirps)

      ! there is no need to keep velocities in z orientation around, so release
      call self%backend%allocator%release_block(u_z)
      call self%backend%allocator%release_block(v_z)
      call self%backend%allocator%release_block(w_z)

      ! gather all the contributions into the x result array
      call self%backend%sum_zintox(du, du_z)
      call self%backend%sum_zintox(dv, dv_z)
      call self%backend%sum_zintox(dw, dw_z)

      ! release all the unnecessary blocks.
      call self%backend%allocator%release_block(du_z)
      call self%backend%allocator%release_block(dv_z)
      call self%backend%allocator%release_block(dw_z)

   end subroutine transeq

   subroutine divergence_v2p(self, div_u, u, v, w)
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: div_u
      class(field_t), intent(in) :: u, v, w

      class(field_t), pointer :: du_x, dv_x, dw_x, &
                                 u_y, v_y, w_y, du_y, dv_y, dw_y, &
                                 u_z, w_z, dw_z

      du_x => self%backend%allocator%get_block()
      dv_x => self%backend%allocator%get_block()
      dw_x => self%backend%allocator%get_block()

      ! Staggared der for u field in x
      ! Interpolation for v field in x
      ! Interpolation for w field in x
      call self%backend%tds_solve(du_x, u, self%xdirps, &
                                  self%xdirps%stagder_v2p)
      call self%backend%tds_solve(dv_x, v, self%xdirps, &
                                  self%xdirps%interpl_v2p)
      call self%backend%tds_solve(dw_x, w, self%xdirps, &
                                  self%xdirps%interpl_v2p)

      ! request fields from the allocator
      u_y => self%backend%allocator%get_block()
      v_y => self%backend%allocator%get_block()
      w_y => self%backend%allocator%get_block()

      ! reorder data from x orientation to y orientation
      call self%backend%reorder(u_y, du_x, RDR_X2Y)
      call self%backend%reorder(v_y, dv_x, RDR_X2Y)
      call self%backend%reorder(w_y, dw_x, RDR_X2Y)

      call self%backend%allocator%release_block(du_x)
      call self%backend%allocator%release_block(dv_x)
      call self%backend%allocator%release_block(dw_x)

      du_y => self%backend%allocator%get_block()
      dv_y => self%backend%allocator%get_block()
      dw_y => self%backend%allocator%get_block()

      ! similar to the x direction, obtain derivatives in y.
      call self%backend%tds_solve(du_y, u_y, self%ydirps, &
                                  self%ydirps%interpl_v2p)
      call self%backend%tds_solve(dv_y, v_y, self%ydirps, &
                                  self%ydirps%stagder_v2p)
      call self%backend%tds_solve(dw_y, w_y, self%ydirps, &
                                  self%ydirps%interpl_v2p)

      ! we don't need the velocities in y orientation any more, so release
      ! them to open up space.
      ! It is important that this doesn't actually deallocate any memory,
      ! it just makes the corresponding memory space available for use.
      call self%backend%allocator%release_block(u_y)
      call self%backend%allocator%release_block(v_y)
      call self%backend%allocator%release_block(w_y)

      ! just like in y direction, get some fields for the z derivatives.
      u_z => self%backend%allocator%get_block()
      w_z => self%backend%allocator%get_block()

      ! dv_y = dv_y + dw_y
      call self%backend%vecadd(1._dp, dw_y, 1._dp, dv_y)

      ! reorder from y to z
      call self%backend%reorder(u_z, du_y, RDR_Y2Z)
      call self%backend%reorder(w_z, dw_y, RDR_Y2Z)

      ! release all the unnecessary blocks.
      call self%backend%allocator%release_block(du_y)
      call self%backend%allocator%release_block(dv_y)
      call self%backend%allocator%release_block(dw_y)

      dw_z => self%backend%allocator%get_block()

      ! get the derivatives in z
      call self%backend%tds_solve(div_u, u_z, self%zdirps, &
                                  self%zdirps%interpl_v2p)
      call self%backend%tds_solve(dw_z, w_z, self%zdirps, &
                                  self%zdirps%stagder_v2p)

      ! div_u = div_u + dw_z
      call self%backend%vecadd(1._dp, dw_z, 1._dp, div_u)

      ! div_u array is in z orientation

      ! there is no need to keep velocities in z orientation around, so release
      call self%backend%allocator%release_block(u_z)
      call self%backend%allocator%release_block(w_z)
      call self%backend%allocator%release_block(dw_z)

   end subroutine divergence_v2p

   subroutine gradient_p2v(self, dpdx, dpdy, dpdz, pressure)
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: dpdx, dpdy, dpdz
      class(field_t), intent(in) :: pressure

      class(field_t), pointer :: p_sxy_z, dpdz_sxy_z, &
                                 p_sxy_y, dpdz_sxy_y, &
                                 p_sx_y, dpdy_sx_y, dpdz_sx_y, &
                                 p_sx_x, dpdy_sx_x, dpdz_sx_x

      p_sxy_z => self%backend%allocator%get_block()
      dpdz_sxy_z => self%backend%allocator%get_block()

      ! Staggared der for pressure field in z
      ! Interpolation for pressure field in z
      call self%backend%tds_solve(p_sxy_z, pressure, self%zdirps, &
                                  self%zdirps%interpl_p2v)
      call self%backend%tds_solve(dpdz_sxy_z, pressure, self%zdirps, &
                                  self%zdirps%stagder_p2v)

      ! request fields from the allocator
      p_sxy_y => self%backend%allocator%get_block()
      dpdz_sxy_y => self%backend%allocator%get_block()

      ! reorder data from z orientation to y orientation
      call self%backend%reorder(p_sxy_y, p_sxy_z, RDR_Z2Y)
      call self%backend%reorder(dpdz_sxy_y, dpdz_sxy_z, RDR_Z2Y)

      call self%backend%allocator%release_block(p_sxy_z)
      call self%backend%allocator%release_block(dpdz_sxy_z)

      p_sx_y => self%backend%allocator%get_block()
      dpdy_sx_y => self%backend%allocator%get_block()
      dpdz_sx_y => self%backend%allocator%get_block()

      ! similar to the z direction, obtain derivatives in y.
      call self%backend%tds_solve(p_sx_y, p_sxy_y, self%ydirps, &
                                  self%ydirps%interpl_p2v)
      call self%backend%tds_solve(dpdy_sx_y, p_sxy_y, self%ydirps, &
                                  self%ydirps%stagder_p2v)
      call self%backend%tds_solve(dpdz_sx_y, dpdz_sxy_y, self%ydirps, &
                                  self%ydirps%interpl_p2v)

      ! release memory
      call self%backend%allocator%release_block(p_sxy_y)
      call self%backend%allocator%release_block(dpdz_sxy_y)

      ! just like in y direction, get some fields for the x derivatives.
      p_sx_x => self%backend%allocator%get_block()
      dpdy_sx_x => self%backend%allocator%get_block()
      dpdz_sx_x => self%backend%allocator%get_block()

      ! reorder from y to x
      call self%backend%reorder(p_sx_x, p_sx_y, RDR_Y2X)
      call self%backend%reorder(dpdy_sx_x, dpdy_sx_y, RDR_Y2X)
      call self%backend%reorder(dpdz_sx_x, dpdz_sx_y, RDR_Y2X)

      ! release all the y directional fields.
      call self%backend%allocator%release_block(p_sx_y)
      call self%backend%allocator%release_block(dpdy_sx_y)
      call self%backend%allocator%release_block(dpdz_sx_y)

      ! get the derivatives in x
      call self%backend%tds_solve(dpdx, p_sx_x, self%xdirps, &
                                  self%xdirps%stagder_p2v)
      call self%backend%tds_solve(dpdy, dpdy_sx_x, self%xdirps, &
                                  self%xdirps%interpl_v2p)
      call self%backend%tds_solve(dpdz, dpdz_sx_x, self%xdirps, &
                                  self%xdirps%interpl_v2p)

      ! release temporary x fields
      call self%backend%allocator%release_block(p_sx_x)
      call self%backend%allocator%release_block(dpdy_sx_x)
      call self%backend%allocator%release_block(dpdz_sx_x)

   end subroutine gradient_p2v

   subroutine curl(self, o_x, o_y, o_z, u, v, w)
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: o_x, o_y, o_z !! omega_x/_y/_z
      class(field_t), intent(in) :: u, v, w

      class(field_t), pointer :: u_y, w_y, du_y, u_z, v_z, du_z, dv_z

      ! o_x = dw/dy - dv/dz
      ! o_y = du/dz - dw/dx
      ! o_z = dv/dx - du/dy

      ! obtain dw/dx, dv/dx and store them directly in omega_y, omega_z
      call self%backend%tds_solve(o_y, w, self%xdirps, self%xdirps%der1st)
      call self%backend%tds_solve(o_z, v, self%xdirps, self%xdirps%der1st)

      u_y => self%backend%allocator%get_block()
      w_y => self%backend%allocator%get_block()

      call self%backend%reorder(u_y, u, RDR_X2Y)
      call self%backend%reorder(w_y, w, RDR_X2Y)

      du_y => self%backend%allocator%get_block()

      ! obtain du/dy, dw/dy
      ! store du/dy in a temporary field to add into omega_z later
      ! dw/dy can be stored directly in omega_x as it is empty
      call self%backend%tds_solve(du_y, u_y, self%ydirps, self%ydirps%der1st)
      call self%backend%tds_solve(o_x, w_y, self%ydirps, self%ydirps%der1st)

      call self%backend%allocator%release_block(u_y)
      call self%backend%allocator%release_block(w_y)

      ! omega_z = dv/dz - du/dy
      call self%backend%vecadd(-1._dp, du_y, 1._dp, o_z)

      call self%backend%allocator%release_block(du_y)

      u_z => self%backend%allocator%get_block()
      v_z => self%backend%allocator%get_block()
      du_z => self%backend%allocator%get_block()
      dv_z => self%backend%allocator%get_block()

      ! obtain du/dz, dv/dz and store them in temporary fields
      call self%backend%tds_solve(du_z, u_z, self%zdirps, self%zdirps%der1st)
      call self%backend%tds_solve(dv_z, v_z, self%zdirps, self%zdirps%der1st)

      ! omega_x = dw/dy - dv/dz
      call self%backend%vecadd(-1._dp, dv_z, 1._dp, o_x)
      ! omega_y = du/dz - dw/dx
      call self%backend%vecadd(1._dp, du_z, -1._dp, o_y)

      call self%backend%allocator%release_block(u_z)
      call self%backend%allocator%release_block(v_z)
      call self%backend%allocator%release_block(du_z)
      call self%backend%allocator%release_block(dv_z)

   end subroutine curl

   subroutine run(self, n_iter, u_out, v_out, w_out)
      implicit none

      class(solver_t), intent(in) :: self
      integer, intent(in) :: n_iter
      real(dp), dimension(:, :, :), intent(inout) :: u_out, v_out, w_out

      class(field_t), pointer :: du, dv, dw, div_u, pressure, dpdx, dpdy, dpdz

      integer :: i

      print*, 'start run'

      do i = 1, n_iter
         du => self%backend%allocator%get_block()
         dv => self%backend%allocator%get_block()
         dw => self%backend%allocator%get_block()

         call self%transeq(du, dv, dw, self%u, self%v, self%w)

         ! time integration
         call self%time_integrator%step(self%u, self%v, self%w, &
                                        du, dv, dw, self%dt)

         call self%backend%allocator%release_block(du)
         call self%backend%allocator%release_block(dv)
         call self%backend%allocator%release_block(dw)

         ! pressure
         div_u => self%backend%allocator%get_block()

         call self%divergence_v2p(div_u, self%u, self%v, self%w)

         pressure => self%backend%allocator%get_block()

         !call self%poisson(pressure, div_u)

         call self%backend%allocator%release_block(div_u)

         dpdx => self%backend%allocator%get_block()
         dpdy => self%backend%allocator%get_block()
         dpdz => self%backend%allocator%get_block()

         call self%gradient_p2v(dpdx, dpdy, dpdz, pressure)

         call self%backend%allocator%release_block(pressure)

         ! velocity correction
         call self%backend%vecadd(-1._dp, dpdx, 1._dp, self%u)
         call self%backend%vecadd(-1._dp, dpdy, 1._dp, self%v)
         call self%backend%vecadd(-1._dp, dpdz, 1._dp, self%w)

         call self%backend%allocator%release_block(dpdx)
         call self%backend%allocator%release_block(dpdy)
         call self%backend%allocator%release_block(dpdz)
      end do

      print*, 'run end'

      call self%backend%get_fields( &
         u_out, v_out, w_out, self%u, self%v, self%w &
      )

   end subroutine run

end module m_solver
