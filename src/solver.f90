module m_solver
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      RDR_Z2C, RDR_C2Z, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, CELL
  use m_tdsops, only: tdsops_t, dirps_t
  use m_time_integrator, only: time_intg_t
  use m_mesh, only: mesh_t

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
    integer :: n_iters, n_output
    integer :: ngrid

    class(field_t), pointer :: u, v, w

    class(base_backend_t), pointer :: backend
    class(mesh_t), pointer :: mesh
    type(time_intg_t) :: time_integrator
    type(allocator_t), pointer :: host_allocator
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    procedure(poisson_solver), pointer :: poisson => null()
  contains
    procedure :: transeq
    procedure :: divergence_v2p
    procedure :: gradient_p2v
    procedure :: curl
    procedure :: output
    procedure :: run
  end type solver_t

  abstract interface
    subroutine poisson_solver(self, pressure, div_u)
      import :: solver_t
      import :: field_t
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: pressure
      class(field_t), intent(in) :: div_u
    end subroutine poisson_solver
  end interface

  interface solver_t
    module procedure init
  end interface solver_t

contains

  function init(backend, mesh, host_allocator) result(solver)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(solver_t) :: solver

    class(field_t), pointer :: u_init, v_init, w_init

    character(len=200) :: input_file
    real(dp) :: Re, dt
    integer :: n_iters, n_output
    character(3) :: poisson_solver_type, time_intg
    character(30) :: der1st_scheme, der2nd_scheme, &
                     interpl_scheme, stagder_scheme
    namelist /solver_params/ Re, dt, n_iters, n_output, poisson_solver_type, &
      time_intg, der1st_scheme, der2nd_scheme, &
      interpl_scheme, stagder_scheme

    real(dp) :: x, y, z
    integer :: i, j, k
    integer, dimension(3) :: dims
    real(dp), dimension(3) :: xloc

    solver%backend => backend
    solver%mesh => mesh
    solver%host_allocator => host_allocator

    allocate (solver%xdirps, solver%ydirps, solver%zdirps)
    solver%xdirps%dir = DIR_X
    solver%ydirps%dir = DIR_Y
    solver%zdirps%dir = DIR_Z

    solver%u => solver%backend%allocator%get_block(DIR_X, VERT)
    solver%v => solver%backend%allocator%get_block(DIR_X, VERT)
    solver%w => solver%backend%allocator%get_block(DIR_X, VERT)

    ! set defaults
    poisson_solver_type = 'FFT'
    time_intg = 'AB3'
    der1st_scheme = 'compact6'; der2nd_scheme = 'compact6'
    interpl_scheme = 'classic'; stagder_scheme = 'compact6'

    if (command_argument_count() >= 1) then
      call get_command_argument(1, input_file)
      open (100, file=input_file)
      read (100, nml=solver_params)
      close (100)
    else
      error stop 'Input file is not provided.'
    end if

    solver%time_integrator = time_intg_t(solver%backend, &
                                         solver%backend%allocator, &
                                         time_intg)
    if (solver%mesh%par%is_root()) then
      print *, time_intg//' time integrator instantiated'
    end if

    solver%dt = dt
    solver%backend%nu = 1._dp/Re
    solver%n_iters = n_iters
    solver%n_output = n_output
    solver%ngrid = product(solver%mesh%get_global_dims(VERT))

    dims = solver%mesh%get_dims(VERT)
    u_init => solver%host_allocator%get_block(DIR_C)
    v_init => solver%host_allocator%get_block(DIR_C)
    w_init => solver%host_allocator%get_block(DIR_C)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = solver%mesh%get_coordinates(i, j, k)
          x = xloc(1)
          y = xloc(2)
          z = xloc(3)

          u_init%data(i, j, k) = sin(x)*cos(y)*cos(z)
          v_init%data(i, j, k) = -cos(x)*sin(y)*cos(z)
          w_init%data(i, j, k) = 0
        end do
      end do
    end do

    call solver%backend%set_field_data(solver%u, u_init%data)
    call solver%backend%set_field_data(solver%v, v_init%data)
    call solver%backend%set_field_data(solver%w, w_init%data)

    call solver%host_allocator%release_block(u_init)
    call solver%host_allocator%release_block(v_init)
    call solver%host_allocator%release_block(w_init)

    ! Allocate and set the tdsops
    call allocate_tdsops(solver%xdirps, solver%backend, &
                         der1st_scheme, der2nd_scheme, &
                         interpl_scheme, stagder_scheme)
    call allocate_tdsops(solver%ydirps, solver%backend, &
                         der1st_scheme, der2nd_scheme, &
                         interpl_scheme, stagder_scheme)
    call allocate_tdsops(solver%zdirps, solver%backend, &
                         der1st_scheme, der2nd_scheme, &
                         interpl_scheme, stagder_scheme)

    select case (trim(poisson_solver_type))
    case ('FFT')
      if (solver%mesh%par%is_root()) print *, 'Poisson solver: FFT'
      call solver%backend%init_poisson_fft(solver%mesh, solver%xdirps, &
                                           solver%ydirps, solver%zdirps)
      solver%poisson => poisson_fft
    case ('CG')
      if (solver%mesh%par%is_root()) &
        print *, 'Poisson solver: CG, not yet implemented'
      solver%poisson => poisson_cg
    case default
      error stop 'poisson_solver_type is not valid. Use "FFT" or "CG".'
    end select

  end function init

  subroutine allocate_tdsops(dirps, backend, der1st_scheme, der2nd_scheme, &
                             interpl_scheme, stagder_scheme)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    character(*), intent(in) :: der1st_scheme, der2nd_scheme, &
                                interpl_scheme, stagder_scheme

    call backend%alloc_tdsops(dirps%der1st, dirps%dir, &
                              'first-deriv', der1st_scheme)
    call backend%alloc_tdsops(dirps%der1st_sym, dirps%dir, &
                              'first-deriv', der1st_scheme)
    call backend%alloc_tdsops(dirps%der2nd, dirps%dir, &
                              'second-deriv', der2nd_scheme)
    call backend%alloc_tdsops(dirps%der2nd_sym, dirps%dir, &
                              'second-deriv', der2nd_scheme)
    call backend%alloc_tdsops(dirps%interpl_v2p, dirps%dir, &
                              'interpolate', interpl_scheme, from_to='v2p')
    call backend%alloc_tdsops(dirps%interpl_p2v, dirps%dir, &
                              'interpolate', interpl_scheme, from_to='p2v')
    call backend%alloc_tdsops(dirps%stagder_v2p, dirps%dir, &
                              'stag-deriv', stagder_scheme, from_to='v2p')
    call backend%alloc_tdsops(dirps%stagder_p2v, dirps%dir, &
                              'stag-deriv', stagder_scheme, from_to='p2v')

  end subroutine

  subroutine transeq(self, du, dv, dw, u, v, w)
      !! Skew-symmetric form of convection-diffusion terms in the
      !! incompressible Navier-Stokes momemtum equations, excluding
      !! pressure terms.
      !! Inputs from velocity grid and outputs to velocity grid.
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
    u_y => self%backend%allocator%get_block(DIR_Y, VERT)
    v_y => self%backend%allocator%get_block(DIR_Y, VERT)
    w_y => self%backend%allocator%get_block(DIR_Y, VERT)
    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

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
    u_z => self%backend%allocator%get_block(DIR_Z, VERT)
    v_z => self%backend%allocator%get_block(DIR_Z, VERT)
    w_z => self%backend%allocator%get_block(DIR_Z, VERT)
    du_z => self%backend%allocator%get_block(DIR_Z)
    dv_z => self%backend%allocator%get_block(DIR_Z)
    dw_z => self%backend%allocator%get_block(DIR_Z)

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
      !! Divergence of a vector field (u, v, w).
      !! Inputs from velocity grid and outputs to pressure grid.
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: div_u
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: du_x, dv_x, dw_x, &
      u_y, v_y, w_y, du_y, dv_y, dw_y, &
      u_z, w_z, dw_z

    du_x => self%backend%allocator%get_block(DIR_X)
    dv_x => self%backend%allocator%get_block(DIR_X)
    dw_x => self%backend%allocator%get_block(DIR_X)

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
    u_y => self%backend%allocator%get_block(DIR_Y, VERT)
    v_y => self%backend%allocator%get_block(DIR_Y, VERT)
    w_y => self%backend%allocator%get_block(DIR_Y, VERT)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, du_x, RDR_X2Y)
    call self%backend%reorder(v_y, dv_x, RDR_X2Y)
    call self%backend%reorder(w_y, dw_x, RDR_X2Y)

    call self%backend%allocator%release_block(du_x)
    call self%backend%allocator%release_block(dv_x)
    call self%backend%allocator%release_block(dw_x)

    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

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
    u_z => self%backend%allocator%get_block(DIR_Z, VERT)
    w_z => self%backend%allocator%get_block(DIR_Z, VERT)

    ! du_y = dv_y + du_y
    call self%backend%vecadd(1._dp, dv_y, 1._dp, du_y)

    ! reorder from y to z
    call self%backend%reorder(u_z, du_y, RDR_Y2Z)
    call self%backend%reorder(w_z, dw_y, RDR_Y2Z)

    ! release all the unnecessary blocks.
    call self%backend%allocator%release_block(du_y)
    call self%backend%allocator%release_block(dv_y)
    call self%backend%allocator%release_block(dw_y)

    dw_z => self%backend%allocator%get_block(DIR_Z)

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
      !! Gradient of a scalar field 'pressure'.
      !! Inputs from pressure grid and outputs to velocity grid.
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: dpdx, dpdy, dpdz
    class(field_t), intent(in) :: pressure

    class(field_t), pointer :: p_sxy_z, dpdz_sxy_z, &
      p_sxy_y, dpdz_sxy_y, &
      p_sx_y, dpdy_sx_y, dpdz_sx_y, &
      p_sx_x, dpdy_sx_x, dpdz_sx_x

    p_sxy_z => self%backend%allocator%get_block(DIR_Z)
    dpdz_sxy_z => self%backend%allocator%get_block(DIR_Z)

    ! Staggared der for pressure field in z
    ! Interpolation for pressure field in z
    call self%backend%tds_solve(p_sxy_z, pressure, self%zdirps, &
                                self%zdirps%interpl_p2v)
    call self%backend%tds_solve(dpdz_sxy_z, pressure, self%zdirps, &
                                self%zdirps%stagder_p2v)

    ! request fields from the allocator
    p_sxy_y => self%backend%allocator%get_block(DIR_Y)
    dpdz_sxy_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from z orientation to y orientation
    call self%backend%reorder(p_sxy_y, p_sxy_z, RDR_Z2Y)
    call self%backend%reorder(dpdz_sxy_y, dpdz_sxy_z, RDR_Z2Y)

    call self%backend%allocator%release_block(p_sxy_z)
    call self%backend%allocator%release_block(dpdz_sxy_z)

    p_sx_y => self%backend%allocator%get_block(DIR_Y)
    dpdy_sx_y => self%backend%allocator%get_block(DIR_Y)
    dpdz_sx_y => self%backend%allocator%get_block(DIR_Y)

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
    p_sx_x => self%backend%allocator%get_block(DIR_X)
    dpdy_sx_x => self%backend%allocator%get_block(DIR_X)
    dpdz_sx_x => self%backend%allocator%get_block(DIR_X)

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
                                self%xdirps%interpl_p2v)
    call self%backend%tds_solve(dpdz, dpdz_sx_x, self%xdirps, &
                                self%xdirps%interpl_p2v)

    ! release temporary x fields
    call self%backend%allocator%release_block(p_sx_x)
    call self%backend%allocator%release_block(dpdy_sx_x)
    call self%backend%allocator%release_block(dpdz_sx_x)

  end subroutine gradient_p2v

  subroutine curl(self, o_i_hat, o_j_hat, o_k_hat, u, v, w)
      !! Curl of a vector field (u, v, w).
      !! Inputs from velocity grid and outputs to velocity grid.
    implicit none

    class(solver_t) :: self
    !> Vector components of the output vector field Omega
    class(field_t), intent(inout) :: o_i_hat, o_j_hat, o_k_hat
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: u_y, u_z, v_z, w_y, dwdy_y, dvdz_z, dvdz_x, &
      dudz_z, dudz_x, dudy_y, dudy_x

    ! omega_i_hat = dw/dy - dv/dz
    ! omega_j_hat = du/dz - dw/dx
    ! omega_k_hat = dv/dx - du/dy

    ! omega_i_hat
    ! dw/dy
    w_y => self%backend%allocator%get_block(DIR_Y, VERT)
    dwdy_y => self%backend%allocator%get_block(DIR_Y)
    call self%backend%reorder(w_y, w, RDR_X2Y)
    call self%backend%tds_solve(dwdy_y, w_y, self%ydirps, self%ydirps%der1st)

    call self%backend%reorder(o_i_hat, dwdy_y, RDR_Y2X)

    call self%backend%allocator%release_block(w_y)
    call self%backend%allocator%release_block(dwdy_y)

    ! dv/dz
    v_z => self%backend%allocator%get_block(DIR_Z)
    dvdz_z => self%backend%allocator%get_block(DIR_Z)
    call self%backend%reorder(v_z, v, RDR_X2Z)
    call self%backend%tds_solve(dvdz_z, v_z, self%zdirps, self%zdirps%der1st)

    dvdz_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dvdz_x, dvdz_z, RDR_Z2X)

    call self%backend%allocator%release_block(v_z)
    call self%backend%allocator%release_block(dvdz_z)

    ! omega_i_hat = dw/dy - dv/dz
    call self%backend%vecadd(-1._dp, dvdz_x, 1._dp, o_i_hat)

    call self%backend%allocator%release_block(dvdz_x)

    ! omega_j_hat
    ! du/dz
    u_z => self%backend%allocator%get_block(DIR_Z, VERT)
    dudz_z => self%backend%allocator%get_block(DIR_Z)
    call self%backend%reorder(u_z, u, RDR_X2Z)
    call self%backend%tds_solve(dudz_z, u_z, self%zdirps, self%zdirps%der1st)

    dudz_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dudz_x, dudz_z, RDR_Z2X)

    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(dudz_z)

    ! dw/dx
    call self%backend%tds_solve(o_j_hat, w, self%xdirps, self%xdirps%der1st)

    ! omega_j_hat = du/dz - dw/dx
    call self%backend%vecadd(1._dp, dudz_x, -1._dp, o_j_hat)

    call self%backend%allocator%release_block(dudz_x)

    ! omega_k_hat
    ! dv/dx
    call self%backend%tds_solve(o_k_hat, v, self%xdirps, self%xdirps%der1st)

    ! du/dy
    u_y => self%backend%allocator%get_block(DIR_Y, VERT)
    dudy_y => self%backend%allocator%get_block(DIR_Y)
    call self%backend%reorder(u_y, u, RDR_X2Y)
    call self%backend%tds_solve(dudy_y, u_y, self%ydirps, self%ydirps%der1st)

    dudy_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dudy_x, dudy_y, RDR_Y2X)

    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(dudy_y)

    ! omega_k_hat = dv/dx - du/dy
    call self%backend%vecadd(-1._dp, dudy_x, 1._dp, o_k_hat)

    call self%backend%allocator%release_block(dudy_x)

  end subroutine curl

  subroutine poisson_fft(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

    class(field_t), pointer :: p_temp

    ! reorder into 3D Cartesian data structure
    p_temp => self%backend%allocator%get_block(DIR_C, CELL)
    call self%backend%reorder(p_temp, div_u, RDR_Z2C)

    ! call forward FFT
    ! output array in spectral space is stored at poisson_fft class
    call self%backend%poisson_fft%fft_forward(p_temp)

    ! postprocess
    call self%backend%poisson_fft%fft_postprocess

    ! call backward FFT
    call self%backend%poisson_fft%fft_backward(p_temp)

    ! reorder back to our specialist data structure from 3D Cartesian
    call self%backend%reorder(pressure, p_temp, RDR_C2Z)

    call self%backend%allocator%release_block(p_temp)

  end subroutine poisson_fft

  subroutine poisson_cg(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

  end subroutine poisson_cg

  subroutine output(self, t)
    implicit none

    class(solver_t), intent(in) :: self
    real(dp), intent(in) :: t

    class(field_t), pointer :: du, dv, dw, div_u
    class(field_t), pointer :: u_out
    real(dp) :: enstrophy, div_u_max, div_u_mean
    integer :: ierr

    if (self%mesh%par%is_root()) print *, 'time = ', t

    du => self%backend%allocator%get_block(DIR_X, VERT)
    dv => self%backend%allocator%get_block(DIR_X, VERT)
    dw => self%backend%allocator%get_block(DIR_X, VERT)

    call self%curl(du, dv, dw, self%u, self%v, self%w)
    enstrophy = 0.5_dp*(self%backend%scalar_product(du, du) &
                        + self%backend%scalar_product(dv, dv) &
                        + self%backend%scalar_product(dw, dw))/self%ngrid
    if (self%mesh%par%is_root()) print *, 'enstrophy:', enstrophy

    call self%backend%allocator%release_block(du)
    call self%backend%allocator%release_block(dv)
    call self%backend%allocator%release_block(dw)

    div_u => self%backend%allocator%get_block(DIR_Z)

    call self%divergence_v2p(div_u, self%u, self%v, self%w)

    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, div_u)

    call self%backend%allocator%release_block(div_u)

    div_u_max = maxval(abs(u_out%data))
    div_u_mean = sum(abs(u_out%data))/self%ngrid

    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'div u max mean:', div_u_max, div_u_mean

  end subroutine output

  subroutine run(self)
    implicit none

    class(solver_t), intent(inout) :: self

    class(field_t), pointer :: du, dv, dw, div_u, pressure, dpdx, dpdy, dpdz
    class(field_t), pointer :: u_out, v_out, w_out

    real(dp) :: t
    integer :: i, j

    if (self%mesh%par%is_root()) print *, 'initial conditions'
    t = 0._dp
    call self%output(t)

    if (self%mesh%par%is_root()) print *, 'start run'

    do i = 1, self%n_iters
      do j = 1, self%time_integrator%nstage
        du => self%backend%allocator%get_block(DIR_X)
        dv => self%backend%allocator%get_block(DIR_X)
        dw => self%backend%allocator%get_block(DIR_X)

        call self%transeq(du, dv, dw, self%u, self%v, self%w)

        ! time integration
        call self%time_integrator%step(self%u, self%v, self%w, &
                                       du, dv, dw, self%dt)

        call self%backend%allocator%release_block(du)
        call self%backend%allocator%release_block(dv)
        call self%backend%allocator%release_block(dw)

        ! pressure
        div_u => self%backend%allocator%get_block(DIR_Z)

        call self%divergence_v2p(div_u, self%u, self%v, self%w)

        pressure => self%backend%allocator%get_block(DIR_Z, CELL)

        call self%poisson(pressure, div_u)

        call self%backend%allocator%release_block(div_u)

        dpdx => self%backend%allocator%get_block(DIR_X)
        dpdy => self%backend%allocator%get_block(DIR_X)
        dpdz => self%backend%allocator%get_block(DIR_X)

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

      if (mod(i, self%n_output) == 0) then
        t = i*self%dt
        call self%output(t)
      end if
    end do

    if (self%mesh%par%is_root()) print *, 'run end'

    ! Below is for demonstrating purpuses only, to be removed when we have
    ! proper I/O in place.
    u_out => self%host_allocator%get_block(DIR_C)
    v_out => self%host_allocator%get_block(DIR_C)
    w_out => self%host_allocator%get_block(DIR_C)

    call self%backend%get_field_data(u_out%data, self%u)
    call self%backend%get_field_data(v_out%data, self%v)
    call self%backend%get_field_data(w_out%data, self%w)

    if (self%mesh%par%is_root()) then
      print *, 'norms', norm2(u_out%data), norm2(v_out%data), norm2(w_out%data)
    end if

    call self%host_allocator%release_block(u_out)
    call self%host_allocator%release_block(v_out)
    call self%host_allocator%release_block(w_out)

  end subroutine run

end module m_solver
