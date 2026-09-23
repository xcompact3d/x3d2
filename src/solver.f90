module m_solver
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, get_argument, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      RDR_Z2C, RDR_C2Z, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, CELL, &
                      X_FACE, Y_FACE, BC_NEUMANN, BC_DIRICHLET
  use m_config, only: solver_config_t, les_config_t
  use m_field, only: field_t, flist_t
  use m_ibm, only: ibm_t
  use m_les, only: les_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t, tdsops_t
  use m_time_integrator, only: time_intg_t
  use m_vector_calculus, only: vector_calculus_t

  implicit none

  type :: solver_t
      !! solver class defines the Incompact3D algorithm at a very high level.
      !!
      !! Procedures defined here that are part of the Incompact3D algorithm
      !! are: transeq, divergence, poisson, and gradient.
      !!
      !! The operations these high level procedures require are provided by
      !! the relevant backend implementations.
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
    real(dp), dimension(:), allocatable :: nu_species
    integer :: n_iters, n_output
    integer :: current_iter = 0
    integer :: ngrid
    integer :: nvars = 3
    integer :: nspecies = 0

    class(field_t), pointer :: u, v, w
    class(field_t), pointer :: pressure => null()      !! Pressure on CELL grid (DIR_Z)
    class(field_t), pointer :: pressure_vert => null() !! Pressure on VERT grid (DIR_X)
    logical :: keep_pressure = .false.                 !! If true, persist pressure for output
    !> Last projection's velocity correction, kept when the mesh has a
    !> Dirichlet face so precorrect_walls can reuse it
    class(field_t), pointer :: dpdx_last => null()
    class(field_t), pointer :: dpdy_last => null()
    class(field_t), pointer :: dpdz_last => null()
    class(field_t), pointer :: vort => null()  !! Vorticity magnitude on VERT grid
    class(field_t), pointer :: qcrit => null() !! Q-criterion on VERT grid
    type(flist_t), dimension(:), pointer :: species => null()

    class(base_backend_t), pointer :: backend
    type(mesh_t), pointer :: mesh
    type(time_intg_t) :: time_integrator
    type(allocator_t), pointer :: host_allocator
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    type(vector_calculus_t) :: vector_calculus
    type(ibm_t) :: ibm
    type(les_t) :: les
    logical :: ibm_on
    logical :: spatial_filter = .false. !! explicit low-pass filter on velocity
    procedure(poisson_solver), pointer :: poisson => null()
    procedure(transport_equation), pointer :: transeq => null()
  contains
    procedure :: transeq_species
    procedure :: apply_les
    procedure :: apply_spatial_filter
    procedure :: finalise
    procedure :: pressure_correction
    procedure :: precorrect_walls
    procedure :: divergence_v2p
    procedure :: gradient_p2v
    procedure :: curl
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

    subroutine transport_equation(self, rhs, variables)
      import :: solver_t
      import :: flist_t
      implicit none

      class(solver_t) :: self
      type(flist_t), intent(inout) :: rhs(:), variables(:)
    end subroutine transport_equation
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

    type(solver_config_t) :: solver_cfg
    type(les_config_t) :: les_cfg
    integer :: i

    solver%backend => backend
    solver%mesh => mesh
    solver%host_allocator => host_allocator

    allocate (solver%xdirps, solver%ydirps, solver%zdirps)
    solver%xdirps%dir = DIR_X
    solver%ydirps%dir = DIR_Y
    solver%zdirps%dir = DIR_Z

    solver%vector_calculus = vector_calculus_t(solver%backend)

    solver%u => solver%backend%allocator%get_block(DIR_X)
    solver%v => solver%backend%allocator%get_block(DIR_X)
    solver%w => solver%backend%allocator%get_block(DIR_X)

    call solver_cfg%read(nml_file=get_argument(1))
    call les_cfg%read(nml_file=get_argument(1))
    solver%les = les_t(les_cfg)
    if (solver%mesh%par%is_root()) &
      print *, 'LES model: ', trim(solver%les%model)

    ! Add transported species
    solver%nspecies = solver_cfg%n_species
    if (solver%nspecies > 0) then
      ! Increase the number of variables
      solver%nvars = solver%nvars + solver%nspecies

      ! Init the diffusivity coefficients for species
      allocate (solver%nu_species(solver%nspecies))
      solver%nu_species = 1._dp/solver_cfg%Re/solver_cfg%pr_species

      ! Get blocks for the species
      allocate (solver%species(solver%nspecies))
      do i = 1, solver%nspecies
        solver%species(i)%ptr => solver%backend%allocator%get_block(DIR_X)
      end do
    end if

    solver%time_integrator = time_intg_t(solver%backend, &
                                         solver%backend%allocator, &
                                         solver_cfg%time_intg, solver%nvars)
    if (solver%mesh%par%is_root()) then
      print *, solver_cfg%time_intg//' time integrator instantiated'
    end if

    solver%dt = solver_cfg%dt
    solver%nu = 1._dp/solver_cfg%Re
    solver%n_iters = solver_cfg%n_iters
    solver%n_output = solver_cfg%n_output
    solver%ngrid = product(solver%mesh%get_global_dims(VERT))

    ! Allocate and set the tdsops. The filter operators are only built when
    ! the case asks for them, so cases that do not filter carry no extra state.
    solver%spatial_filter = solver_cfg%spatial_filter
    if (solver%spatial_filter) then
      if (solver%mesh%par%is_root()) &
        print *, 'Spatial filter on, alpha =', solver_cfg%filter_alpha
      call allocate_tdsops( &
        solver%xdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme, filter_alpha=solver_cfg%filter_alpha &
        )
      call allocate_tdsops( &
        solver%ydirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme, filter_alpha=solver_cfg%filter_alpha &
        )
      call allocate_tdsops( &
        solver%zdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme, filter_alpha=solver_cfg%filter_alpha &
        )
    else
      call allocate_tdsops( &
        solver%xdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme &
        )
      call allocate_tdsops( &
        solver%ydirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme &
        )
      call allocate_tdsops( &
        solver%zdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
        solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
        solver_cfg%stagder_scheme &
        )
    end if

    select case (trim(solver_cfg%poisson_solver_type))
    case ('FFT')
      if (solver%mesh%par%is_root()) print *, 'Poisson solver: FFT'
      call solver%backend%init_poisson_fft(solver%mesh, solver%xdirps, &
                                           solver%ydirps, solver%zdirps, &
                                           solver_cfg%lowmem_fft)
      solver%poisson => poisson_fft
    case ('CG')
      if (solver%mesh%par%is_root()) &
        print *, 'Poisson solver: CG, not yet implemented'
      solver%poisson => poisson_cg
    case default
      error stop 'poisson_solver_type is not valid. Use "FFT" or "CG".'
    end select

    ! Initialize the IBM module
    solver%ibm_on = solver_cfg%ibm_on
    if (solver%ibm_on) &
      solver%ibm = ibm_t(backend, mesh, host_allocator)

    if (solver_cfg%lowmem_transeq) then
      solver%transeq => transeq_lowmem
    else
      solver%transeq => transeq_default
    end if

  end function init

  subroutine allocate_tdsops(dirps, backend, mesh, der1st_scheme, &
                             der2nd_scheme, interpl_scheme, stagder_scheme, &
                             filter_alpha)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    character(*), intent(in) :: der1st_scheme, der2nd_scheme, &
                                interpl_scheme, stagder_scheme
    !! When present, also build the low-pass filter operators for this
    !! direction. Absent means the case does not filter.
    real(dp), optional, intent(in) :: filter_alpha

    integer :: dir, bc_start, bc_end, bc_mp_start, bc_mp_end, n_vert, n_cell, i
    real(dp) :: d

    dir = dirps%dir
    d = mesh%geo%d(dir)

    bc_start = mesh%grid%BCs(dir, 1)
    bc_end = mesh%grid%BCs(dir, 2)

    ! For the FFT based poisson solver, the BC for the pressure has to be
    ! Neumann. This is not strictly compatible with Navier-Stokes equations,
    ! but it does not affect the quality of the simulation.
    ! Thus, if the BC is Dirichlet, we enforce Neumann for the midpoint &
    ! operators. (It could be BC_HALO too, if so we just keep it as is)
    if (bc_start == BC_DIRICHLET) then
      bc_mp_start = BC_NEUMANN
    else
      bc_mp_start = bc_start
    end if
    if (bc_end == BC_DIRICHLET) then
      bc_mp_end = BC_NEUMANN
    else
      bc_mp_end = bc_end
    end if

    n_vert = mesh%get_n(dir, VERT)
    n_cell = mesh%get_n(dir, CELL)

    call backend%alloc_tdsops( &
      dirps%der1st, n_vert, d, 'first-deriv', der1st_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%der1st_sym, n_vert, d, 'first-deriv', der1st_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds(1:n_vert, dir), &
      sym=.true. &
      )
    call backend%alloc_tdsops( &
      dirps%der2nd, n_vert, d, 'second-deriv', der2nd_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds2(1:n_vert, dir), &
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%der2nd_sym, n_vert, d, 'second-deriv', der2nd_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds2(1:n_vert, dir), &
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir), &
      sym=.true. &
      )
    call backend%alloc_tdsops( &
      dirps%stagder_v2p, n_cell, d, 'stag-deriv', stagder_scheme, &
      bc_mp_start, bc_mp_end, from_to='v2p', &
      stretch=mesh%geo%midp_ds(1:n_cell, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%stagder_p2v, n_vert, d, 'stag-deriv', stagder_scheme, &
      bc_mp_start, bc_mp_end, from_to='p2v', &
      stretch=mesh%geo%vert_ds(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%interpl_v2p, n_cell, d, 'interpolate', interpl_scheme, &
      bc_mp_start, bc_mp_end, from_to='v2p', stretch=[(1._dp, i=1, n_cell)] &
      )
    call backend%alloc_tdsops( &
      dirps%interpl_p2v, n_vert, d, 'interpolate', interpl_scheme, &
      bc_mp_start, bc_mp_end, from_to='p2v', stretch=[(1._dp, i=1, n_vert)] &
      )

    if (present(filter_alpha)) then
      ! Same parity split as the first derivatives: the component along this
      ! direction is odd across a free-slip boundary, the other two even.
      call backend%alloc_tdsops( &
        dirps%lowpass, n_vert, d, 'filter', der1st_scheme, &
        bc_start, bc_end, filter_alpha=filter_alpha &
        )
      call backend%alloc_tdsops( &
        dirps%lowpass_sym, n_vert, d, 'filter', der1st_scheme, &
        bc_start, bc_end, sym=.true., filter_alpha=filter_alpha &
        )
      ! The filter's system is only marginally diagonally dominant as
      ! filter_alpha approaches 0.5, so the distributed solver's truncation
      ! is no longer negligible: at Incompact3d's 0.49 it leaves a 1e-3 error
      ! at zero wavenumber on 64 points. Solve it exactly with Thomas where
      ! the direction is not decomposed; Thomas is serial along the line.
      if (mesh%par%nproc_dir(dir) == 1) then
        dirps%lowpass%prefer_thomas = .true.
        dirps%lowpass_sym%prefer_thomas = .true.
      end if
    end if

  end subroutine

  subroutine transeq_default(self, rhs, variables)
    !! Skew-symmetric form of convection-diffusion terms in the
    !! incompressible Navier-Stokes momentum equations, excluding
    !! pressure terms.
    !! Inputs from velocity grid and outputs to velocity grid.
    implicit none

    class(solver_t) :: self
    type(flist_t), intent(inout) :: rhs(:)
    type(flist_t), intent(inout) :: variables(:)

    class(field_t), pointer :: u_y, v_y, w_y, u_z, v_z, w_z, &
      du_y, dv_y, dw_y, du_z, dv_z, dw_z, &
      du, dv, dw, u, v, w

    du => rhs(1)%ptr
    dv => rhs(2)%ptr
    dw => rhs(3)%ptr
    u => variables(1)%ptr
    v => variables(2)%ptr
    w => variables(3)%ptr

    ! -1/2(nabla u curl u + u nabla u) + nu nablasq u

    ! call derivatives in x direction. Based on the run time arguments this
    ! executes a distributed algorithm or the Thomas algorithm.
    call self%backend%transeq_x(du, dv, dw, u, v, w, self%nu, self%xdirps)

    ! request fields from the allocator
    u_y => self%backend%allocator%get_block(DIR_Y)
    v_y => self%backend%allocator%get_block(DIR_Y)
    w_y => self%backend%allocator%get_block(DIR_Y)
    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, u, RDR_X2Y)
    call self%backend%reorder(v_y, v, RDR_X2Y)
    call self%backend%reorder(w_y, w, RDR_X2Y)

    ! similar to the x direction, obtain derivatives in y.
    call self%backend%transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, &
                                self%nu, self%ydirps)

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
    u_z => self%backend%allocator%get_block(DIR_Z)
    v_z => self%backend%allocator%get_block(DIR_Z)
    w_z => self%backend%allocator%get_block(DIR_Z)
    du_z => self%backend%allocator%get_block(DIR_Z)
    dv_z => self%backend%allocator%get_block(DIR_Z)
    dw_z => self%backend%allocator%get_block(DIR_Z)

    ! reorder from x to z
    call self%backend%reorder(u_z, u, RDR_X2Z)
    call self%backend%reorder(v_z, v, RDR_X2Z)
    call self%backend%reorder(w_z, w, RDR_X2Z)

    ! get the derivatives in z
    call self%backend%transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, &
                                self%nu, self%zdirps)

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

    call self%apply_les(rhs, variables)

    ! Convection-diffusion for species
    if (self%nspecies > 0) then
      call self%transeq_species(rhs(4:), variables)
    end if

  end subroutine transeq_default

  subroutine transeq_lowmem(self, rhs, variables)
    !! low memory version of the transport equation, roughly %2 slower overall
    implicit none

    class(solver_t) :: self
    type(flist_t), intent(inout) :: rhs(:)
    type(flist_t), intent(inout) :: variables(:)

    class(field_t), pointer :: u_y, v_y, w_y, u_z, v_z, w_z, &
      du_y, dv_y, dw_y, du_z, dv_z, dw_z, du, dv, dw, u, v, w

    du => rhs(1)%ptr
    dv => rhs(2)%ptr
    dw => rhs(3)%ptr
    u => variables(1)%ptr
    v => variables(2)%ptr
    w => variables(3)%ptr

    ! -1/2(nabla u curl u + u nabla u) + nu nablasq u

    ! call derivatives in x direction. Based on the run time arguments this
    ! executes a distributed algorithm or the Thomas algorithm.
    call self%backend%transeq_x(du, dv, dw, u, v, w, self%nu, self%xdirps)

    ! request fields from the allocator
    u_y => self%backend%allocator%get_block(DIR_Y)
    v_y => self%backend%allocator%get_block(DIR_Y)
    w_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, u, RDR_X2Y)
    call self%backend%reorder(v_y, v, RDR_X2Y)
    call self%backend%reorder(w_y, w, RDR_X2Y)

    ! now release the x-directional fields for saving memory
    call self%backend%allocator%release_block(u)
    call self%backend%allocator%release_block(v)
    call self%backend%allocator%release_block(w)

    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

    ! similar to the x direction, obtain derivatives in y.
    call self%backend%transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, &
                                self%nu, self%ydirps)

    call self%backend%sum_yintox(du, du_y)
    call self%backend%sum_yintox(dv, dv_y)
    call self%backend%sum_yintox(dw, dw_y)

    call self%backend%allocator%release_block(du_y)
    call self%backend%allocator%release_block(dv_y)
    call self%backend%allocator%release_block(dw_y)

    ! just like in y direction, get some fields for the z derivatives.
    u_z => self%backend%allocator%get_block(DIR_Z)
    v_z => self%backend%allocator%get_block(DIR_Z)
    w_z => self%backend%allocator%get_block(DIR_Z)

    ! reorder from y to z
    call self%backend%reorder(u_z, u_y, RDR_Y2Z)
    call self%backend%reorder(v_z, v_y, RDR_Y2Z)
    call self%backend%reorder(w_z, w_y, RDR_Y2Z)

    ! we don't need the velocities in y orientation any more, so release
    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(v_y)
    call self%backend%allocator%release_block(w_y)

    du_z => self%backend%allocator%get_block(DIR_Z)
    dv_z => self%backend%allocator%get_block(DIR_Z)
    dw_z => self%backend%allocator%get_block(DIR_Z)

    ! get the derivatives in z
    call self%backend%transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, &
                                self%nu, self%zdirps)

    ! gather all the contributions into the x result array
    call self%backend%sum_zintox(du, du_z)
    call self%backend%sum_zintox(dv, dv_z)
    call self%backend%sum_zintox(dw, dw_z)

    ! release all the unnecessary blocks.
    call self%backend%allocator%release_block(du_z)
    call self%backend%allocator%release_block(dv_z)
    call self%backend%allocator%release_block(dw_z)

    u => self%backend%allocator%get_block(DIR_X)
    v => self%backend%allocator%get_block(DIR_X)
    w => self%backend%allocator%get_block(DIR_X)

    ! reorder from z to x
    call self%backend%reorder(u, u_z, RDR_Z2X)
    call self%backend%reorder(v, v_z, RDR_Z2X)
    call self%backend%reorder(w, w_z, RDR_Z2X)

    ! there is no need to keep velocities in z orientation around, so release
    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(v_z)
    call self%backend%allocator%release_block(w_z)

    variables(1)%ptr => u
    variables(2)%ptr => v
    variables(3)%ptr => w
    self%u => u
    self%v => v
    self%w => w

    call self%apply_les(rhs, variables)

    ! Convection-diffusion for species
    if (self%nspecies > 0) then
      call self%transeq_species(rhs(4:), variables)
    end if

  end subroutine transeq_lowmem

  subroutine apply_spatial_filter(self)
    !! Apply the explicit low-pass filter to the velocity, as Incompact3d does
    !! for wall-modelled ABL runs (ifilter=1, C_filter).
    !!
    !! Compact schemes cannot dissipate the 2*dx mode and the staggered
    !! pressure projection cannot see it, yet the collocated derivatives in
    !! the skew-symmetric convection do. Left alone that mode drives the
    !! collocated divergence away from zero, and the u*div(u) half of the
    !! convection then acts as a spurious momentum source.
    class(solver_t), intent(inout) :: self

    if (.not. self%spatial_filter) return

    call filter_field(self, self%u, self%xdirps, self%ydirps, self%zdirps, 1)
    call filter_field(self, self%v, self%xdirps, self%ydirps, self%zdirps, 2)
    call filter_field(self, self%w, self%xdirps, self%ydirps, self%zdirps, 3)

  end subroutine apply_spatial_filter

  subroutine filter_field(self, f, xdirps, ydirps, zdirps, component)
    !! Filter one velocity component in all three directions, in place.
    !!
    !! Across a free-slip boundary the component along that direction is odd
    !! and the other two are even, the same parity split the first derivatives
    !! use, so `component` selects which operator each direction applies.
    class(solver_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, intent(in) :: component

    call filter_in_dir(self, f, xdirps, component == 1)
    call filter_in_dir(self, f, ydirps, component == 2)
    call filter_in_dir(self, f, zdirps, component == 3)

  end subroutine filter_field

  subroutine filter_in_dir(self, f, dirps, is_normal)
    !! Filter a DIR_X field along one direction, in place.
    class(solver_t), intent(inout) :: self
    class(field_t), intent(inout) :: f
    type(dirps_t), target, intent(in) :: dirps
    !! .true. when f is the component along dirps%dir, which is the odd one
    !! across a free-slip boundary.
    logical, intent(in) :: is_normal

    class(field_t), pointer :: f_dir, filtered_dir, filtered
    class(tdsops_t), pointer :: op

    if (is_normal) then
      op => dirps%lowpass
    else
      op => dirps%lowpass_sym
    end if

    select case (dirps%dir)
    case (DIR_X)
      filtered => self%backend%allocator%get_block(DIR_X, f%data_loc)
      call self%backend%tds_solve(filtered, f, op)
      call self%backend%veccopy(f, filtered)
      call self%backend%allocator%release_block(filtered)
    case (DIR_Y, DIR_Z)
      f_dir => self%backend%allocator%get_block(dirps%dir)
      filtered_dir => self%backend%allocator%get_block(dirps%dir)
      filtered => self%backend%allocator%get_block(DIR_X, f%data_loc)
      if (dirps%dir == DIR_Y) then
        call self%backend%reorder(f_dir, f, RDR_X2Y)
      else
        call self%backend%reorder(f_dir, f, RDR_X2Z)
      end if
      call self%backend%tds_solve(filtered_dir, f_dir, op)
      if (dirps%dir == DIR_Y) then
        call self%backend%reorder(filtered, filtered_dir, RDR_Y2X)
      else
        call self%backend%reorder(filtered, filtered_dir, RDR_Z2X)
      end if
      call self%backend%veccopy(f, filtered)
      call self%backend%allocator%release_block(f_dir)
      call self%backend%allocator%release_block(filtered_dir)
      call self%backend%allocator%release_block(filtered)
    case default
      error stop 'Invalid direction in spatial filter.'
    end select

  end subroutine filter_in_dir

  subroutine apply_les(self, rhs, variables)
    !! Add the configured explicit SGS closure to the momentum RHS.
    class(solver_t), intent(inout) :: self
    type(flist_t), intent(inout) :: rhs(:)
    type(flist_t), intent(in) :: variables(:)

    call self%les%apply_sgs_stress( &
      self%backend, self%mesh, &
      rhs(1)%ptr, rhs(2)%ptr, rhs(3)%ptr, &
      variables(1)%ptr, variables(2)%ptr, variables(3)%ptr, &
      self%xdirps, self%ydirps, self%zdirps)
  end subroutine apply_les

  subroutine finalise(self)
    !! Release resources owned by the solver and its runtime models.
    class(solver_t), intent(inout) :: self

    call self%les%finalise(self%backend)
    if (associated(self%dpdx_last)) then
      call self%backend%allocator%release_block(self%dpdx_last)
      call self%backend%allocator%release_block(self%dpdy_last)
      call self%backend%allocator%release_block(self%dpdz_last)
    end if
  end subroutine finalise

  subroutine transeq_species(self, rhs, variables)
    !! Skew-symmetric form of convection-diffusion terms in the
    !! species equation.
    !! Inputs from velocity grid and outputs to velocity grid.
    implicit none

    class(solver_t) :: self
    type(flist_t), intent(inout) :: rhs(:)
    type(flist_t), intent(in) :: variables(:)

    integer :: i
    class(field_t), pointer :: u, v, w, &
      v_y, spec_y, dspec_y, &
      w_z, spec_z, dspec_z

    ! Map the velocity vector
    u => variables(1)%ptr
    v => variables(2)%ptr
    w => variables(3)%ptr

    ! FIXME later
    ! Minor optimization
    ! species could start with z convection-diffusion
    ! velocity components are ready to use in the z dir.

    ! derivatives in x
    do i = 1, size(rhs)
      call self%backend%transeq_species(rhs(i)%ptr, u, &
                                        variables(3 + i)%ptr, &
                                        self%nu_species(i), &
                                        self%xdirps, &
                                        i <= 1)
    end do

    ! Request blocks
    v_y => self%backend%allocator%get_block(DIR_Y)
    spec_y => self%backend%allocator%get_block(DIR_Y)
    dspec_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder velocity
    call self%backend%reorder(v_y, v, RDR_X2Y)

    do i = 1, size(rhs)

      ! reorder spec in y
      call self%backend%reorder(spec_y, variables(3 + i)%ptr, RDR_X2Y)

      ! y-derivatives
      call self%backend%transeq_species(dspec_y, v_y, &
                                        spec_y, &
                                        self%nu_species(i), &
                                        self%ydirps, &
                                        i <= 1)

      ! sum_yintox
      call self%backend%sum_yintox(rhs(i)%ptr, dspec_y)

    end do

    ! Release blocks
    call self%backend%allocator%release_block(v_y)
    call self%backend%allocator%release_block(spec_y)
    call self%backend%allocator%release_block(dspec_y)

    ! Request blocks
    w_z => self%backend%allocator%get_block(DIR_Z)
    spec_z => self%backend%allocator%get_block(DIR_Z)
    dspec_z => self%backend%allocator%get_block(DIR_Z)

    ! reorder velocity
    call self%backend%reorder(w_z, w, RDR_X2Z)

    do i = 1, size(rhs)

      ! reorder spec in z
      call self%backend%reorder(spec_z, variables(3 + i)%ptr, RDR_X2Z)

      ! z-derivatives
      call self%backend%transeq_species(dspec_z, w_z, &
                                        spec_z, &
                                        self%nu_species(i), &
                                        self%zdirps, &
                                        i <= 1)

      ! sum_zintox
      call self%backend%sum_zintox(rhs(i)%ptr, dspec_z)

    end do

    ! Release blocks
    call self%backend%allocator%release_block(w_z)
    call self%backend%allocator%release_block(spec_z)
    call self%backend%allocator%release_block(dspec_z)

  end subroutine transeq_species

  subroutine divergence_v2p(self, div_u, u, v, w)
    !! Wrapper for divergence_v2p
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: div_u
    class(field_t), intent(in) :: u, v, w

    call self%vector_calculus%divergence_v2c( &
      div_u, u, v, w, &
      self%xdirps%stagder_v2p, self%xdirps%interpl_v2p, &
      self%ydirps%stagder_v2p, self%ydirps%interpl_v2p, &
      self%zdirps%stagder_v2p, self%zdirps%interpl_v2p &
      )

  end subroutine divergence_v2p

  subroutine gradient_p2v(self, dpdx, dpdy, dpdz, pressure)
    !! Wrapper for gradient_p2v
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: dpdx, dpdy, dpdz
    class(field_t), intent(in) :: pressure

    call self%vector_calculus%gradient_c2v( &
      dpdx, dpdy, dpdz, pressure, &
      self%xdirps%stagder_p2v, self%xdirps%interpl_p2v, &
      self%ydirps%stagder_p2v, self%ydirps%interpl_p2v, &
      self%zdirps%stagder_p2v, self%zdirps%interpl_p2v &
      )

  end subroutine gradient_p2v

  subroutine curl(self, o_i_hat, o_j_hat, o_k_hat, u, v, w)
    !! Wrapper for curl
    implicit none

    class(solver_t) :: self
    !> Vector components of the output vector field Omega
    class(field_t), intent(inout) :: o_i_hat, o_j_hat, o_k_hat
    class(field_t), intent(in) :: u, v, w

    call self%vector_calculus%curl( &
      o_i_hat, o_j_hat, o_k_hat, u, v, w, &
      self%xdirps%der1st, self%ydirps%der1st, self%zdirps%der1st &
      )

  end subroutine curl

  subroutine poisson_fft(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

    class(field_t), pointer :: p_temp, temp

    ! reorder into 3D Cartesian data structure
    p_temp => self%backend%allocator%get_block(DIR_C)
    call self%backend%reorder(p_temp, div_u, RDR_Z2C)

    temp => self%backend%allocator%get_block(DIR_C)

    ! solve poisson equation with FFT based approach
    call self%backend%poisson_fft%solve_poisson(p_temp, temp)

    call self%backend%allocator%release_block(temp)

    ! reorder back to our specialist data structure from 3D Cartesian
    call self%backend%reorder(pressure, p_temp, RDR_C2Z)

    call self%backend%allocator%release_block(p_temp)

  end subroutine poisson_fft

  subroutine poisson_cg(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

    ! set the pressure field to 0 so that we can do performance tests easily
    ! this will be removed once the CG solver is implemented of course
    call pressure%fill(0._dp)

  end subroutine poisson_cg

  subroutine pressure_correction(self, u, v, w)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: u, v, w

    class(field_t), pointer :: div_u, p, dpdx, dpdy, dpdz

    div_u => self%backend%allocator%get_block(DIR_Z)

    call self%divergence_v2p(div_u, u, v, w)

    if (self%keep_pressure) then
      ! Persist pressure for snapshot output
      if (.not. associated(self%pressure)) then
        self%pressure => self%backend%allocator%get_block(DIR_Z, CELL)
      end if
      p => self%pressure
    else
      ! Temporary pressure, released after use
      p => self%backend%allocator%get_block(DIR_Z)
    end if

    call self%poisson(p, div_u)

    call self%backend%allocator%release_block(div_u)

    dpdx => self%backend%allocator%get_block(DIR_X)
    dpdy => self%backend%allocator%get_block(DIR_X)
    dpdz => self%backend%allocator%get_block(DIR_X)

    call self%gradient_p2v(dpdx, dpdy, dpdz, p)

    if (.not. self%keep_pressure) then
      call self%backend%allocator%release_block(p)
    end if

    ! velocity correction
    call self%backend%vecadd(-1._dp, dpdx, 1._dp, u)
    call self%backend%vecadd(-1._dp, dpdy, 1._dp, v)
    call self%backend%vecadd(-1._dp, dpdz, 1._dp, w)

    if (any(self%mesh%grid%BCs_global == BC_DIRICHLET)) then
      ! Keep this correction for the next precorrect_walls
      if (associated(self%dpdx_last)) then
        call self%backend%allocator%release_block(self%dpdx_last)
        call self%backend%allocator%release_block(self%dpdy_last)
        call self%backend%allocator%release_block(self%dpdz_last)
      end if
      self%dpdx_last => dpdx
      self%dpdy_last => dpdy
      self%dpdz_last => dpdz
    else
      call self%backend%allocator%release_block(dpdx)
      call self%backend%allocator%release_block(dpdy)
      call self%backend%allocator%release_block(dpdz)
    end if

  end subroutine pressure_correction

  subroutine precorrect_walls(self, u, v, w)
    !! Incompact3d's pre_correc: on each Dirichlet face, add the previous
    !! projection's tangential correction to the value set by the case. The
    !! projection that follows then leaves the face at its prescribed value
    !! plus dt*(grad p_old - grad p_new), instead of -dt*grad p_new.
    !! Exact for AB, where every step uses the same dt; RK stages would need
    !! Incompact3d's gdt ratio. Before the first projection (or after a
    !! restart, as in Incompact3d) there is no correction to add.
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: u, v, w

    integer :: bcs(3, 2)

    if (.not. associated(self%dpdx_last)) return

    ! Local BCs, so only the ranks that own a face touch it
    bcs = self%mesh%grid%BCs
    if (any(bcs(1, :) == BC_DIRICHLET)) then
      call self%backend%field_add_face_from_field( &
        v, self%dpdy_last, X_FACE, bcs(1, 1), bcs(1, 2))
      call self%backend%field_add_face_from_field( &
        w, self%dpdz_last, X_FACE, bcs(1, 1), bcs(1, 2))
    end if
    if (any(bcs(2, :) == BC_DIRICHLET)) then
      call self%backend%field_add_face_from_field( &
        u, self%dpdx_last, Y_FACE, bcs(2, 1), bcs(2, 2))
      call self%backend%field_add_face_from_field( &
        w, self%dpdz_last, Y_FACE, bcs(2, 1), bcs(2, 2))
    end if
    if (any(bcs(3, :) == BC_DIRICHLET)) &
      error stop 'precorrect_walls: Dirichlet z-faces are not supported.'

  end subroutine precorrect_walls

end module m_solver
