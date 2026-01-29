module m_omp_backend
  !! OpenMP/CPU backend implementation for X3D2 solver operations.
  !!
  !! This module provides the CPU-based backend using OpenMP for shared-memory
  !! parallelism and MPI for distributed-memory parallelism. It implements all
  !! abstract backend operations defined in `base_backend_t`.
  !!
  !! **Parallelisation Strategy:**
  !!
  !! - **MPI**: Domain decomposition across nodes/processes
  !! - **OpenMP**: Thread parallelism within each MPI rank
  !! - **Hybrid MPI+OpenMP**: Enables efficient use of multi-core clusters
  !!
  !! **Key Features:**
  !!
  !! - Compact finite difference operators (tridiagonal solves)
  !! - Halo exchange for distributed derivatives
  !! - FFT-based Poisson solver integration
  !! - Vectorised array operations
  !! - Optimised data reordering between decomposition directions
  !!
  !! **Memory Management:**
  !!
  !! - Send/receive buffers for MPI halo exchange (`u`, `v`, `w`, `du`, `dud`, `d2u`)
  !! - Buffers sized based on largest decomposition direction
  !! - Persistent buffers to avoid repeated allocation
  !!
  !! **Solver Operations:**
  !!
  !! - `transeq`: Transport equation terms with halo exchange
  !! - `tds_solve`: Tridiagonal system solves (Thomas algorithm)
  !! - `reorder`: Data layout transformations (`DIR_X`, `DIR_Y`, `DIR_Z`)
  !! - Field operations: copy, add, multiply, integrate, etc.
  !!
  !! **Note:** This backend requires 2DECOMP&FFT library for FFT operations
  !! when using the spectral Poisson solver.
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, MPI_X3D2_DP, get_dirs_from_rdr, move_data_loc, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, NULL_LOC, &
                      X_FACE, Y_FACE, Z_FACE, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_ordering, only: get_index_reordering
  use m_tdsops, only: dirps_t, tdsops_t

  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_tds_compact, exec_dist_transeq_compact
  use m_omp_sendrecv, only: sendrecv_fields

  implicit none

  private :: transeq_halo_exchange, transeq_dist_component

  type, extends(base_backend_t) :: omp_backend_t
    !! OpenMP/CPU backend for solver operations.
    !!
    !! Extends `base_backend_t` with CPU-specific implementations using
    !! OpenMP for threading and MPI for distributed parallelism. Maintains
    !! communication buffers for halo exchange operations.
    !!
    !! **Communication Buffers:**
    !! Arrays sized (SZ, n_halo, n_groups) where:
    !! - SZ: maximum pencil size for data reordering
    !! - n_halo: halo region depth (typically 4 for compact schemes)
    !! - n_groups: maximum number of groups across all directions
    !!
    !! Buffer naming convention: {field}_{send/recv}_{s/e}
    !! - field: u, v, w (velocity), du, dud, d2u (derivatives)
    !! - send/recv: data direction
    !! - s/e: start/end of domain (neighbouring ranks)
    !character(len=*), parameter :: name = 'omp' !! Backend identifier
    real(dp), allocatable, dimension(:, :, :) :: &
      u_recv_s, u_recv_e, u_send_s, u_send_e, &     !! Velocity u halo buffers
      v_recv_s, v_recv_e, v_send_s, v_send_e, &     !! Velocity v halo buffers
      w_recv_s, w_recv_e, w_send_s, w_send_e, &     !! Velocity w halo buffers
      du_send_s, du_send_e, du_recv_s, du_recv_e, & !! First derivative buffers
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, & !! Mixed derivative buffers
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e    !! Second derivative buffers
  contains
    procedure :: alloc_tdsops => alloc_omp_tdsops          !! Allocate tridiagonal operators
    procedure :: transeq_x => transeq_x_omp                !! Transport equation in X
    procedure :: transeq_y => transeq_y_omp                !! Transport equation in Y
    procedure :: transeq_z => transeq_z_omp                !! Transport equation in Z
    procedure :: transeq_species => transeq_species_omp    !! Transport for species/scalars
    procedure :: tds_solve => tds_solve_omp                !! Tridiagonal solve
    procedure :: reorder => reorder_omp                    !! Data reordering
    procedure :: sum_yintox => sum_yintox_omp              !! Sum Y data into X
    procedure :: sum_zintox => sum_zintox_omp              !! Sum Z data into X
    procedure :: veccopy => veccopy_omp                    !! Vector copy
    procedure :: vecadd => vecadd_omp                      !! Vector add
    procedure :: vecmult => vecmult_omp                    !! Vector multiply
    procedure :: scalar_product => scalar_product_omp      !! Scalar product
    procedure :: field_max_mean => field_max_mean_omp      !! Compute max and mean
    procedure :: field_scale => field_scale_omp            !! Scale field
    procedure :: field_shift => field_shift_omp            !! Shift field values
    procedure :: field_set_face => field_set_face_omp      !! Set face values
    procedure :: field_volume_integral => field_volume_integral_omp !! Volume integral
    procedure :: copy_data_to_f => copy_data_to_f_omp      !! Copy data to field
    procedure :: copy_f_to_data => copy_f_to_data_omp      !! Copy field to data
    procedure :: init_poisson_fft => init_omp_poisson_fft  !! Initialise FFT Poisson
    procedure :: transeq_omp_dist                          !! Distributed transeq (internal)
  end type omp_backend_t

  interface omp_backend_t
    module procedure init
  end interface omp_backend_t

contains

  function init(mesh, allocator) result(backend)
    !! Initialise OpenMP backend with mesh and allocator.
    !!
    !! Sets up the backend by:
    !! 1. Calling base initialisation
    !! 2. Linking mesh and allocator
    !! 3. Determining maximum number of groups across directions
    !! 4. Allocating communication buffers for halo exchange
    !!
    !! **Buffer Sizing:** Buffers are sized based on the largest decomposition
    !! direction to handle all reordering operations efficiently.
    implicit none

    type(mesh_t), target, intent(inout) :: mesh         !! Mesh with decomposition
    class(allocator_t), target, intent(inout) :: allocator !! Memory allocator
    type(omp_backend_t) :: backend                       !! Initialised backend instance

    integer :: n_groups

    call backend%base_init()

    select type (allocator)
    type is (allocator_t)
      ! class level access to the allocator
      backend%allocator => allocator
    end select

    backend%mesh => mesh
    n_groups = maxval([backend%allocator%get_n_groups(DIR_X), &
                       backend%allocator%get_n_groups(DIR_Y), &
                       backend%allocator%get_n_groups(DIR_Z)])

    allocate (backend%u_send_s(SZ, backend%n_halo, n_groups))
    allocate (backend%u_send_e(SZ, backend%n_halo, n_groups))
    allocate (backend%u_recv_s(SZ, backend%n_halo, n_groups))
    allocate (backend%u_recv_e(SZ, backend%n_halo, n_groups))
    allocate (backend%v_send_s(SZ, backend%n_halo, n_groups))
    allocate (backend%v_send_e(SZ, backend%n_halo, n_groups))
    allocate (backend%v_recv_s(SZ, backend%n_halo, n_groups))
    allocate (backend%v_recv_e(SZ, backend%n_halo, n_groups))
    allocate (backend%w_send_s(SZ, backend%n_halo, n_groups))
    allocate (backend%w_send_e(SZ, backend%n_halo, n_groups))
    allocate (backend%w_recv_s(SZ, backend%n_halo, n_groups))
    allocate (backend%w_recv_e(SZ, backend%n_halo, n_groups))

    allocate (backend%du_send_s(SZ, 1, n_groups))
    allocate (backend%du_send_e(SZ, 1, n_groups))
    allocate (backend%du_recv_s(SZ, 1, n_groups))
    allocate (backend%du_recv_e(SZ, 1, n_groups))
    allocate (backend%dud_send_s(SZ, 1, n_groups))
    allocate (backend%dud_send_e(SZ, 1, n_groups))
    allocate (backend%dud_recv_s(SZ, 1, n_groups))
    allocate (backend%dud_recv_e(SZ, 1, n_groups))
    allocate (backend%d2u_send_s(SZ, 1, n_groups))
    allocate (backend%d2u_send_e(SZ, 1, n_groups))
    allocate (backend%d2u_recv_s(SZ, 1, n_groups))
    allocate (backend%d2u_recv_e(SZ, 1, n_groups))

  end function init

  subroutine alloc_omp_tdsops( &
    self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    )
    !! Allocate and initialise tridiagonal operator for OMP backend.
    !!
    !! Creates a `tdsops_t` object configured for the specified operation
    !! (derivative, interpolation) with chosen compact scheme and boundary
    !! conditions. Handles grid stretching and viscous corrections.
    implicit none

    class(omp_backend_t) :: self                        !! Backend instance
    class(tdsops_t), allocatable, intent(inout) :: tdsops !! Tridiagonal operator to allocate
    integer, intent(in) :: n_tds                        !! Number of points in direction
    real(dp), intent(in) :: delta                       !! Grid spacing
    character(*), intent(in) :: operation, scheme       !! Operation type and scheme name
    integer, intent(in) :: bc_start, bc_end             !! Boundary condition codes
    real(dp), optional, intent(in) :: stretch(:), stretch_correct(:) !! Grid stretching
    integer, optional, intent(in) :: n_halo             !! Halo depth
    character(*), optional, intent(in) :: from_to       !! Data location transition
    logical, optional, intent(in) :: sym                !! Symmetry flag
    real(dp), optional, intent(in) :: c_nu, nu0_nu      !! Viscous correction parameters

    allocate (tdsops_t :: tdsops)

    select type (tdsops)
    type is (tdsops_t)
      tdsops = tdsops_t(n_tds, delta, operation, scheme, bc_start, bc_end, &
                        stretch, stretch_correct, n_halo, from_to, sym, &
                        c_nu, nu0_nu)
    end select

  end subroutine alloc_omp_tdsops

  subroutine transeq_x_omp(self, du, dv, dw, u, v, w, nu, dirps)
    !! Compute transport equation RHS in X direction.
    !!
    !! Evaluates convection and diffusion terms for momentum equations:
    !! \( du/dt = -u \cdot \nabla u + \nu \nabla^2 u \)
    !!
    !! Delegates to `transeq_omp_dist` which handles halo exchange and
    !! distributed compact schemes.
    implicit none

    class(omp_backend_t) :: self                  !! Backend instance
    class(field_t), intent(inout) :: du, dv, dw  !! Output: velocity RHS
    class(field_t), intent(in) :: u, v, w        !! Input: velocity fields
    real(dp), intent(in) :: nu                    !! Kinematic viscosity
    type(dirps_t), intent(in) :: dirps            !! Directional operators

    call self%transeq_omp_dist(du, dv, dw, u, v, w, nu, dirps)

  end subroutine transeq_x_omp

  subroutine transeq_y_omp(self, du, dv, dw, u, v, w, nu, dirps)
    !! Compute transport equation RHS in Y direction.
    !!
    !! Calculates convective and viscous terms for Y-pencil decomposition.
    !! Velocity components are reordered (v, u, w) to align primary
    !! direction with pencil orientation before calling distributed kernel.
    !!
    !! See [[transeq_x_omp]] for transport equation formulation.
    implicit none

    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(inout) :: du, dv, dw  !! Time derivatives (output)
    class(field_t), intent(in) :: u, v, w      !! Velocity components
    real(dp), intent(in) :: nu                 !! Kinematic viscosity
    type(dirps_t), intent(in) :: dirps         !! Spectral operators

    ! u, v, w is reordered so that we pass v, u, w
    call self%transeq_omp_dist(dv, du, dw, v, u, w, nu, dirps)

  end subroutine transeq_y_omp

  subroutine transeq_z_omp(self, du, dv, dw, u, v, w, nu, dirps)
    !! Compute transport equation RHS in Z direction.
    !!
    !! Calculates convective and viscous terms for Z-pencil decomposition.
    !! Velocity components are reordered (w, u, v) to align primary
    !! direction with pencil orientation before calling distributed kernel.
    !!
    !! See [[transeq_x_omp]] for transport equation formulation.
    implicit none

    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(inout) :: du, dv, dw  !! Time derivatives (output)
    class(field_t), intent(in) :: u, v, w      !! Velocity components
    real(dp), intent(in) :: nu                 !! Kinematic viscosity
    type(dirps_t), intent(in) :: dirps         !! Spectral operators

    ! u, v, w is reordered so that we pass w, u, v
    call self%transeq_omp_dist(dw, du, dv, w, u, v, nu, dirps)

  end subroutine transeq_z_omp

  subroutine transeq_species_omp(self, dspec, uvw, spec, nu, dirps, sync)
    !! Compute transport equation RHS for scalar species.
    !!
    !! Calculates convective and diffusive terms for a passive scalar
    !! (temperature, concentration, etc.) transported by velocity field.
    !!
    !! **Equation:** `$\partial\phi/\partial t = -\mathbf{u}\cdot\nabla\phi + \nu\nabla^2\phi$` where $\phi$ is the scalar species.
    !!
    !! **Synchronisation:** When `sync=.true.`, performs halo exchange
    !! for velocity field before computation. Always exchanges scalar halos.
    implicit none

    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(inout) :: dspec     !! Time derivative of species (output)
    class(field_t), intent(in) :: uvw          !! Velocity component in pencil direction
    class(field_t), intent(in) :: spec         !! Species concentration/temperature
    real(dp), intent(in) :: nu                 !! Diffusivity coefficient
    type(dirps_t), intent(in) :: dirps         !! Spectral operators
    logical, intent(in) :: sync                !! Perform velocity halo exchange if true

    integer :: n_groups

    n_groups = self%allocator%get_n_groups(dirps%dir)

    ! Halo exchange for momentum if needed
    if (sync) then
      call copy_into_buffers(self%u_send_s, self%u_send_e, uvw%data, &
                             dirps%der1st%n_tds, n_groups)
      call sendrecv_fields(self%u_recv_s, self%u_recv_e, &
                           self%u_send_s, self%u_send_e, &
                           SZ*self%n_halo*n_groups, &
                           self%mesh%par%nproc_dir(dirps%dir), &
                           self%mesh%par%pprev(dirps%dir), &
                           self%mesh%par%pnext(dirps%dir))
    end if

    ! Halo exchange for the given field
    call copy_into_buffers(self%v_send_s, self%v_send_e, spec%data, &
                           dirps%der1st%n_tds, n_groups)
    call sendrecv_fields(self%v_recv_s, self%v_recv_e, &
                         self%v_send_s, self%v_send_e, &
                         SZ*self%n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dirps%dir), &
                         self%mesh%par%pprev(dirps%dir), &
                         self%mesh%par%pnext(dirps%dir))

    ! combine convection and diffusion
    call transeq_dist_component(self, dspec, spec, uvw, nu, &
                                self%v_recv_s, self%v_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st, dirps%der1st_sym, &
                                dirps%der2nd, dirps%dir)

  end subroutine transeq_species_omp

  subroutine transeq_omp_dist(self, du, dv, dw, u, v, w, nu, dirps)
    !! Internal: Distributed transport equation implementation.
    !!
    !! Orchestrates the complete transport equation calculation for
    !! all velocity components. First performs halo exchange for
    !! distributed compact derivatives, then computes each component's
    !! RHS using transeq_dist_component.
    !!
    !! **Called by:** transeq_x/y/z_omp after velocity reordering
    implicit none

    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(inout) :: du, dv, dw  !! Time derivatives (output)
    class(field_t), intent(in) :: u, v, w      !! Velocity components (reordered for pencil direction)
    real(dp), intent(in) :: nu                 !! Kinematic viscosity
    type(dirps_t), intent(in) :: dirps         !! Spectral operators

    call transeq_halo_exchange(self, u, v, w, dirps%dir)

    call transeq_dist_component(self, du, u, u, nu, &
                                self%u_recv_s, self%u_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st, dirps%der1st_sym, &
                                dirps%der2nd, dirps%dir)
    call transeq_dist_component(self, dv, v, u, nu, &
                                self%v_recv_s, self%v_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st_sym, dirps%der1st, &
                                dirps%der2nd_sym, dirps%dir)
    call transeq_dist_component(self, dw, w, u, nu, &
                                self%w_recv_s, self%w_recv_e, &
                                self%u_recv_s, self%u_recv_e, &
                                dirps%der1st_sym, dirps%der1st, &
                                dirps%der2nd_sym, dirps%dir)

  end subroutine transeq_omp_dist

  subroutine transeq_halo_exchange(self, u, v, w, dir)
    !! Internal: Perform halo exchange for all velocity components.
    !!
    !! Exchanges 4-point halos between neighbouring MPI processes for
    !! distributed compact finite difference stencils. Copies boundary
    !! data into send buffers, performs MPI sendrecv, stores in receive
    !! buffers for use in derivative calculations.
    !!
    !! **Operation:** Copy to buffers $\rightarrow$ MPI_Sendrecv $\rightarrow$ Store halos
    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(in) :: u, v, w      !! Velocity components
    integer, intent(in) :: dir                 !! Communication direction
    integer :: n, nproc_dir, pprev, pnext
    integer :: n_groups

    n_groups = self%allocator%get_n_groups(dir)
    n = self%mesh%get_n(u)
    nproc_dir = self%mesh%par%nproc_dir(dir)
    pprev = self%mesh%par%pprev(dir)
    pnext = self%mesh%par%pnext(dir)

    call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, &
                           n, n_groups)
    call copy_into_buffers(self%v_send_s, self%v_send_e, v%data, &
                           n, n_groups)
    call copy_into_buffers(self%w_send_s, self%w_send_e, w%data, &
                           n, n_groups)

    call sendrecv_fields(self%u_recv_s, self%u_recv_e, &
                         self%u_send_s, self%u_send_e, &
                         SZ*self%n_halo*n_groups, &
                         nproc_dir, pprev, pnext)
    call sendrecv_fields(self%v_recv_s, self%v_recv_e, &
                         self%v_send_s, self%v_send_e, &
                         SZ*self%n_halo*n_groups, &
                         nproc_dir, pprev, pnext)
    call sendrecv_fields(self%w_recv_s, self%w_recv_e, &
                         self%w_send_s, self%w_send_e, &
                         SZ*self%n_halo*n_groups, &
                         nproc_dir, pprev, pnext)

  end subroutine transeq_halo_exchange

  subroutine transeq_dist_component(self, rhs_du, u, conv, nu, &
                                    u_recv_s, u_recv_e, &
                                    conv_recv_s, conv_recv_e, &
                                    tdsops_du, tdsops_dud, tdsops_d2u, dir)
    !! Internal: Compute single component of transport equation RHS.
    !!
    !! Calculates RHS for one velocity component using skew-symmetric form:
    !!
    !! **Formula:** `rhs = -0.5*(conv*du/dx + d(u*conv)/dx) + nu*d2u/dx2`
    !!
    !! Uses distributed compact FD kernels with halo data from neighbours.
    !! Allocates temporary storage for derivatives and releases after use.
    !!
    !! **Skew-symmetric:** Reduces aliasing errors in nonlinear convection.
    class(omp_backend_t) :: self               !! Backend instance
    !> The result field, it is also used as temporary storage
    class(field_t), intent(inout) :: rhs_du   !! RHS output (also temp storage)
    class(field_t), intent(in) :: u, conv     !! Velocity component and convecting velocity
    real(dp), intent(in) :: nu                !! Kinematic viscosity
    real(dp), dimension(:, :, :), intent(in) :: u_recv_s, u_recv_e, &
                                                conv_recv_s, conv_recv_e  !! Halo data from neighbours
    class(tdsops_t), intent(in) :: tdsops_du      !! First derivative operator
    class(tdsops_t), intent(in) :: tdsops_dud     !! Product derivative operator
    class(tdsops_t), intent(in) :: tdsops_d2u     !! Second derivative operator
    integer, intent(in) :: dir                    !! Direction index
    class(field_t), pointer :: d2u, dud

    dud => self%allocator%get_block(dir)
    d2u => self%allocator%get_block(dir)

    call exec_dist_transeq_compact( &
      rhs_du%data, dud%data, d2u%data, &
      self%du_send_s, self%du_send_e, self%du_recv_s, self%du_recv_e, &
      self%dud_send_s, self%dud_send_e, self%dud_recv_s, self%dud_recv_e, &
      self%d2u_send_s, self%d2u_send_e, self%d2u_recv_s, self%d2u_recv_e, &
      u%data, u_recv_s, u_recv_e, &
      conv%data, conv_recv_s, conv_recv_e, &
      tdsops_du, tdsops_dud, tdsops_d2u, nu, &
      self%mesh%par%nproc_dir(dir), self%mesh%par%pprev(dir), &
      self%mesh%par%pnext(dir), self%allocator%get_n_groups(dir))

    call rhs_du%set_data_loc(u%data_loc)

    call self%allocator%release_block(dud)
    call self%allocator%release_block(d2u)

  end subroutine transeq_dist_component

  subroutine tds_solve_omp(self, du, u, tdsops)
    !! Solve tridiagonal system for compact finite difference operation.
    !!
    !! Applies compact scheme operator to field using Thomas algorithm.
    !! Handles both local (single-process) and distributed (multi-process)
    !! solves depending on decomposition configuration.
    !!
    !! **Data Location:** Updates output data location based on operator's
    !! `move` specification (e.g., CELL to VERT for interpolation).
    implicit none

    class(omp_backend_t) :: self              !! Backend instance
    class(field_t), intent(inout) :: du      !! Output field
    class(field_t), intent(in) :: u          !! Input field
    class(tdsops_t), intent(in) :: tdsops    !! Tridiagonal operator

    ! Check if direction matches for both in/out fields
    if (u%dir /= du%dir) then
      error stop 'DIR mismatch between fields in tds_solve.'
    end if

    if (u%data_loc /= NULL_LOC) then
      call du%set_data_loc(move_data_loc(u%data_loc, u%dir, tdsops%move))
    end if

    call tds_solve_dist(self, du, u, tdsops)

  end subroutine tds_solve_omp

  subroutine tds_solve_dist(self, du, u, tdsops)
    !! Internal: Distributed tridiagonal solve with halo exchange.
    !!
    !! Solves compact finite difference system across multiple MPI processes.
    !! Performs halo exchange before calling distributed Thomas algorithm
    !! kernel. Used when domain decomposition splits the pencil direction.
    !!
    !! **Algorithm:**
    !! 1. Copy boundary data into send buffers
    !! 2. MPI_Sendrecv for halo exchange
    !! 3. Distributed Thomas algorithm with boundary coupling
    implicit none

    class(omp_backend_t) :: self               !! Backend instance
    class(field_t), intent(inout) :: du       !! Solution field (output)
    class(field_t), intent(in) :: u           !! RHS field
    class(tdsops_t), intent(in) :: tdsops     !! Tridiagonal operator
    integer :: n_groups, dir

    dir = u%dir
    n_groups = self%allocator%get_n_groups(dir)

    call copy_into_buffers(self%u_send_s, self%u_send_e, u%data, &
                           tdsops%n_tds, n_groups)

    ! halo exchange
    call sendrecv_fields(self%u_recv_s, self%u_recv_e, &
                         self%u_send_s, self%u_send_e, &
                         SZ*self%n_halo*n_groups, &
                         self%mesh%par%nproc_dir(dir), &
                         self%mesh%par%pprev(dir), &
                         self%mesh%par%pnext(dir))

    call exec_dist_tds_compact( &
      du%data, u%data, self%u_recv_s, self%u_recv_e, &
      self%du_send_s, self%du_send_e, self%du_recv_s, self%du_recv_e, &
      tdsops, self%mesh%par%nproc_dir(dir), &
      self%mesh%par%pprev(dir), self%mesh%par%pnext(dir), &
      n_groups)

  end subroutine tds_solve_dist

  subroutine reorder_omp(self, u_, u, direction)
    !! Reorder field data between different pencil decompositions.
    !!
    !! Transforms field layout from one decomposition direction to another
    !! (e.g., X-pencils to Y-pencils). Uses MPI All-to-All communication
    !! to redistribute data across processes.
    !!
    !! **Directions:** DIR_X, DIR_Y, DIR_Z specify pencil orientations.
    !! Each pencil is contiguous along its direction and distributed in
    !! the other two dimensions.
    !!
    !! **Performance:** Critical operation for multi-dimensional algorithms.
    !! Uses `get_index_reordering` for efficient cache-friendly reordering.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: u_   !! Output field (reordered)
    class(field_t), intent(in) :: u       !! Input field
    integer, intent(in) :: direction      !! Reordering direction code
    integer, dimension(3) :: dims, cart_padded
    integer :: i, j, k
    integer :: out_i, out_j, out_k
    integer :: dir_from, dir_to

    dims = self%allocator%get_padded_dims(u%dir)
    cart_padded = self%allocator%get_padded_dims(DIR_C)
    call get_dirs_from_rdr(dir_from, dir_to, direction)

    !$omp parallel do private(out_i, out_j, out_k) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u_%data(out_i, out_j, out_k) = u%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

    ! reorder keeps the data_loc the same
    call u_%set_data_loc(u%data_loc)

  end subroutine reorder_omp

  subroutine sum_yintox_omp(self, u, u_)
    !! Sum Y-pencils into X-pencils through reordering.
    !!
    !! Performs directional reduction by reordering from Y to X pencils
    !! and summing the result into the destination field. Useful for
    !! integrating quantities along the Y direction.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: u    !! Destination field (X-pencils, accumulates result)
    class(field_t), intent(in) :: u_      !! Source field (Y-pencils)

    call sum_intox_omp(self, u, u_, DIR_Y)

  end subroutine sum_yintox_omp

  subroutine sum_zintox_omp(self, u, u_)
    !! Sum Z-pencils into X-pencils through reordering.
    !!
    !! Performs directional reduction by reordering from Z to X pencils
    !! and summing the result into the destination field. Useful for
    !! integrating quantities along the Z direction.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: u    !! Destination field (X-pencils, accumulates result)
    class(field_t), intent(in) :: u_      !! Source field (Z-pencils)

    call sum_intox_omp(self, u, u_, DIR_Z)

  end subroutine sum_zintox_omp

  subroutine sum_intox_omp(self, u, u_, dir_to)
    !! Internal helper: Sum reordered field into X-pencils.
    !!
    !! Reorders source field from X-pencils to specified direction,
    !! then accumulates into destination field. Called by sum_yintox_omp
    !! and sum_zintox_omp for directional integration.
    !!
    !! **Algorithm:** Reorder with index mapping, accumulate with +=

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: u    !! Destination field (accumulates result)
    class(field_t), intent(in) :: u_      !! Source field
    integer, intent(in) :: dir_to         !! Target direction (DIR_Y or DIR_Z)

    integer :: dir_from
    integer, dimension(3) :: dims, cart_padded
    integer :: i, j, k    ! Working indices
    integer :: ii, jj, kk ! Transpose indices

    dir_from = DIR_X

    dims = self%allocator%get_padded_dims(u%dir)
    cart_padded = self%allocator%get_padded_dims(DIR_C)

    !$omp parallel do private(i, ii, jj, kk) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(ii, jj, kk, i, j, k, &
                                    dir_from, dir_to, SZ, cart_padded)
          u%data(i, j, k) = u%data(i, j, k) + u_%data(ii, jj, kk)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine sum_intox_omp

  subroutine veccopy_omp(self, dst, src)
    !! Copy field data from source to destination.
    !!
    !! Element-wise copy with OpenMP parallelisation. Both fields
    !! must have the same decomposition direction and dimensions.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: dst  !! Destination field
    class(field_t), intent(in) :: src     !! Source field
    integer :: i, j, k                     !! Loop indices

    if (src%dir /= dst%dir) then
      error stop "Called vector copy with incompatible fields"
    end if

    if (dst%dir == DIR_C) error stop 'veccopy does not support DIR_C fields'

    !$omp parallel do
    do k = 1, size(dst%data, 3)
      do j = 1, size(dst%data, 2)
        !$omp simd
        do i = 1, SZ
          dst%data(i, j, k) = src%data(i, j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine veccopy_omp

  subroutine vecadd_omp(self, a, x, b, y)
    !! Vector addition: y = a*x + b*y (in-place AXPBY).
    !!
    !! Scaled in-place vector addition with OpenMP parallelisation
    !! and SIMD vectorisation. Implements the BLAS AXPBY operation.
    !!
    !! **Formula:** `y := a*x + b*y` where a, b are scalars.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    real(dp), intent(in) :: a              !! Scalar multiplier for x
    class(field_t), intent(in) :: x        !! First input field
    real(dp), intent(in) :: b              !! Scalar multiplier for y
    class(field_t), intent(inout) :: y     !! Second input field (overwritten with result)
    integer :: i, j, k

    if (x%dir /= y%dir) then
      error stop "Called vector add with incompatible fields"
    end if

    if (y%dir == DIR_C) error stop 'vecadd does not support DIR_C fields'

    !$omp parallel do
    do k = 1, size(y%data, 3)
      do j = 1, size(y%data, 2)
        !$omp simd
        do i = 1, SZ
          y%data(i, j, k) = a*x%data(i, j, k) + b*y%data(i, j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine vecadd_omp

  subroutine vecmult_omp(self, y, x)
    !! Element-wise multiplication: y = y * x (in-place).
    !!
    !! In-place element-wise multiplication with OpenMP parallelisation
    !! and SIMD vectorisation. Often used for applying masks or
    !! multiplying solution components.
    !! [[m_base_backend(module):vecmult(interface)]]
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: y    !! Field to multiply and store result
    class(field_t), intent(in) :: x        !! Multiplier field
    integer :: i, j, k                     !! Loop indices

    if (x%dir /= y%dir) then
      error stop "Called vector multiply with incompatible fields"
    end if

    if (y%dir == DIR_C) error stop 'vecmult does not support DIR_C fields'

    !$omp parallel do
    do k = 1, size(y%data, 3)
      do j = 1, size(y%data, 2)
        !$omp simd
        do i = 1, SZ
          y%data(i, j, k) = y%data(i, j, k)*x%data(i, j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine vecmult_omp

  real(dp) function scalar_product_omp(self, x, y) result(s)
    !! Compute global scalar product (dot product) of two fields.
    !!
    !! Calculates the dot product $\sum(x_i \times y_i)$ across all grid points
    !! and all MPI processes. Uses OpenMP parallelisation with reduction
    !! and MPI_Allreduce for global sum.
    !!
    !! **Algorithm:** Local parallel reduction $\rightarrow$ MPI_Allreduce
    !! **Data location:** Both fields must be at the same location (CELL/VERT).
    !! [[m_base_backend(module):scalar_product(interface)]]
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(in) :: x, y     !! Input fields
    class(field_t), pointer :: x_, y_      !! Pointers for data access
    integer, dimension(3) :: dims          !! Field dimensions
    integer :: i, j, k, ii                 !! Loop indices
    integer :: nvec, remstart              !! Vectorisation variables
    integer :: ierr                        !! MPI error code

    if ((x%data_loc == NULL_LOC) .or. (y%data_loc == NULL_LOC)) then
      error stop "You must set the data_loc before calling scalar product"
    end if
    if (x%data_loc /= y%data_loc) then
      error stop "Called scalar product with incompatible fields"
    end if

    ! Reorient data into temporary DIR_C storage
    x_ => self%allocator%get_block(DIR_C, x%data_loc)
    call self%get_field_data(x_%data, x)
    y_ => self%allocator%get_block(DIR_C, y%data_loc)
    call self%get_field_data(y_%data, y)

    dims = self%mesh%get_dims(x_%data_loc)

    nvec = dims(1)/SZ
    remstart = nvec*SZ + 1

    s = 0.0_dp
    !$omp parallel do reduction(+:s) private(i, ii) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        ! Execute inner vectorised loops
        do ii = 1, nvec
          !$omp simd reduction(+:s)
          do i = 1, SZ
            s = s + x_%data(i + (ii - 1)*SZ, j, k)* &
                y_%data(i + (ii - 1)*SZ, j, k)
          end do
          !$omp end simd
        end do

        ! Remainder loop
        do i = remstart, dims(1)
          s = s + x_%data(i, j, k)*y_%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

    ! Release temporary storage
    call self%allocator%release_block(x_)
    call self%allocator%release_block(y_)

    ! Reduce the result
    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, &
                       ierr)

  end function scalar_product_omp

  subroutine copy_into_buffers(u_send_s, u_send_e, u, n, n_groups)
    !! Internal helper: Copy halo data into send buffers.
    !!
    !! Extracts 4-point halos from start and end of domain for
    !! MPI communication. Used in transeq_halo_exchange to prepare
    !! boundary data for neighbour processes.
    !!
    !! **Buffer layout:** (SZ, 4, n_groups) for cache efficiency
    implicit none

    real(dp), dimension(:, :, :), intent(out) :: u_send_s  !! Send buffer for start boundary
    real(dp), dimension(:, :, :), intent(out) :: u_send_e  !! Send buffer for end boundary
    real(dp), dimension(:, :, :), intent(in) :: u          !! Field data
    integer, intent(in) :: n                               !! Domain size in communication direction
    integer, intent(in) :: n_groups                        !! Number of pencil groups
    integer :: i, j, k                                     !! Loop indices
    integer :: n_halo = 4                                  !! Halo width (compact scheme stencil)

    !$omp parallel do
    do k = 1, n_groups
      do j = 1, n_halo
        !$omp simd
        do i = 1, SZ
          u_send_s(i, j, k) = u(i, j, k)
          u_send_e(i, j, k) = u(i, n - n_halo + j, k)
        end do
        !$omp end simd
      end do
    end do
    !$omp end parallel do

  end subroutine copy_into_buffers

  subroutine field_max_mean_omp(self, max_val, mean_val, f, enforced_data_loc)
    !! Compute global maximum and mean of a field.
    !!
    !! Calculates maximum and mean values across all grid points and
    !! MPI processes. Uses data location (CELL/VERT) to determine
    !! valid domain extents, excluding padding and ghost cells.
    !!
    !! **Algorithm:**
    !! 1. Local parallel max/sum reduction with OpenMP
    !! 2. MPI_Allreduce for global max/sum
    !! 3. Mean = global_sum / global_count
    !!
    !! **Data location:** Can be enforced or read from field metadata.
    !! [[m_base_backend(module):field_max_mean(interface)]]
    implicit none

    class(omp_backend_t) :: self                   !! Backend instance
    real(dp), intent(out) :: max_val, mean_val     !! Global maximum and mean values
    class(field_t), intent(in) :: f                !! Input field
    integer, optional, intent(in) :: enforced_data_loc  !! Override data location if provided

    real(dp) :: val, max_p, sum_p, max_pncl, sum_pncl
    integer :: data_loc, dims(3), dims_padded(3), n, n_i, n_i_pad, n_j
    integer :: i, j, k, k_i, k_j, ierr

    if (f%data_loc == NULL_LOC .and. (.not. present(enforced_data_loc))) then
      error stop 'The input field to omp::field_max_mean does not have a &
                  &valid f%data_loc. You may enforce a data_loc of your &
                  &choice as last argument to carry on at your own risk!'
    end if

    if (present(enforced_data_loc)) then
      data_loc = enforced_data_loc
    else
      data_loc = f%data_loc
    end if

    dims = self%mesh%get_dims(data_loc)
    dims_padded = self%allocator%get_padded_dims(DIR_C)

    if (f%dir == DIR_X) then
      n = dims(1); n_j = dims(2); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Y) then
      n = dims(2); n_j = dims(1); n_i = dims(3); n_i_pad = dims_padded(3)
    else if (f%dir == DIR_Z) then
      n = dims(3); n_j = dims(1); n_i = dims(2); n_i_pad = dims_padded(2)
    else
      error stop 'field_max_mean does not support DIR_C fields!'
    end if

    sum_p = 0._dp
    max_p = 0._dp
    !$omp parallel do collapse(2) reduction(+:sum_p) reduction(max:max_p) &
    !$omp private(k, val, sum_pncl, max_pncl)
    do k_j = 1, (n_j - 1)/SZ + 1 ! loop over stacked groups
      do k_i = 1, n_i
        k = k_j + (k_i - 1)*((n_j - 1)/SZ + 1)
        sum_pncl = 0._dp
        max_pncl = 0._dp
        do j = 1, n
          ! loop over only non-padded entries in the present group
          do i = 1, min(SZ, n_j - (k_j - 1)*SZ)
            val = abs(f%data(i, j, k))
            sum_pncl = sum_pncl + val
            max_pncl = max(max_pncl, val)
          end do
        end do
        sum_p = sum_p + sum_pncl
        max_p = max(max_p, max_pncl)
      end do
    end do
    !$omp end parallel do

    ! rank-local values
    max_val = max_p
    mean_val = sum_p/product(self%mesh%get_global_dims(data_loc))

    ! make sure all ranks have final values
    call MPI_Allreduce(MPI_IN_PLACE, max_val, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mean_val, 1, MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine field_max_mean_omp

  subroutine field_scale_omp(self, f, a)
    !! Scale field by constant: f = a * f.
    !!
    !! Multiplies all field values by scalar a in-place.
    !! Uses Fortran array syntax for simplicity.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(in) :: f        !! Field to scale (modified in-place)
    real(dp), intent(in) :: a              !! Scaling factor

    f%data = a*f%data
  end subroutine field_scale_omp

  subroutine field_shift_omp(self, f, a)
    !! Shift field by constant: f = f + a.
    !!
    !! Adds scalar a to all field values in-place.
    !! Uses Fortran array syntax for simplicity.
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(in) :: f        !! Field to shift (modified in-place)
    real(dp), intent(in) :: a              !! Shift amount

    f%data = f%data + a
  end subroutine field_shift_omp

  subroutine field_set_face_omp(self, f, c_start, c_end, face)
    !! Set boundary face values to specified constants.
    !!
    !! Sets values on a specified domain face (X/Y/Z start/end)
    !! to given constants. Used for boundary condition enforcement.
    !!
    !! **Faces:** VERT_START_FACE, VERT_END_FACE, etc.
    !! [[m_base_backend(module):field_set_face(subroutine)]]
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(inout) :: f     !! Field to modify
    real(dp), intent(in) :: c_start        !! Value for start side of face
    real(dp), intent(in) :: c_end          !! Value for end side of face
    integer, intent(in) :: face            !! Face identifier constant

    integer :: dims(3), k, j, i_mod, k_end

    if (f%dir /= DIR_X) then
      error stop 'Setting a field face is only supported for DIR_X fields.'
    end if

    if (f%data_loc == NULL_LOC) then
      error stop 'field_set_face require a valid data_loc.'
    end if

    dims = self%mesh%get_dims(f%data_loc)

    select case (face)
    case (X_FACE)
      error stop 'Setting X_FACE is not yet supported.'
    case (Y_FACE)
      i_mod = mod(dims(2) - 1, SZ) + 1
      !$omp parallel do private(k_end)
      do k = 1, dims(3)
        k_end = k + (dims(2) - 1)/SZ*dims(3)
        do j = 1, dims(1)
          f%data(1, j, k) = c_start
          f%data(i_mod, j, k_end) = c_end
        end do
      end do
      !$omp end parallel do
    case (Z_FACE)
      error stop 'Setting Z_FACE is not yet supported.'
    case default
      error stop 'face is undefined.'
    end select

  end subroutine field_set_face_omp

  real(dp) function field_volume_integral_omp(self, f) result(s)
    !! Compute volume integral of field over domain.
    !!
    !! Calculates $\int f \,dV$ by summing all field values (at cell centres)
    !! and multiplying by grid cell volumes. Uses MPI_Allreduce for
    !! global sum across all processes.
    !!
    !! **Formula:** $\int f \,dV = \sum(f_i \times \Delta V_i)$ where $\Delta V$ from mesh
    !! **Assumption:** Field at cell centres (data_loc = CELL)
    implicit none

    class(omp_backend_t) :: self           !! Backend instance
    class(field_t), intent(in) :: f        !! Field to integrate

    real(dp) :: sum_p, sum_pncl
    integer :: dims(3), stacked, i, j, k, k_i, k_j, ierr

    if (f%data_loc == NULL_LOC) then
      error stop 'You must set the data_loc before calling volume integral.'
    end if
    if (f%dir /= DIR_X) then
      error stop 'Volume integral can only be called on DIR_X fields.'
    end if

    dims = self%mesh%get_dims(f%data_loc)
    stacked = (dims(2) - 1)/SZ + 1

    sum_p = 0._dp
    !$omp parallel do collapse(2) reduction(+:sum_p) private(k, sum_pncl)
    do k_j = 1, stacked ! loop over stacked groups
      do k_i = 1, dims(3)
        k = k_j + (k_i - 1)*stacked
        sum_pncl = 0._dp
        do j = 1, dims(1)
          ! loop over only non-padded entries in the present group
          do i = 1, min(SZ, dims(2) - (k_j - 1)*SZ)
            sum_pncl = sum_pncl + f%data(i, j, k)
          end do
        end do
        sum_p = sum_p + sum_pncl
      end do
    end do
    !$omp end parallel do

    ! rank-local values
    s = sum_p

    call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

  end function field_volume_integral_omp

  subroutine copy_data_to_f_omp(self, f, data)
    !! Copy raw array into field structure.
    !!
    !! Simple wrapper for field initialisation from external data.
    !! Uses Fortran array assignment for efficiency.
    class(omp_backend_t), intent(inout) :: self  !! Backend instance
    class(field_t), intent(inout) :: f           !! Target field
    real(dp), dimension(:, :, :), intent(in) :: data  !! Source data array

    f%data = data
  end subroutine copy_data_to_f_omp

  subroutine copy_f_to_data_omp(self, data, f)
    !! Copy field structure into raw array.
    !!
    !! Simple wrapper for field extraction to external data.
    !! Uses Fortran array assignment for efficiency.
    class(omp_backend_t), intent(inout) :: self  !! Backend instance
    real(dp), dimension(:, :, :), intent(out) :: data  !! Destination data array
    class(field_t), intent(in) :: f              !! Source field

    data = f%data
  end subroutine copy_f_to_data_omp

  subroutine init_omp_poisson_fft(self, mesh, xdirps, ydirps, zdirps, lowmem)
    !! Initialise FFT-based Poisson solver for OMP backend.
    !!
    !! Creates and configures omp_poisson_fft_t solver for pressure
    !! correction step. Uses 2DECOMP&FFT library for parallel FFTs
    !! in pencil decomposition.
    !!
    !! **Requirement:** WITH_2DECOMPFFT must be defined at compile time.
    !! **Low-memory mode:** Optional flag to reduce memory footprint.
#ifdef WITH_2DECOMPFFT
    use m_omp_poisson_fft, only: omp_poisson_fft_t
#endif

    implicit none

    class(omp_backend_t) :: self                   !! Backend instance
    type(mesh_t), intent(in) :: mesh               !! Mesh with grid spacing
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps  !! Spectral operators for each direction
    logical, optional, intent(in) :: lowmem        !! Enable low-memory mode

#ifdef WITH_2DECOMPFFT
    allocate (omp_poisson_fft_t :: self%poisson_fft)

    select type (poisson_fft => self%poisson_fft)
    type is (omp_poisson_fft_t)
      poisson_fft = omp_poisson_fft_t(mesh, xdirps, ydirps, zdirps, lowmem)
    end select
#else
    error stop 'This build does not support FFT based Poisson solver &
                &on the OpenMP backend!'
#endif

  end subroutine init_omp_poisson_fft

end module m_omp_backend

