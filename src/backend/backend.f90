module m_base_backend
  !! Abstract base backend defining the computational interface for X3D2 solver.
  !!
  !! This module defines the `base_backend_t` abstract type, which establishes
  !! the interface for all backend implementations (CUDA GPU, OpenMP CPU, etc.).
  !! The solver operates exclusively through these abstract interfaces, enabling
  !! complete architecture independence.
  !!
  !! **Architecture Pattern:**
  !!
  !! The backend abstraction follows the Strategy design pattern:
  !!
  !! - **Abstract interface** (`base_backend_t`): Defines deferred procedures for
  !!   all computational operations required by the solver
  !! - **Concrete implementations**: CUDA backend (`m_cuda_backend`) and OMP
  !!   backend (`m_omp_backend`) extend this base and provide architecture-specific
  !!   implementations
  !! - **Solver independence**: The solver (`m_solver`) calls backend methods
  !!   through the abstract interface without knowing the underlying implementation
  !!
  !! **Key Operations Defined:**
  !!
  !! - **Transport equation derivatives**: `transeq_x`, `transeq_y`, `transeq_z`
  !!   compute directional derivatives with halo exchange for distributed compact schemes
  !! - **Tridiagonal solves**: `tds_solve` applies compact finite difference operators
  !! - **Data reordering**: `reorder` transforms data between pencil decomposition
  !!   orientations (X, Y, Z directions)
  !! - **Field operations**: Vector arithmetic (`veccopy`, `vecadd`, `vecmult`),
  !!   reductions (`scalar_product`, `field_volume_integral`), and utilities
  !!   (`field_scale`, `field_shift`, `field_set_face`)
  !! - **Summation**: `sum_yintox`, `sum_zintox` for integrating fields along
  !!   specific directions
  !!
  !! **Backend Implementations:**
  !!
  !! - **CUDA backend** (`src/backend/cuda/backend.f90`): GPU-accelerated using
  !!   NVIDIA CUDA with device memory management and kernel launches
  !! - **OMP backend** (`src/backend/omp/backend.f90`): CPU parallelism via
  !!   OpenMP threading and MPI domain decomposition
  !!
  !! **Usage:**
  !!
  !! Backends are instantiated at runtime based on compile-time configuration and
  !! passed to the solver as a polymorphic pointer (`class(base_backend_t), pointer`).
  use mpi

  use m_allocator, only: allocator_t
  use m_common, only: dp, DIR_C, get_rdr_from_dirs
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: tdsops_t, dirps_t

  implicit none

  type, abstract :: base_backend_t
      !! Abstract base type defining the computational backend interface.
      !!
      !! This type encapsulates all architecture-specific operations required
      !! by the solver, enabling transparent execution on different hardware
      !! platforms (GPU via CUDA, CPU via OpenMP) without modifying solver code.
      !!
      !! **Design Philosophy:**
      !!
      !! The solver executes high-level operations (compute transport equation,
      !! solve tridiagonal systems, reorder data, etc.) through deferred procedures
      !! defined in this abstract interface. Each backend (CUDA, OMP) extends this
      !! type and implements these procedures using architecture-specific kernels,
      !! libraries, and memory management strategies.
      !!
      !! **Example Workflow:**
      !!
      !! When computing the transport equation, the solver calls:
      !!
      !! 1. `transeq_x`, `transeq_y`, `transeq_z` to compute directional derivatives
      !! 2. `reorder` to transform data between pencil orientations
      !! 3. `vecadd` to combine derivatives into divergence of \(U^*\)
      !!
      !! Each call dispatches to the appropriate backend implementation at runtime
      !! via dynamic polymorphism.
      !!
      !! **Components:**
      !!
      !! - `n_halo`: Number of halo layers for distributed compact schemes (fixed at 4)
      !! - `mesh`: Pointer to mesh object (grid dimensions, boundary conditions, decomposition)
      !! - `allocator`: Memory allocator for field storage (host for OMP, device for CUDA)
      !! - `poisson_fft`: FFT-based Poisson solver for pressure correction

    !> DistD2 implementation is hardcoded for 4 halo layers for all backends
    integer :: n_halo = 4
    type(mesh_t), pointer :: mesh
    class(allocator_t), pointer :: allocator
    class(poisson_fft_t), pointer :: poisson_fft
  contains
    procedure(transeq_ders), deferred :: transeq_x
    procedure(transeq_ders), deferred :: transeq_y
    procedure(transeq_ders), deferred :: transeq_z
    procedure(transeq_ders_spec), deferred :: transeq_species
    procedure(tds_solve), deferred :: tds_solve
    procedure(reorder), deferred :: reorder
    procedure(sum_intox), deferred :: sum_yintox
    procedure(sum_intox), deferred :: sum_zintox
    procedure(veccopy), deferred :: veccopy
    procedure(vecadd), deferred :: vecadd
    procedure(vecmult), deferred :: vecmult
    procedure(scalar_product), deferred :: scalar_product
    procedure(field_max_mean), deferred :: field_max_mean
    procedure(field_ops), deferred :: field_scale
    procedure(field_ops), deferred :: field_shift
    procedure(field_reduce), deferred :: field_volume_integral
    procedure(field_set_face), deferred :: field_set_face
    procedure(copy_data_to_f), deferred :: copy_data_to_f
    procedure(copy_f_to_data), deferred :: copy_f_to_data
    procedure(alloc_tdsops), deferred :: alloc_tdsops
    procedure(init_poisson_fft), deferred :: init_poisson_fft
    procedure :: base_init
    procedure :: get_field_data
    procedure :: set_field_data
  end type base_backend_t

  abstract interface
    subroutine transeq_ders(self, du, dv, dw, u, v, w, nu, dirps)
      !! Compute transport equation derivatives for velocity components.
      !!
      !! This is the core computational kernel for the transport equation,
      !! computing the advection-diffusion terms in one coordinate direction:
      !!
      !! \[
      !! \frac{\partial u_i}{\partial t} = -u \frac{\partial u_i}{\partial x_j}
      !!                                    - v \frac{\partial u_i}{\partial x_j}
      !!                                    - w \frac{\partial u_i}{\partial x_j}
      !!                                    + \nu \nabla^2 u_i
      !! \]
      !!
      !! (where the direction \(x_j\) is specified by `dirps`).
      !!
      !! **Runtime algorithm selection:**
      !!
      !! The exact algorithm used to obtain the derivatives is decided at runtime
      !! by the backend implementation. Backend implementations are responsible
      !! for directing calls to the appropriate algorithm based on:
      !!
      !! - Operator configuration in `dirps` (distributed vs local compact schemes)
      !! - Domain decomposition (number of processes in current direction)
      !! - Boundary conditions (periodic vs non-periodic)
      !!
      !! The implementation routes to either:
      !!
      !! - **Distributed algorithm** (`exec_dist_transeq_3fused`): For distributed
      !!   compact schemes with MPI halo exchange
      !! - **Thomas algorithm** (`exec_thom_transeq`): For localized/periodic operators
      import :: base_backend_t
      import :: field_t
      import :: dirps_t
      import :: dp
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw  !! Derivative outputs (momentum equation RHS)
      class(field_t), intent(in) :: u, v, w        !! Velocity components
      real(dp), intent(in) :: nu                   !! Kinematic viscosity
      type(dirps_t), intent(in) :: dirps           !! Directional derivative operators
    end subroutine transeq_ders
  end interface

  abstract interface
    subroutine transeq_ders_spec(self, dspec, uvw, spec, nu, dirps, sync)
      !! Compute transport equation derivatives for passive scalar species.
      !!
      !! Similar to `transeq_ders` but for passive scalar transport:
      !!
      !! \[
      !! \frac{\partial \phi}{\partial t} = -u \frac{\partial \phi}{\partial x_j}
      !!                                     - v \frac{\partial \phi}{\partial x_j}
      !!                                     - w \frac{\partial \phi}{\partial x_j}
      !!                                     + \nu \nabla^2 \phi
      !! \]
      !!
      !! where \(\phi\) is the scalar concentration and \(x_j\) is the direction
      !! specified by `dirps`.
      !!
      !! **Synchronization:**
      !!
      !! The `sync` flag controls whether to synchronize device-to-host memory
      !! transfers (CUDA backend) after computation. Set `.false.` when chaining
      !! multiple operations to avoid unnecessary transfers.
      import :: base_backend_t
      import :: field_t
      import :: dirps_t
      import :: dp
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: dspec  !! Scalar derivative output
      class(field_t), intent(in) :: uvw       !! Velocity component in current direction
      class(field_t), intent(in) :: spec      !! Scalar species concentration
      real(dp), intent(in) :: nu              !! Diffusion coefficient
      type(dirps_t), intent(in) :: dirps      !! Directional derivative operators
      logical, intent(in) :: sync             !! Synchronize device transfers (CUDA only)
    end subroutine transeq_ders_spec
  end interface

  abstract interface
    subroutine tds_solve(self, du, u, tdsops)
      !! Apply a tridiagonal operator to a field (compact finite difference operation).
      !!
      !! Solves the tridiagonal system arising from compact finite difference
      !! schemes:
      !!
      !! \[
      !! A f' = B f
      !! \]
      !!
      !! where \(A\) is the implicit (tridiagonal) operator, \(B\) is the explicit
      !! stencil, and \(f'\) is the derivative (or interpolated value).
      !!
      !! **Backend dispatch:**
      !!
      !! Routes to the appropriate tridiagonal solver:
      !!
      !! - **Distributed compact**: Uses `exec_dist_tds_compact` with MPI communication
      !!   for boundary coupling between processes
      !! - **Thomas algorithm**: Uses `exec_thom_tds_compact` for local/periodic systems
      !! - **GPU**: Uses batched tridiagonal solvers (cuSPARSE or custom kernels)
      !!
      !! **Operations supported:**
      !!
      !! First derivative, second derivative, interpolation, staggered derivatives
      !! (configured in `tdsops`).
      import :: base_backend_t
      import :: field_t
      import :: tdsops_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: du     !! Output field (derivative or interpolated values)
      class(field_t), intent(in) :: u         !! Input field
      class(tdsops_t), intent(in) :: tdsops   !! Tridiagonal operator (preprocessed)
    end subroutine tds_solve
  end interface

  abstract interface
    subroutine reorder(self, u_, u, direction)
      !! Reorder field data between pencil decomposition orientations.
      !!
      !! Transforms field layout from one pencil orientation to another to enable
      !! efficient tridiagonal solves in different coordinate directions:
      !!
      !! - **DIR_X**: X-pencils (data contiguous in X, decomposed in Y-Z)
      !! - **DIR_Y**: Y-pencils (data contiguous in Y, decomposed in X-Z)
      !! - **DIR_Z**: Z-pencils (data contiguous in Z, decomposed in X-Y)
      !! - **DIR_C**: Special compact orientation
      !!
      !! The `direction` parameter specifies the target orientation using reorder
      !! constants (`RDR_X2Y`, `RDR_Y2Z`, etc.).
      !!
      !! **Backend implementation:**
      !!
      !! - **CUDA**: GPU transpose kernels with coalesced memory access
      !! - **OMP**: MPI all-to-all communication with OpenMP threading
      !!
      !! **Performance note:** This is a bandwidth-intensive operation requiring
      !! global data movement (MPI or device memory transfers).
      import :: base_backend_t
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: u_      !! Output field (reordered)
      class(field_t), intent(in) :: u          !! Input field
      integer, intent(in) :: direction         !! Reorder direction (RDR_X2Y, RDR_Y2Z, etc.)
    end subroutine reorder
  end interface

  abstract interface
    subroutine sum_intox(self, u, u_)
      !! Sum directional derivatives back into X-oriented fields.
      !!
      !! Combines derivative contributions computed in different pencil orientations
      !! (Y-pencils, Z-pencils) back into the X-pencil orientation:
      !!
      !! \[
      !! u = u + u'
      !! \]
      !!
      !! This operation accumulates terms when computing composite derivatives
      !! like divergence:
      !!
      !! \[
      !! \nabla \cdot \mathbf{u} = \frac{\partial u}{\partial x}
      !!                          + \frac{\partial v}{\partial y}
      !!                          + \frac{\partial w}{\partial z}
      !! \]
      !!
      !! Each directional derivative is computed in its respective pencil orientation,
      !! then summed into X-pencils via `sum_yintox` and `sum_zintox`.
      !!
      !! **Note:** The input field `u_` must be in a Y or Z pencil orientation;
      !! the output `u` is always in X-pencil orientation.
      import :: base_backend_t
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: u   !! Accumulated field (X-pencils, updated in-place)
      class(field_t), intent(in) :: u_     !! Contribution to add (Y or Z pencils)
    end subroutine sum_intox
  end interface

  abstract interface
    subroutine veccopy(self, dst, src)
      !! Copy one field to another: `dst = src`.
      !!
      !! Performs an element-wise copy of all field data from `src` to `dst`.
      !! Both fields must have compatible dimensions and memory layout.
      !!
      !! **Backend implementation:**
      !!
      !! - **CUDA**: Device-to-device memory copy (cudaMemcpy)
      !! - **OMP**: Host memory copy (array assignment or memcpy)
      !!
      !! **Note:** This is a deep copy operation; the fields remain independent
      !! after the copy.
      import :: base_backend_t
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: dst  !! Destination field
      class(field_t), intent(in) :: src     !! Source field
    end subroutine veccopy
  end interface

  abstract interface
    subroutine vecadd(self, a, x, b, y)
      !! Compute linear combination of two fields (AXPBY operation).
      !!
      !! Performs the vector operation: \(y = a \cdot x + b \cdot y\)
      !!
      !! This is equivalent to the BLAS AXPBY operation, computing a scaled
      !! sum of two vectors. The result is stored in-place in `y`.
      !!
      !! **Common use cases:**
      !!
      !! - **Vector addition**: `vecadd(self, 1.0_dp, x, 1.0_dp, y)` \(\rightarrow\) \(y = x + y\)
      !! - **Scaled addition**: `vecadd(self, alpha, x, 1.0_dp, y)` \(\rightarrow\) \(y = \alpha x + y\)
      !! - **Replacement**: `vecadd(self, 1.0_dp, x, 0.0_dp, y)` \(\rightarrow\) \(y = x\)
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      real(dp), intent(in) :: a         !! Scaling factor for x
      class(field_t), intent(in) :: x   !! Input field
      real(dp), intent(in) :: b         !! Scaling factor for y
      class(field_t), intent(inout) :: y !! Input/output field (modified in-place)
    end subroutine vecadd
  end interface

  abstract interface
    subroutine vecmult(self, y, x)
      !! Element-wise (pointwise) multiplication of two fields.
      !!
      !! Performs the element-wise product: \(y = y \odot x\)
      !!
      !! Each element of `y` is multiplied by the corresponding element of `x`.
      !! The result is stored in-place in `y`. This is also known as the
      !! Hadamard product or pointwise multiplication.
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: y !! Input/output field (modified in-place)
      class(field_t), intent(in) :: x    !! Multiplier field
    end subroutine vecmult
  end interface

  abstract interface
    real(dp) function scalar_product(self, x, y) result(s)
      !! Compute the global scalar (dot) product of two fields.
      !!
      !! Calculates: \(s = \sum_{i} x_i \cdot y_i\)
      !!
      !! This computes the inner product (dot product) of two fields across
      !! all grid points. For distributed memory systems (MPI), partial sums
      !! from each process are accumulated via MPI reduction to produce the
      !! global sum.
      !!
      !! **Note:** The result includes contributions from all MPI ranks.
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(in) :: x !! First field
      class(field_t), intent(in) :: y !! Second field
    end function scalar_product
  end interface

  abstract interface
    subroutine field_ops(self, f, a)
      !! Generic interface for in-place field operations with a scalar constant.
      !!
      !! This abstract interface is implemented by two operations:
      !!
      !! - **field_scale**: Multiply field by constant: \(f = a \cdot f\)
      !! - **field_shift**: Add constant to field: \(f = f + a\)
      !!
      !! Both operations modify the field in-place and are backend-specific
      !! (GPU kernels for CUDA, array operations for OMP).
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(in) :: f  !! Field to operate on (modified in-place)
      real(dp), intent(in) :: a        !! Scalar constant (scaling factor or shift amount)
    end subroutine field_ops
  end interface

  abstract interface
    real(dp) function field_reduce(self, f) result(s)
      !! Reduce a field to a single scalar value via global summation.
      !!
      !! This abstract interface is currently implemented by:
      !!
      !! - **field_volume_integral**: Computes the volume integral \(\int f \,dV\)
      !!
      !! **Algorithm:**
      !!
      !! 1. **Local summation**: Each MPI process sums its local field values
      !!    (optionally weighted by cell volumes for volume integration)
      !! 2. **Global reduction**: MPI_Allreduce combines partial sums from all
      !!    processes to produce the global result
      !!
      !! **Backend implementations:**
      !!
      !! - **CUDA**: GPU reduction kernel followed by MPI_Allreduce
      !! - **OMP**: OpenMP parallel reduction followed by MPI_Allreduce
      !!
      !! **Requirements:**
      !!
      !! - Field must have `data_loc` set (cannot be `NULL_LOC`)
      !! - Field must be in X-pencil orientation (`dir = DIR_X`)
      !!
      !! **Use cases:**
      !!
      !! - Volume integrals for conservation checks
      !! - Global norms (L1, L2) for convergence monitoring
      !! - Total mass/energy calculations
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(in) :: f  !! Field to reduce
    end function field_reduce
  end interface

  abstract interface
    subroutine field_max_mean(self, max_val, mean_val, f, enforced_data_loc)
      !! Obtains maximum and mean values in a field
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      real(dp), intent(out) :: max_val, mean_val
      class(field_t), intent(in) :: f
      integer, optional, intent(in) :: enforced_data_loc
    end subroutine field_max_mean
  end interface

  abstract interface
    subroutine field_set_face(self, f, c_start, c_end, face)
      !! A field is a subdomain with a rectangular cuboid shape.
      !! It has 6 faces, and these faces are either a subdomain boundary
      !! or a global domain boundary based on the location of the subdomain.
      !! This subroutine allows us to set any of these faces to a value,
      !! 'c_start' and 'c_end' for faces at opposite sides.
      !! 'face' is one of `X_FACE`, `Y_FACE`, `Z_FACE` from `common.f90`
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: f
      real(dp), intent(in) :: c_start, c_end
      integer, intent(in) :: face
    end subroutine field_set_face
  end interface

  abstract interface
    subroutine copy_data_to_f(self, f, data)
         !! Copy the specialist data structure from device or host back
         !! to a regular 3D data array in host memory.
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t), intent(inout) :: self
      class(field_t), intent(inout) :: f
      real(dp), dimension(:, :, :), intent(in) :: data
    end subroutine copy_data_to_f

    subroutine copy_f_to_data(self, data, f)
         !! Copy a regular 3D array in host memory into the specialist
         !! data structure field that lives on device or host
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t), intent(inout) :: self
      real(dp), dimension(:, :, :), intent(out) :: data
      class(field_t), intent(in) :: f
    end subroutine copy_f_to_data
  end interface

  abstract interface
    subroutine alloc_tdsops( &
      !! Allocate and initialise a backend-specific tridiagonal operator.
      !!
      !! This deferred procedure creates a `tdsops_t` object configured for
      !! compact finite difference operations (derivatives, interpolation, etc.).
      !! The backend implementation allocates the appropriate subtype:
      !!
      !! - **CUDA backend**: Allocates `cuda_tdsops_t` with device memory pointers
      !!   for GPU execution
      !! - **OMP backend**: Allocates `omp_tdsops_t` with host memory for CPU execution
      !!
      !! The operator is fully preprocessed and ready for repeated application via
      !! `tds_solve`.
      !!
      !! **Required arguments:**
      !!
      !! - `n_tds`: System size (number of grid points in the operator direction)
      !! - `delta`: Grid spacing
      !! - `operation`: Operation type (`'first-deriv'`, `'second-deriv'`,
      !!   `'interpolate'`, `'stag-deriv'`)
      !! - `scheme`: Numerical scheme name (e.g., `'compact6'`, `'compact4'`)
      !! - `bc_start`, `bc_end`: Boundary condition flags (`BC_PERIODIC`,
      !!   `BC_NEUMANN`, `BC_DIRICHLET`)
      !!
      !! **Optional arguments:**
      !!
      !! - `stretch`: Stretching coefficients for non-uniform grids
      !! - `stretch_correct`: Correction for second derivatives on stretched grids
      !! - `n_halo`: Number of halo layers (default from backend)
      !! - `from_to`: Staggered grid direction (`'v2p'`, `'p2v'`)
      !! - `sym`: Field symmetry at Neumann boundaries (`.true.` = symmetric/even,
      !!   `.false.` = anti-symmetric/odd)
      !! - `c_nu`, `nu0_nu`: Hyperviscosity parameters for compact6-hyperviscous
      !!   second derivatives
      self, tdsops, n_tds, delta, operation, scheme, bc_start, bc_end, &
      stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
      )
      import :: base_backend_t
      import :: dp
      import :: tdsops_t
      implicit none

      class(base_backend_t) :: self
      class(tdsops_t), allocatable, intent(inout) :: tdsops
      integer, intent(in) :: n_tds
      real(dp), intent(in) :: delta
      character(*), intent(in) :: operation, scheme
      integer, intent(in) :: bc_start, bc_end
      real(dp), optional, intent(in) :: stretch(:), stretch_correct(:)
      integer, optional, intent(in) :: n_halo
      character(*), optional, intent(in) :: from_to
      logical, optional, intent(in) :: sym
      real(dp), optional, intent(in) :: c_nu, nu0_nu
    end subroutine alloc_tdsops
  end interface

  abstract interface
    subroutine init_poisson_fft(self, mesh, xdirps, ydirps, zdirps, lowmem)
      !! Initialise the backend-specific FFT-based Poisson solver.
      !!
      !! This deferred procedure creates and configures the Poisson solver object
      !! (`self%poisson_fft`) for solving the pressure Poisson equation:
      !! \(\nabla^2 \phi = f\)
      !!
      !! The backend implementation allocates the appropriate solver subtype:
      !!
      !! - **CUDA backend**: Allocates `cuda_poisson_fft_t` using cuFFT library
      !!   for GPU-accelerated FFT transforms
      !! - **OMP backend**: Allocates `omp_poisson_fft_t` using 2DECOMP&FFT library
      !!   for CPU FFT transforms with MPI parallelisation
      !!
      !! The solver requires directional derivative operators (`xdirps`, `ydirps`,
      !! `zdirps`) to construct spectral equivalence constants for handling:
      !!
      !! - Non-uniform grid spacing (stretching) in the Y-direction
      !! - Mixed boundary conditions (e.g., periodic in X/Z, Dirichlet in Y)
      !!
      !! **Arguments:**
      !!
      !! - `mesh`: Mesh object containing grid dimensions, boundary conditions,
      !!   and parallel decomposition information
      !! - `xdirps`, `ydirps`, `zdirps`: Second-derivative operators in each direction,
      !!   used to compute spectral equivalence constants for the modified wavenumbers
      !! - `lowmem` (optional): Low-memory mode flag. When `.true.`, reduces memory
      !!   footprint by deallocating temporary arrays after initialisation (CUDA only)
      !!
      !! **Note:** The Poisson solver is stored in `self%poisson_fft` and accessed
      !! by the solver during the pressure correction step of the fractional-step
      !! method.
      import :: base_backend_t
      import :: dirps_t
      import :: mesh_t
      implicit none

      class(base_backend_t) :: self
      type(mesh_t), intent(in) :: mesh
      type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
      logical, optional, intent(in) :: lowmem
    end subroutine init_poisson_fft
  end interface

contains

  subroutine base_init(self)
    implicit none

    class(base_backend_t) :: self

  end subroutine base_init

  subroutine get_field_data(self, data, f, dir)
   !! Extract data from field `f` optionally reordering into `dir` orientation.
   !! To output in same orientation as `f`, use `call ...%get_field_data(data, f, f%dir)`
    implicit none

    class(base_backend_t) :: self
    real(dp), dimension(:, :, :), intent(out) :: data !! Output array
    class(field_t), intent(in) :: f !! Field
    integer, optional, intent(in) :: dir !! Desired orientation of output array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(f%dir, direction)

    ! Carry out a reorder if we need, and copy from field to data array
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction)
      call self%reorder(f_temp, f, rdr_dir)
      call self%copy_f_to_data(data, f_temp)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_f_to_data(data, f)
    end if

  end subroutine get_field_data

  subroutine set_field_data(self, f, data, dir)
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f !! Field
    real(dp), dimension(:, :, :), intent(in) :: data !! Input array
    integer, optional, intent(in) :: dir !! Orientation of input array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(direction, f%dir)

    ! Carry out a reorder if we need, and copy from data array to field
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction, f%data_loc)
      call self%copy_data_to_f(f_temp, data)
      call self%reorder(f, f_temp, rdr_dir)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_data_to_f(f, data)
    end if

  end subroutine set_field_data

end module m_base_backend
