!!! PETSc-based implementation of the iterative Poisson solver

module m_cg_types
  !! Types module providing the context type required by the PETSc matrix-free
  !! operator.

  use m_base_backend, only: base_backend_t
  use m_allocator, only: field_t
  use m_base_poisson_cg, only: laplace_operator_t

  implicit none

  private

  type, public :: mat_ctx_t
    class(base_backend_t), pointer :: backend
    type(laplace_operator_t) :: lapl
    class(field_t), pointer :: xfield
    class(field_t), pointer :: ffield
  end type mat_ctx_t

  interface mat_ctx_t
    module procedure init_ctx
  end interface mat_ctx_t

contains

  function init_ctx(backend, lapl, dir) result(ctx)

    class(base_backend_t), target, intent(in) :: backend
    class(laplace_operator_t), intent(in) :: lapl
    integer, intent(in) :: dir
    type(mat_ctx_t) :: ctx

    ctx%xfield => backend%allocator%get_block(dir)
    ctx%ffield => backend%allocator%get_block(dir)
    ctx%backend => backend
    ctx%lapl = lapl

  end function init_ctx

end module m_cg_types

module m_poisson_cg_backend
  !! Module implementing a Poisson solver based on the (preconditioned)
  !! Conjugate Gradient method using PETSc.

  use m_base_poisson_cg, only: poisson_solver_t, poisson_precon_t, laplace_operator_t

  use petsc

  use m_cg_types

  use m_common, only: dp, DIR_X, DIR_C, CELL
  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: init_solver
  public :: init_precon

  type, extends(poisson_solver_t), public :: poisson_solver_impl
    !! Conjugate Gradient based Poisson solver using PETSc as a backend.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.
    type(mat_ctx_t) :: ctx
    type(tKSP) :: ksp  ! The solver
    type(tMat) :: Amat ! The operator matrix
    type(tVec) :: fvec ! The RHS vector
    type(tVec) :: pvec ! The solution vector
    type(laplace_operator_t) :: lapl
    class(poisson_precon_t), allocatable :: precon
  contains
    procedure, public :: solve => solve_petsc
    procedure :: create_operator
    procedure :: create_vectors
    procedure :: create_solver
  end type poisson_solver_impl

  type, extends(poisson_precon_t), public :: poisson_precon_impl
    !! The PETSc implementation of the Poisson preconditioner, implements a 2nd
    !! order finite difference approximation of the Laplacian.
    type(tMat) :: Pmat ! The preconditioner matrix
  contains
    procedure :: init => init_precon_petsc
    procedure :: apply => petsc_apply_precon
  end type poisson_precon_impl

  interface MatCreateShell
    !! Defines the interface to the external (PETSc) function to create a
    !! matrix-free operator.
    subroutine MatCreateShell(comm, nrow_l, ncol_l, nrow_g, ncol_g, ctx, M, &
                              ierr)
      use petsc
      use m_cg_types
      integer :: comm
      integer :: nrow_l ! Local number of rows
      integer :: ncol_l ! Local number of columns
      integer :: nrow_g ! Global number of rows
      integer :: ncol_g ! Global number of columns
      type(mat_ctx_t) :: ctx ! The shell matrix context
      type(tMat) :: M   ! The matrix object
      integer :: ierr
    end subroutine MatCreateShell
  end interface MatCreateShell

  interface MatShellSetContext
    !! Defines the interface to the external (PETSc) function to store
    !! application-dependent information required by the matrix-free operator.
    subroutine MatShellSetContext(M, ctx, ierr)
      use petsc
      use m_cg_types
      type(tMat) :: M      ! The matrix object
      type(mat_ctx_t) :: ctx ! The shell matrix context
      integer :: ierr
    end subroutine MatShellSetContext
  end interface MatShellSetContext

  interface MatShellGetContext
    !! Defines the interface to the external (PETSc) function to retrieve
    !! application-dependent information from the matrix-free operator.
    subroutine MatShellGetContext(M, ctx, ierr)
      use petsc
      use m_cg_types
      type(tMat) :: M      ! The matrix object
      type(mat_ctx_t) :: ctx ! The shell matrix context
      integer :: ierr
    end subroutine MatShellGetContext
  end interface MatShellGetContext

  interface MatShellSetOperation
    !! Defines the interface to the external (PETSc) function to set the
    !! matrix-free operator procedure that evaluates `f = Mx`.
    subroutine MatShellSetOperation(M, OP, fn, ierr)
      use petsc
      type(tMat) :: M
      integer :: OP
      interface
        subroutine fn(M, x, f, ierr)
          use petsc
          type(tMat) :: M ! The operator
          type(tVec) :: x ! The input vector
          type(tVec) :: f ! The output vector
          integer :: ierr ! The error code
        end subroutine fn
      end interface
      integer :: ierr
    end subroutine MatShellSetOperation
  end interface MatShellSetOperation

  type(mat_ctx_t), save :: ctx_global ! XXX: This sucks!

contains

  function init_precon(backend) result(precon)
    ! Constructs the PETSc preconditioner implementation
    class(base_backend_t), intent(in) :: backend
    class(poisson_precon_t), allocatable :: precon

    allocate (poisson_precon_impl :: precon)
    select type (precon)
    type is (poisson_precon_impl)
      call precon%init(backend)
    class default
      error stop "IMPOSSIBLE"
    end select

  end function init_precon

  subroutine init_precon_petsc(self, backend)
    ! Initialise the PETSc implementation of the preconditioner object

    class(poisson_precon_impl), intent(out) :: self
    class(base_backend_t), intent(in) :: backend

    integer :: n

    type(tMatNullSpace) :: nsp
    integer :: ierr

    integer, parameter :: nnb = 6 ! Number of neighbours (7-point star has 6 neighbours)

    integer, dimension(3) :: dims
    integer :: i, j, k
    real(dp) :: dx, dy, dz

    real(dp), dimension(nnb + 1) :: coeffs
    integer, dimension(nnb + 1) :: cols
    integer :: row

    logical :: initialised

    integer :: istep, jstep, kstep

    ! Ensure PETSc is initialised
    call PetscInitialized(initialised, ierr)
    if (.not. initialised) then
      if (backend%mesh%par%nrank == 0) then
        print *, "Initialising PETSc"
      end if
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    end if
    if (backend%mesh%par%nrank == 0) then
      print *, "PETSc Initialised"
    end if

    ! Create an explicit preconditioner matrix
    n = product(backend%mesh%get_dims(CELL))

    call MatCreate(PETSC_COMM_WORLD, self%Pmat, ierr)
    call MatSetSizes(self%Pmat, n, n, PETSC_DECIDE, PETSC_DECIDE, ierr)
    call MatSetFromOptions(self%Pmat, ierr)
    call MatSeqAIJSetPreallocation(self%Pmat, nnb + 1, PETSC_NULL_INTEGER, &
                                   ierr)
    call MatMPIAIJSetPreallocation(self%Pmat, nnb + 1, PETSC_NULL_INTEGER, &
                                   nnb, PETSC_NULL_INTEGER, &
                                   ierr)
    call MatSetUp(self%Pmat, ierr)

    ! Set the Poisson coefficients
    associate (mesh => backend%mesh)
      call build_index_map(mesh, self%Pmat)

      dims = mesh%get_dims(CELL)
      dx = mesh%geo%d(1); dy = mesh%geo%d(2); dz = mesh%geo%d(3)

      istep = 1
      jstep = dims(1) + 2
      kstep = (dims(1) + 2)*(dims(2) + 2)

      row = kstep + jstep + istep + 1
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            coeffs = 0
            cols = -1 ! Set null (simplifies BCs)
            cols(1) = row

            ! d2pdx2
            coeffs(1) = coeffs(1) - 2/dx**2
            coeffs(2) = 1/dx**2
            coeffs(3) = 1/dx**2
            cols(2) = cols(1) - istep
            cols(3) = cols(1) + istep

            ! d2pdy2
            coeffs(1) = coeffs(1) - 2/dy**2
            coeffs(4) = 1/dy**2
            coeffs(5) = 1/dy**2
            cols(4) = cols(1) - jstep
            cols(5) = cols(1) + jstep

            ! d2pdz2
            coeffs(1) = coeffs(1) - 2/dz**2
            coeffs(6) = 1/dz**2
            coeffs(7) = 1/dz**2
            cols(6) = cols(1) - kstep
            cols(7) = cols(1) + kstep

            ! Push to matrix
            ! Recall Fortran (1-based) -> C (0-based) indexing
            call MatSetValuesLocal(self%Pmat, 1, row - 1, nnb + 1, cols - 1, &
                                   coeffs, INSERT_VALUES, ierr)

            ! Advance row counter
            row = row + 1
          end do
          ! Step in j
          row = row + 2
        end do
        ! Step in k
        row = row + 2*(dims(1) + 2)
      end do
    end associate

    call MatAssemblyBegin(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Pmat, MAT_FINAL_ASSEMBLY, ierr)

    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_VEC, nsp, ierr)
    call MatSetnullSpace(self%Pmat, nsp, ierr)
    call MatNullSpaceDestroy(nsp, ierr)

  end subroutine init_precon_petsc

  module subroutine petsc_apply_precon(self, p, b, backend)
    class(poisson_precon_impl) :: self
    class(field_t), intent(in) :: p
    class(field_t), intent(inout) :: b
    class(base_backend_t), intent(in) :: backend

    type(tVec) :: pVec, bVec
    integer :: ierr

    integer :: n

    n = product(backend%mesh%get_dims(CELL))
    call create_vec(pVec, n)
    call create_vec(bVec, n)

    call copy_field_to_vec(pVec, p, backend)
    call MatMult(self%PMat, pVec, bVec, ierr)
    call copy_vec_to_field(b, bVec, backend)

    call VecDestroy(pVec, ierr)
    call VecDestroy(bVec, ierr)

  end subroutine petsc_apply_precon

  module subroutine solve_petsc(self, p, f, backend)
    class(poisson_solver_impl) :: self
    class(field_t), intent(inout) :: p ! Pressure solution
    class(field_t), intent(in) :: f    ! Poisson RHS
    class(base_backend_t), intent(in) :: backend

    integer :: ierr

    ctx_global = mat_ctx_t(backend, self%lapl, DIR_X)
    call copy_field_to_vec(self%fvec, f, backend)
    call KSPSolve(self%ksp, self%fvec, self%pvec, ierr)
    call copy_vec_to_field(p, self%pvec, backend)

  end subroutine solve_petsc

  module subroutine init_solver(solver, backend, mesh)
    !! Public constructor for the poisson_cg_t type.
    class(poisson_solver_t), allocatable, intent(out) :: solver
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), intent(in) :: mesh

    allocate (poisson_solver_impl :: solver)

    select type (solver)
    type is (poisson_solver_impl)
      solver%precon = init_precon(backend)
      solver%lapl = laplace_operator_t(backend, mesh)
      call init_petsc_cg(solver, backend)
    class default
      ! This should be impossible
      error stop "Failure in allocating PETSc Poisson solver -- this indicates a serious problem"
    end select
  end subroutine init_solver

  subroutine init_petsc_cg(self, backend)
    !! Private constructor for the poisson_cg_t type.
    type(poisson_solver_impl), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend

    integer :: n ! Local problem size

    integer :: ierr

    logical :: initialised

    call PetscInitialized(initialised, ierr)
    if (.not. initialised) then
      if (backend%mesh%par%nrank == 0) then
        print *, "Initialising PETSc"
      end if
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    end if
    if (backend%mesh%par%nrank == 0) then
      print *, "PETSc Initialised"
    end if

    ! Determine local problem size
    n = product(backend%mesh%get_dims(CELL))

    ! Initialise preconditioner and operator matrices
    ! XXX: Add option to use preconditioner as operator (would imply low-order
    !      solution)?
    call self%create_operator(n)

    ! Initialise RHS and solution vectors
    call self%create_vectors(n)

    ! Create the linear system
    call self%create_solver()

  end subroutine init_petsc_cg

  subroutine create_operator(self, n)
    ! Set the PETSc MATVEC to use the x3d2 high-order Laplacian operator
    class(poisson_solver_impl) :: self
    integer, intent(in) :: n ! The local problem size

    type(tMatNullSpace) :: nsp
    integer :: ierr

    call MatCreateShell(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, &
                        PETSC_DETERMINE, self%ctx, self%Amat, ierr)
    call MatShellSetContext(self%Amat, self%ctx, ierr) ! Is this necessary?
    call MatShellSetOperation(self%Amat, MATOP_MULT, poissmult_petsc, ierr)
    call MatSetUp(self%Amat, ierr)

    call MatAssemblyBegin(self%Amat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(self%Amat, MAT_FINAL_ASSEMBLY, ierr)

    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_VEC, nsp, ierr)
    call MatSetnullSpace(self%Amat, nsp, ierr)
    call MatNullSpaceDestroy(nsp, ierr)

  end subroutine create_operator

  subroutine build_index_map(mesh, P)
    ! Builds the map from local indices to the global (equation ordering) index for the preconditioner
#include "petsc/finclude/petscis.h"
    type(mesh_t), intent(in) :: mesh ! The mesh for the simulation
    type(tMat) :: P                  ! The preconditioner matrix
    ISLocalToGlobalMapping map
    integer :: ierr

    integer, dimension(3) :: dims
    integer :: n
    integer :: global_start

    integer, dimension(:), allocatable :: idx

    dims = mesh%get_dims(CELL)
    n = product(dims + 2) ! Size of domain + 1 deep halo

    allocate (idx(n))
    idx(:) = 0

    ! Determine global start point based on PETSc decomposition
    ! | 0 ... P0 ... P0_n-1 | P0_n ... P1 ... P0_n+P1_n-1 | P0_n+P1_n ... P2
    ! i.e. an exclusive scan sum of each rank's allocation
    call MPI_Exscan(product(dims), global_start, 1, MPI_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr)
    if (mesh%par%nrank == 0) then
      global_start = 0
    end if
    global_start = global_start + 1 ! C->F

    ! Build the local->global index map
    call build_interior_index_map(idx, mesh, global_start)
    call build_neighbour_index_map(idx, mesh, global_start)
    idx = idx - 1 ! F->C

    call ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, n, idx, PETSC_COPY_VALUES, map, ierr)
    call MatSetLocalToGlobalMapping(P, map, map, ierr)

    deallocate (idx)

  end subroutine build_index_map

  subroutine build_interior_index_map(idx, mesh, global_start)
    ! Builds the map from local indices to the global (equation ordering) index
    ! for internal (local) indices

    ! The local->global index map, local indices are offset into the domain by 1 in each axis to simplify finite differences
    integer, dimension(:), intent(inout) :: idx
    type(mesh_t), intent(in) :: mesh    ! The mesh for the simulation
    integer, intent(in) :: global_start ! The global offset for this rank in the global equation system

    integer, dimension(3) :: dims
    integer :: nx, ny, nz

    integer :: i, j, k
    integer :: ctr, local_ctr

    dims = mesh%get_dims(CELL)
    nx = dims(1); ny = dims(2); nz = dims(3)

    ctr = 0
    local_ctr = ((nx + 2)*(ny + 2)) + (nx + 2) + 2
    do k = 2, nz + 1
      do j = 2, ny + 1
        do i = 2, nx + 1
          idx(local_ctr) = global_start + ctr
          local_ctr = local_ctr + 1
          ctr = ctr + 1
        end do
        local_ctr = local_ctr + 2
      end do
      local_ctr = local_ctr + 2*(nx + 2)
    end do

  end subroutine build_interior_index_map

  subroutine build_neighbour_index_map(idx, mesh, global_start)
    ! Builds the map from local indices to the global (equation ordering) index
    ! for partition boundary indices

    use mpi

    integer, dimension(:), intent(inout) :: idx
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: global_start

    integer, dimension(3) :: dims
    integer :: nx, ny, nz
    integer, dimension(4) :: myinfo
    integer, dimension(4, 2, 3) :: info
    integer :: d

    integer, dimension(:, :, :), allocatable :: halobuf_x, halobuf_y, halobuf_z
    integer :: ctr
    integer :: i, j, k

    integer, dimension(12) :: requests
    integer :: tag, nbtag, nproc
    integer :: ierr

    dims = mesh%get_dims(CELL)
    nx = dims(1); ny = dims(2); nz = dims(3)

    ! Create and fill halobuffers
    allocate (halobuf_x(2, ny, nz))
    allocate (halobuf_y(nx, 2, nz))
    allocate (halobuf_z(nx, ny, 2))
    halobuf_x = -1; halobuf_y = -1; halobuf_z = -1

    myinfo = [global_start, dims(1), dims(2), dims(3)]
    nproc = mesh%par%nproc
    tag = mesh%par%nrank*max(nproc, 6)
    do d = 1, 3
      !! Recv left and right
      !  Right Recv
      nbtag = mesh%par%pnext(d)*max(nproc, 6)
      call MPI_IRecv(info(:, 2, d), 4, MPI_INTEGER, mesh%par%pnext(d), nbtag + 2 * (d - 1) + 1, & 
                     MPI_COMM_WORLD, requests(4*(d - 1) + 1), ierr)
      !  Left Recv
      nbtag = mesh%par%pprev(d)*max(nproc, 6)
      call MPI_IRecv(info(:, 1, d), 4, MPI_INTEGER, mesh%par%pprev(d), nbtag + 2 * (d - 1) + 2, & 
                     MPI_COMM_WORLD, requests(4*(d - 1) + 2), ierr)

      !! Send left and right
      !  Left Send
      call MPI_ISend(myinfo, 4, MPI_INTEGER, mesh%par%pprev(d), tag + 2 * (d - 1) + 1, &
                     MPI_COMM_WORLD, requests(4*(d - 1) + 3), ierr)
      !  Right Send
      call MPI_ISend(myinfo, 4, MPI_INTEGER, mesh%par%pnext(d), tag + 2 * (d - 1) + 2, &
                     MPI_COMM_WORLD, requests(4*(d - 1) + 4), ierr)
    end do
    call MPI_Waitall(12, requests, MPI_STATUS_IGNORE, ierr)

    !! X halos
    ! Left halo
    associate (offset_left => info(1, 1, 1), &
               nx_left => info(2, 1, 1), &
               ny_left => info(3, 1, 1), &
               nz_left => info(4, 1, 1))
      ctr = offset_left + (nx_left - 1) ! Global starting index -> xend
      do k = 1, nz_left
        do j = 1, ny_left
          halobuf_x(1, j, k) = ctr
          ctr = ctr + nx_left ! Step in j
        end do
      end do
    end associate
    ! Right halo
    associate (offset_right => info(1, 2, 1), &
               nx_right => info(2, 2, 1), &
               ny_right => info(3, 2, 1), &
               nz_right => info(4, 2, 1))
      ctr = offset_right ! Global starting index == xstart
      do k = 1, nz_right
        do j = 1, ny_right
          halobuf_x(2, j, k) = ctr
          ctr = ctr + nx_right ! Step in j
        end do
      end do
    end associate

    !! Y halos
    ! Bottom halo
    associate (offset_down => info(1, 1, 2), &
               nx_down => info(2, 1, 2), &
               ny_down => info(3, 1, 2), &
               nz_down => info(4, 1, 2))
      ctr = offset_down + (ny_down - 1)*nx_down ! Global starting index -> yend
      do k = 1, nz_down
        do i = 1, nx_down
          halobuf_y(i, 1, k) = ctr
          ctr = ctr + 1 ! Step in i
        end do
        ctr = ctr - nx_down  ! Reset counter to start of line
        ctr = ctr + (nx_down*ny_down) ! Step in k
      end do
    end associate
    ! Top halo
    associate (offset_up => info(1, 2, 2), &
               nx_up => info(2, 2, 2), &
               ny_up => info(3, 2, 2), &
               nz_up => info(4, 2, 2))
      ctr = offset_up ! Global starting index == ystart
      do k = 1, nz_up
        do i = 1, nx_up
          halobuf_y(i, 2, k) = ctr
          ctr = ctr + 1 ! Step in i
        end do
        ctr = ctr - nx_up  ! Reset counter to start of line
        ctr = ctr + (nx_up*ny_up) ! Step in k
      end do
    end associate

    !! Z halos
    ! Back halo
    associate (offset_back => info(1, 1, 3), &
               nx_back => info(2, 1, 3), &
               ny_back => info(3, 1, 3), &
               nz_back => info(4, 1, 3))
      ctr = offset_back + (nz_back - 1)*ny_back*nx_back ! Global starting index -> zend
      do j = 1, ny_back
        do i = 1, nx_back
          halobuf_z(i, j, 1) = ctr
          ctr = ctr + 1 ! Step in i
        end do
      end do
    end associate
    ! Front halo
    associate (offset_front => info(1, 2, 3), &
               nx_front => info(2, 2, 3), &
               ny_front => info(3, 2, 3), &
               nz_front => info(4, 2, 3))
      ctr = offset_front ! Global startin index == zstart
      do j = 1, ny_front
        do i = 1, nx_front
          halobuf_z(i, j, 2) = ctr
          ctr = ctr + 1 ! Step in i
        end do
      end do
    end associate

    !! Map my neighbours indices into my halos
    ctr = 1
    do k = 1, nz + 2
      do j = 1, ny + 2

        ! Left halo
        if ((j > 1) .and. (j < (ny + 2)) .and. &
            (k > 1) .and. (k < (nz + 2))) then
          idx(ctr) = halobuf_x(1, j - 1, k - 1)
        end if
        ctr = ctr + 1

        do i = 2, nx + 1

          if ((k > 1) .and. (k < (nz + 2))) then
            if (j == 1) then
              ! Bottom halo
              idx(ctr) = halobuf_y(i - 1, 1, k - 1)
            else if (j == (ny + 2)) then
              ! Top halo
              idx(ctr) = halobuf_y(i - 1, 2, k - 1)
            end if
          else if ((j > 1) .and. (j < (ny + 2))) then
            if (k == 1) then
              ! Back halo
              idx(ctr) = halobuf_z(i - 1, j - 1, 1)
            else if (k == (nz + 2)) then
              ! Front halo
              idx(ctr) = halobuf_z(i - 1, j - 1, 2)
            end if
          end if

          ctr = ctr + 1
        end do

        ! Right halo
        if ((j > 1) .and. (j < (ny + 2)) .and. &
            (k > 1) .and. (k < (nz + 2))) then
          idx(ctr) = halobuf_x(2, j - 1, k - 1)
        end if
        ctr = ctr + 1

      end do
    end do

  end subroutine build_neighbour_index_map

  subroutine create_vectors(self, n)
    ! Allocates the pressure and forcing vectors.

    class(poisson_solver_impl) :: self ! The Poisson solver
    integer, intent(in) :: n          ! The local vector size

    call create_vec(self%fvec, n)
    call create_vec(self%pvec, n)

  end subroutine create_vectors

  subroutine create_vec(v, n)
    ! Utility subroutine to allocate a PETSc vector.

    type(tVec), intent(out) :: v ! The vector
    integer, intent(in) :: n     ! The local vector size

    integer :: ierr

    call VecCreate(PETSC_COMM_WORLD, v, ierr)
    call VecSetSizes(v, n, PETSC_DETERMINE, ierr)
    call VecSetFromOptions(v, ierr)

    call VecAssemblyBegin(v, ierr)
    call VecAssemblyEnd(v, ierr)

  end subroutine create_vec

  subroutine create_solver(self)
    ! Sets up the PETSc linear solver.

    class(poisson_solver_impl) :: self

    integer :: ierr

    associate (precon => self%precon)
      select type (precon)
      type is (poisson_precon_impl)
        call KSPCreate(PETSC_COMM_WORLD, self%ksp, ierr)
        call KSPSetOperators(self%ksp, self%Amat, precon%Pmat, ierr)
        call KSPSetFromOptions(self%ksp, ierr)
        call KSPSetInitialGuessNonzero(self%ksp, PETSC_TRUE, ierr)
      class default
        error stop "Poisson preconditioner type is wrong"
      end select
    end associate

  end subroutine create_solver

  subroutine poissmult_petsc(M, x, f, ierr)
    !! Computes the action of the Poisson operator, i.e. `f = Mx` where `M` is
    !! the discrete Laplacian.
    type(tMat) :: M ! The operator
    type(tVec) :: x ! The input vector
    type(tVec) :: f ! The output vector
    integer :: ierr ! The error code

    type(mat_ctx_t) :: ctx

    ierr = 0; 
    ! XXX: Fixme
    ! call MatShellGetContext(M, ctx, ierr)
    associate (matrix => M); end associate ! Silence unused argument
    ! print *, ctx%foo
    ctx = ctx_global

    call copy_vec_to_field(ctx%xfield, x, ctx%backend)
    call ctx%lapl%apply(ctx%ffield, ctx%xfield, ctx%backend)
    call copy_field_to_vec(f, ctx%ffield, ctx%backend)

  end subroutine poissmult_petsc

  subroutine copy_vec_to_field(f, v, backend)
    !! Copies the contents of a PETSc vector into an x3d2 field object
    ! XXX: This can be avoided if a field can wrap the vector memory
    class(field_t), intent(inout) :: f ! The destination field
    type(tVec) :: v                    ! The source vector
    class(base_backend_t), intent(in) :: backend

    real(dp), dimension(:), pointer :: vdata
    real(dp), dimension(:, :, :), pointer :: vdata3d
    integer :: ierr

    integer, dimension(3) :: dims
    integer :: nx, ny, nz

    dims = backend%mesh%get_dims(CELL)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    ! Local copy
    call VecGetArrayReadF90(v, vdata, ierr)
    if (nx*ny*nz /= size(vdata)) then
      print *, "Vector and field sizes are incompatible (padding?)"
      stop 1
    end if
    vdata3d(1:nx, 1:ny, 1:nz) => vdata(:) ! Get a 3D representation of the vector
    call backend%set_field_data(f, vdata3d, DIR_C)
    call VecRestoreArrayReadF90(v, vdata, ierr)
    nullify (vdata3d)

    ! Halo exchange

  end subroutine copy_vec_to_field

  subroutine copy_field_to_vec(v, f, backend)
    !! Copies the contents of an x3d2 field object into a PETSc vector
    ! XXX: This can be avoided if a field can wrap the vector memory
    type(tVec) :: v                 ! The destination vector.
    class(field_t), intent(in) :: f ! The source field.
    class(base_backend_t), intent(in) :: backend

    real(dp), dimension(:), pointer :: vdata
    real(dp), dimension(:, :, :), pointer :: vdata3d
    integer :: ierr

    integer, dimension(3) :: dims
    integer :: nx, ny, nz

    dims = backend%mesh%get_dims(CELL)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    ! Local copy
    call VecGetArrayF90(v, vdata, ierr)
    if (nx*ny*nz /= size(vdata)) then
      print *, "Vector and field sizes are incompatible (padding?)"
      stop 1
    end if
    vdata3d(1:nx, 1:ny, 1:nz) => vdata(:) ! Get a 3D representation of the vector
    call backend%get_field_data(vdata3d, f, DIR_C)
    call VecRestoreArrayF90(v, vdata, ierr)
    nullify (vdata3d)

    ! Halo exchange
    call VecAssemblyBegin(v, ierr)
    call VecAssemblyEnd(v, ierr)

  end subroutine copy_field_to_vec

end module m_poisson_cg_backend
