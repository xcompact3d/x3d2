!!! test_poisson_cg_eval.f90
!!
!! SPDX-License-Identifier: BSD-3-Clause

program test_poisson_cg_eval
  !! Tests evaluating the Poisson/Laplace operator Lapl(p) on the pressure grid,
  !! used by the iterative Poisson solver.

  use MPI
  
  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t
#ifdef CUDA
#else
  use m_omp_backend
#endif
  use m_common, only: globs_t, DIR_X, DIR_Y, DIR_Z, CELL, POISSON_SOLVER_CG
  use m_poisson_cg, only: laplace_operator_t
  use m_tdsops, only: dirps_t
  
  implicit none

  type(globs_t) :: globs
  class(allocator_t), allocatable :: allocator
  class(base_backend_t), allocatable :: backend
  class(field_t), pointer :: pressure
  class(field_t), pointer :: f
  type(laplace_operator_t) :: lapl
  
  integer :: irank, nproc
  integer :: ierr
  
  logical :: test_pass

  integer, parameter :: nx = 32
  integer, parameter :: ny = 32
  integer, parameter :: nz = 32
  real(dp), parameter :: Lx = 1.0_dp
  real(dp), parameter :: Ly = 1.0_dp
  real(dp), parameter :: Lz = 1.0_dp
  

  test_pass = .true.

  call run_test([nx, ny, nz])

  if (irank == 0) then
    if (.not. test_pass) then
      error stop "Test failed"
    end if
  end if

contains

  subroutine run_test(n)

    integer, dimension(3), intent(in) :: n
    
    type(mesh_t) :: mesh

    call initialise_test(mesh, n)

    call test_constant_field(f, lapl, pressure, mesh)
    call test_variable_field(f, lapl, pressure, mesh)

    ! Finalise test
    call backend%allocator%release_block(pressure)
  
    call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, &
                       MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
                       ierr)
    call MPI_Finalize(ierr)

  end subroutine run_test

  subroutine initialise_test(mesh, n)

    type(mesh_t), intent(out) :: mesh
    integer, dimension(3), intent(in) :: n

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (irank == 0) then
      print *, "Testing the pressure Laplacian operator"
      print *, "Parallel run with", nproc, "ranks"
    end if

    call init_globs(globs, mesh, n, nproc, [Lx, Ly, Lz])
#ifdef CUDA
    if (irank == 0) then
      error stop "CUDA iterative solver not currently supported"
    end if
#else
    allocate(allocator)
    allocator = allocator_t(mesh, SZ)
    allocate(omp_backend_t :: backend)
    backend = omp_backend_t(mesh, allocator)
#endif
    call init_backend(backend)

    pressure => backend%allocator%get_block(DIR_X)
    f => backend%allocator%get_block(DIR_X)
 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (irank == 0) then
      print *, "Initialisation complete"
    end if

  end subroutine initialise_test
  
  subroutine init_globs(globs, mesh, n, nproc, L)
    !! Initialisation for the globs object
    
    type(globs_t), intent(out) :: globs
    type(mesh_t), intent(out) :: mesh
    integer, dimension(3), intent(in) :: n  ! The grid sizes
    integer, intent(in) :: nproc            ! The number of processors
    real(dp), dimension(3), intent(in) :: L ! The domain dimensions

    mesh = mesh_t(n, [1, 1, nproc], L)

    globs%poisson_solver_type = POISSON_SOLVER_CG

  end subroutine init_globs

  subroutine init_backend(backend)

    class(base_backend_t), intent(inout) :: backend

    allocate(backend%xdirps)
    allocate(backend%ydirps)
    allocate(backend%zdirps)
    call init_dirps(backend%xdirps, backend%ydirps, backend%zdirps)

    call backend%alloc_tdsops(backend%xdirps%der2nd, DIR_X, &
                              "second-deriv", "compact6")
    call backend%alloc_tdsops(backend%ydirps%der2nd, DIR_Y, &
                              "second-deriv", "compact6")
    call backend%alloc_tdsops(backend%zdirps%der2nd, DIR_Z, &
                              "second-deriv", "compact6")
    
  end subroutine init_backend
  
  subroutine init_dirps(xdirps, ydirps, zdirps)
    
    type(dirps_t), intent(out) :: xdirps, ydirps, zdirps

    xdirps%dir = DIR_X
    ydirps%dir = DIR_Y
    zdirps%dir = DIR_Z
      
  end subroutine init_dirps

  subroutine test_constant_field(f, lapl, pressure, mesh)

    class(field_t), intent(inout) :: f
    type(laplace_operator_t), intent(in) :: lapl
    class(field_t), intent(in) :: pressure
    type(mesh_t), intent(in) :: mesh

    real(dp), dimension(:, :, :), allocatable :: expect
    
    if (irank == 0) then
      print *, "Testing constant field"
    end if

    ! Set pressure field to some constant
    pressure%data = 42
    allocate(expect, mold = f%data)
    expect = 0  ! Correct answer
    f%data = 17 ! Initialise with wrong answer
    
    call lapl%apply(f, pressure, backend)

    ! Check Laplacian evaluation (expect zero)
    call check_soln(mesh, f, expect)

  end subroutine test_constant_field

  subroutine test_variable_field(f, lapl, pressure, mesh)

    use m_common, only: pi
    use m_ordering, only: get_index_dir
    
    class(field_t), intent(inout) :: f
    type(laplace_operator_t), intent(in) :: lapl
    class(field_t), intent(in) :: pressure
    type(mesh_t), intent(in) :: mesh

    real(dp), dimension(:, :, :), allocatable :: expect

    integer, dimension(3) :: n
    real(dp) :: x, y, z
    integer :: i, j, k
    integer :: ii, jj, kk

    real(dp) :: dx, dy, dz
    real(dp) :: Lx, Ly, Lz
    
    if (irank == 0) then
      print *, "Testing variable field"
    end if

    dx = mesh%geo%d(1); dy = mesh%geo%d(2); dz = mesh%geo%d(3)
    Lx = mesh%geo%L(1); Ly = mesh%geo%L(2); Lz = mesh%geo%L(3)
    n = mesh%get_dims(CELL)

    ! Set pressure field to some variable
    allocate(expect, mold = f%data)
    associate(xdirps => backend%xdirps, &
              ydirps => backend%ydirps, &
              zdirps => backend%zdirps)
      do k = 1, n(3)
        do j = 1, n(2)
          do i = 1, n(1)
            x = (mesh%par%n_offset(1) + (i - 1)) * dx
            y = (mesh%par%n_offset(2) + (j - 1)) * dy
            z = (mesh%par%n_offset(3) + (k - 1)) * dz

            ! Need to get Cartesian -> memory layout mapping
            call get_index_dir(ii, jj, kk, i, j, k, &
                               DIR_X, &
                               SZ, n(1), n(2), n(3))

            pressure%data(ii, jj, kk) = cos(2 * pi * (x / Lx)) + &
                                        cos(2 * pi * (y / Ly)) + &
                                        cos(2 * pi * (z / Lz))
            expect(ii, jj, kk) = -((2 * pi / Lx)**2 * cos(2 * pi * (x / Lx)) + &
                                   (2 * pi / Ly)**2 * cos(2 * pi * (y / Ly)) + &
                                   (2 * pi / Lz)**2 * cos(2 * pi * (z / Lz)))
          end do
        end do
      end do
      f%data = 17 ! Initialise with wrong answer
    end associate

    call lapl%apply(f, pressure, backend)

    ! Check Laplacian evaluation
    call check_soln(mesh, f, expect)
    
  end subroutine test_variable_field

  subroutine check_soln(mesh, soln, expect, opttol)

    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(in) :: soln
    real(dp), dimension(:, :, :), intent(in) :: expect
    real(dp), intent(in), optional :: opttol

    real(dp) :: rms
    real(dp) :: tol

    integer, dimension(3) :: ng
    integer :: n
    
    if (present(opttol)) then
      tol = opttol
    else
      tol = 1.0e-8_dp
    end if

    ng = mesh%get_global_dims(CELL)
    n = product(ng)
    
    rms = sum((soln%data - expect)**2)
    call MPI_Allreduce(MPI_IN_PLACE, rms, 1, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, &
                       ierr)
    rms = sqrt(rms / n)

    if (rms /= rms) then ! NAN check
      print *, "- SEVERE ERROR: RMS=NAN"
      test_pass = .false.
    else
      if (rms > tol) then
        test_pass = .false.

        if (irank == 0) then
          print *, "- FAILED RMS(err) = ", rms
        end if
      else
        if (irank == 0) then
          print *, "- PASSED"
        end if
      end if
    end if

    test_pass = .false.

  end subroutine check_soln

end program test_poisson_cg_eval
