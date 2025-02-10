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
  use m_common, only: DIR_X, DIR_Y, DIR_Z, CELL, POISSON_SOLVER_CG
  use m_poisson_cg, only: laplace_operator_t, poisson_precon_t

  implicit none

  class(allocator_t), allocatable :: allocator
  class(base_backend_t), allocatable :: backend
  class(field_t), pointer :: pressure
  class(field_t), pointer :: f

  integer :: irank, nproc
  integer :: ierr

  logical :: test_pass

  real(dp), parameter :: Lx = 1.0_dp
  real(dp), parameter :: Ly = 1.0_dp
  real(dp), parameter :: Lz = 1.0_dp
  integer, parameter :: nref = 4

  real(dp), dimension(nref) :: e

  type(laplace_operator_t) :: lapl
#ifdef HAVE_CG
  type(poisson_precon_t) :: precon
#endif

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (irank == 0) then
    print *, "Testing the pressure Laplacian operator"
    print *, "Parallel run with", nproc, "ranks"
  end if

  test_pass = .true.

  call test_driver(lapl)
#ifdef HAVE_CG
  call test_driver(precon)
#endif

  call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, &
                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
                     ierr)
  call MPI_Finalize(ierr)

  if (irank == 0) then
    if (.not. test_pass) then
      error stop "Test failed"
    end if
  end if

contains

  subroutine test_driver(linear_operator)

    class(*), intent(in) :: linear_operator

    integer :: nx
    integer :: ny
    integer :: nz
    integer :: i

    integer :: order
    character(len=:), allocatable :: opname

    select type (linear_operator)
    type is (laplace_operator_t)
      order = 6
      opname = "High-Order Laplacian"
    type is (poisson_precon_t)
      order = 2
      opname = "Low-Order preconditioner"
    class default
      error stop "Unknown linear operator"
    end select

    nx = 16; ny = 16; nz = 16
    ! nx = nx * nproc
    ny = ny*nproc
    ! nz = nz * nproc
    do i = 1, nref
      if (irank == 0) then
        print *, "---------------------------------"
        print *, "Testing refinement level ", i - 1
        print *, "Using the ", opname
      end if
      e(i) = run_test([nx, ny, nz], linear_operator)
      nx = 2*nx; ny = 2*ny; nz = 2*nz
    end do

    do i = 2, nref
      if (e(i) > 2.2*(e(i - 1)/(2**order))) then
        if (irank == 0) then
          print *, "Error convergence ", i, " failed ", &
            e(i), e(i - 1)/(2**order)
        end if
        test_pass = .false.
      else
        if (irank == 0) then
          print *, "Error convergence ", i, "satisfied ", &
            e(i), e(i - 1)/(2**order)
        end if
      end if
    end do

  end subroutine test_driver

  real(dp) function run_test(n, linear_operator)

    integer, dimension(3), intent(in) :: n
    class(*), intent(in) :: linear_operator

    type(mesh_t) :: mesh
    real(dp) :: rms_err

    call initialise_test(mesh, n)

    call test_constant_field(f, linear_operator, pressure, mesh)
    rms_err = test_variable_field(f, linear_operator, pressure, mesh)

    ! Finalise test
    call backend%allocator%release_block(pressure)

    run_test = rms_err

  end function run_test

  subroutine initialise_test(mesh, n)

    type(mesh_t), intent(out) :: mesh
    integer, dimension(3), intent(in) :: n

    call init_globs(mesh, n, nproc, [Lx, Ly, Lz])
#ifdef CUDA
    if (irank == 0) then
      error stop "CUDA iterative solver not currently supported"
    end if
#else
    if (allocated(allocator)) then
      deallocate (allocator)
    end if
    allocate (allocator)
    allocator = allocator_t(mesh, SZ)

    if (allocated(backend)) then
      deallocate (backend)
    end if
    allocate (omp_backend_t :: backend)
    backend = omp_backend_t(mesh, allocator)
#endif
    lapl = laplace_operator_t(backend, mesh)
#ifdef HAVE_CG
    precon = poisson_precon_t(backend)
#endif

    ! Main solver calls Poisson in the DIR_Z orientation
    pressure => backend%allocator%get_block(DIR_Z)
    f => backend%allocator%get_block(DIR_Z)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (irank == 0) then
      print *, "Initialisation complete"
    end if

  end subroutine initialise_test

  subroutine init_globs(mesh, n, nproc, L)
    !! Initialisation for the globs object

    type(mesh_t), intent(out) :: mesh
    integer, dimension(3), intent(in) :: n  ! The grid sizes
    integer, intent(in) :: nproc            ! The number of processors
    real(dp), dimension(3), intent(in) :: L ! The domain dimensions

    mesh = mesh_t(n, [1, nproc, 1], L, &
                  ["periodic", "periodic"], &
                  ["periodic", "periodic"], &
                  ["periodic", "periodic"])

  end subroutine init_globs

  subroutine test_constant_field(f, linear_operator, pressure, mesh)

    class(field_t), intent(inout) :: f
    class(*), intent(in) :: linear_operator
    class(field_t), intent(in) :: pressure
    type(mesh_t), intent(in) :: mesh

    real(dp), dimension(:, :, :), allocatable :: expect

    real(dp) :: rms_err

    logical :: check_pass

    if (irank == 0) then
      print *, "Testing constant field"
    end if

    ! Set pressure field to some constant
    pressure%data = 0 !42
    allocate (expect, mold=f%data)
    expect = 0  ! Correct answer
    f%data = 17 ! Initialise with wrong answer

    select type (linear_operator)
    type is (laplace_operator_t)
      call linear_operator%apply(f, pressure, backend)
    type is (poisson_precon_t)
      call linear_operator%apply(pressure, f, backend)
    class default
      error stop "Unsupported linear operator type"
    end select

    ! Check Laplacian evaluation (expect zero)
    rms_err = check_soln(check_pass, mesh, f, expect)

    if (.not. check_pass) then
      test_pass = .false.

      if (irank == 0) then
        print *, "- FAILED RMS(err) = ", rms_err
      end if
    else
      if (irank == 0) then
        print *, "- PASS"
      end if
    end if

  end subroutine test_constant_field

  real(dp) function test_variable_field(f, linear_operator, pressure, mesh)

    use m_common, only: pi
    use m_ordering, only: get_index_dir

    class(field_t), intent(inout) :: f
    class(*), intent(in) :: linear_operator
    class(field_t), intent(in) :: pressure
    type(mesh_t), intent(in) :: mesh

    real(dp), dimension(:, :, :), allocatable :: expect

    integer, dimension(3) :: n
    real(dp) :: x, y, z
    real(dp), dimension(3) :: coords
    integer :: i, j, k
    integer :: ii, jj, kk

    real(dp) :: dx, dy, dz
    real(dp) :: Lx, Ly, Lz

    logical :: check_pass

    if (irank == 0) then
      print *, "Testing variable field"
    end if

    dx = mesh%geo%d(1); dy = mesh%geo%d(2); dz = mesh%geo%d(3)
    Lx = mesh%geo%L(1); Ly = mesh%geo%L(2); Lz = mesh%geo%L(3)
    n = mesh%get_dims(CELL)

    ! Set pressure field to some variable
    allocate (expect, mold=f%data)
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          coords = mesh%get_coordinates(i, j, k)
          x = coords(1); y = coords(2); z = coords(3)

          ! Need to get Cartesian -> memory layout mapping
          call get_index_dir(ii, jj, kk, i, j, k, &
                             DIR_Z, &
                             SZ, n(1), n(2), n(3))

          pressure%data(ii, jj, kk) = cos(2*pi*(x/Lx)) + &
                                      cos(2*pi*(y/Ly)) + &
                                      cos(2*pi*(z/Lz))
          expect(ii, jj, kk) = -((2*pi/Lx)**2*cos(2*pi*(x/Lx)) + &
                                 (2*pi/Ly)**2*cos(2*pi*(y/Ly)) + &
                                 (2*pi/Lz)**2*cos(2*pi*(z/Lz)))
          ! pressure%data(ii, jj, kk) = cos(2 * pi * (z / Lz))
          ! expect(ii, jj, kk) = -((2 * pi / Lz)**2 * cos(2 * pi * (z / Lz)))
        end do
      end do
    end do
    f%data = 0 ! Initialise with wrong answer

    select type (linear_operator)
    type is (laplace_operator_t)
      call linear_operator%apply(f, pressure, backend)
    type is (poisson_precon_t)
      call linear_operator%apply(pressure, f, backend)
    class default
      error stop "Unsupported linear operator type"
    end select

    ! Check Laplacian evaluation
    ! XXX: Note had to relax the tolerance, otherwise obtains RMS(err)~=8e-7
    test_variable_field = check_soln(check_pass, mesh, f, expect, &
                                     opttol=1.0e-6_dp)

  end function test_variable_field

  real(dp) function check_soln(check_pass, mesh, soln, expect, opttol)

    logical, intent(out) :: check_pass
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
    rms = sqrt(rms/n)

    if (rms /= rms) then ! NAN check
      print *, "- SEVERE ERROR: RMS=NAN"
      test_pass = .false.
      check_pass = .false.
    else
      if (rms > tol) then
        check_pass = .false.
      else
        check_pass = .true.
      end if
    end if

    check_soln = rms

  end function check_soln

end program test_poisson_cg_eval
