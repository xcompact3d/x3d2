!!! test_poisson_cg_solve.f90
!!
!! SPDX-License-Identifier: BSD-3-Clause

program test_poisson_cg_solve
  !! Tests solving the Poisson problem using CG.

  use mpi
  use petsc

  use m_common, only: dp, DIR_Z, CELL, pi
  use m_mesh, only: mesh_t
  use m_base_backend, only: base_backend_t
#ifdef CUDA
#else
  use m_omp_backend
#endif
  use m_poisson_cg, only: poisson_cg_t, poisson_cg_t
  use m_ordering, only: get_index_dir

  implicit none

  integer :: irank, nproc

  integer, parameter :: nref = 3 ! Number of refinements to perform
  real(dp), parameter :: Lx = 1.0_dp
  real(dp), parameter :: Ly = 1.0_dp
  real(dp), parameter :: Lz = 1.0_dp

  logical :: test_pass

  call setup()
  call run()
  call teardown()

contains

  subroutine setup()
    !! Perform test program setup

    integer :: ierr

    call MPI_Init(ierr)
    call PETScInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    test_pass = .true.

  end subroutine setup

  subroutine run()

    integer :: nx, ny, nz
    integer :: i

    type(mesh_t), allocatable :: mesh
    class(allocator_t), allocatable :: allocator
    class(base_backend_t), allocatable :: backend
    class(poisson_cg_t), allocatable :: poisson_cg
    class(field_t), pointer :: p, f

    real(dp), dimension(nref) :: rms

    ! Test across multiple refinement levels
    nx = 16; ny = 16; nz = 16
    do i = 1, nref
      mesh = mesh_t([nx, ny, nz], [1, 1, nproc], [Lx, Ly, Lz], &
                    ["periodic", "periodic"], &
                    ["periodic", "periodic"], &
                    ["periodic", "periodic"])
#ifdef CUDA
      error stop "CUDA iterative solver not currently supported"
#else
      allocator = allocator_t(mesh, SZ)
      if (allocated(backend)) then
        deallocate (backend)
      end if
      allocate (omp_backend_t :: backend)
      backend = omp_backend_t(mesh, allocator)
#endif
      poisson_cg = poisson_cg_t(backend)

      ! Main solver calls Poisson in the DIR_Z orientation
      p => backend%allocator%get_block(DIR_Z, CELL)
      f => backend%allocator%get_block(DIR_Z, CELL)

      call set_rhs(f, backend%mesh)

      call poisson_cg%solve(p, f, backend)

      rms(i) = calc_err(p, backend%mesh)

      call backend%allocator%release_block(p)
      call backend%allocator%release_block(f)

      nx = 2*nx
      ny = 2*ny
      nz = 2*nz
    end do

    call test_convergence(rms)

  end subroutine run

  subroutine set_rhs(f, mesh)
    class(field_t), intent(inout) :: f
    type(mesh_t), intent(in) :: mesh

    integer, dimension(3) :: n
    real(dp), dimension(3) :: d, L

    integer :: i, j, k
    integer :: ii, jj, kk
    real(dp) :: x, y, z

    n = mesh%get_dims(CELL)
    d = mesh%geo%d
    L = mesh%geo%L

    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          x = (mesh%par%n_offset(1) + (i - 1))*d(1)
          y = (mesh%par%n_offset(2) + (j - 1))*d(2)
          z = (mesh%par%n_offset(3) + (k - 1))*d(3)

          ! Need to get Cartesian -> memory layout mapping
          call get_index_dir(ii, jj, kk, i, j, k, &
                             DIR_Z, &
                             SZ, n(1), n(2), n(3))

          f%data(ii, jj, kk) = -((2*pi/L(1))**2*cos(2*pi*(x/L(1))) + &
                                 (2*pi/L(2))**2*cos(2*pi*(y/L(2))) + &
                                 (2*pi/L(3))**2*cos(2*pi*(z/L(3))))
        end do
      end do
    end do

  end subroutine set_rhs

  real(dp) function calc_err(p, mesh)
    class(field_t), intent(in) :: p
    type(mesh_t), intent(in) :: mesh

    real(dp) :: err_rms
    real(dp) :: p_an

    integer, dimension(3) :: n
    real(dp), dimension(3) :: d, L

    integer :: i, j, k
    integer :: ii, jj, kk
    real(dp) :: x, y, z

    integer :: ierr

    err_rms = 0.0_dp

    n = mesh%get_dims(CELL)
    d = mesh%geo%d
    L = mesh%geo%L

    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          x = (mesh%par%n_offset(1) + (i - 1))*d(1)
          y = (mesh%par%n_offset(2) + (j - 1))*d(2)
          z = (mesh%par%n_offset(3) + (k - 1))*d(3)

          ! Need to get Cartesian -> memory layout mapping
          call get_index_dir(ii, jj, kk, i, j, k, &
                             DIR_Z, &
                             SZ, n(1), n(2), n(3))

          p_an = cos(2*pi*(x/L(1))) + &
                 cos(2*pi*(y/L(2))) + &
                 cos(2*pi*(z/L(3)))

          err_rms = err_rms + (p%data(ii, jj, kk) - p_an)**2
        end do
      end do
    end do

    call MPI_Allreduce(MPI_IN_PLACE, err_rms, 1, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, &
                       ierr)
    n = mesh%get_global_dims(CELL)
    calc_err = sqrt(err_rms/product(n))

  end function calc_err

  subroutine test_convergence(rms_err)
    real(dp), dimension(:), intent(in) :: rms_err

    integer :: i
    real(dp) :: r
    real(dp) :: tol

    do i = 1, size(rms_err)
      if (irank == 0) then
        print *, "RMS = ", rms_err(i)
      end if
    end do
    do i = 2, size(rms_err)
      r = rms_err(i)/rms_err(i - 1)
      if (irank == 0) then
        print *, "Convergence ratio ", i, ": ", r
      end if

      tol = 1.1_dp/(2**6)
      if (r > tol) then
        if (irank == 0) then
          print *, "- Error convergence failed, tolerance = ", tol
        end if

        test_pass = .false.
      end if
    end do

  end subroutine test_convergence

  subroutine teardown()
    !! Perform test program cleanup
    use petsc

    integer :: ierr

    ! Reduce test status
    call MPI_Allreduce(MPI_IN_PLACE, test_pass, 1, &
                       MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
                       ierr)

    ! Cleanup parallelism
    call PetscFinalize(ierr)
    call MPI_Finalize(ierr)

    ! Report pass/fail
    if (irank == 0) then
      if (test_pass) then
        print *, "PASS"
      else
        print *, "FAIL"
        stop 1
      end if
    end if

  end subroutine teardown

end program test_poisson_cg_solve
