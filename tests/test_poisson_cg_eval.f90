!!! test_poisson_cg_eval.f90
!!
!! SPDX-License-Identifier: BSD-3-Clause

program test_poisson_cg_eval
  !! Tests evaluating the Poisson/Laplace operator Lapl(p) on the pressure grid,
  !! used by the iterative Poisson solver.

  use MPI
  
  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
#ifdef CUDA
#else
  use m_omp_backend
#endif
  use m_common, only: globs_t, DIR_X, DIR_Y, DIR_Z
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
  
  call initialise_test()

  test_pass = .true.
 
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (irank == 0) then
    print *, "Initialisation complete"
  end if
  
  ! Run test
  call test_constant_field(pressure, lapl, f)
  call test_variable_field(pressure, lapl, f)
  
  ! Finalise test
  call backend%allocator%release_block(pressure)
  
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

  subroutine initialise_test()

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    if (irank == 0) then
      print *, "Testing the pressure Laplacian operator"
      print *, "Parallel run with", nproc, "ranks"
    end if

    print *, "INIT GLOBS"
    call init_globs(globs, [nx, ny, nz], nproc, [Lx, Ly, Lz])
    print *, "SELECT BACKEND"
#ifdef CUDA
    if (irank == 0) then
      error stop "CUDA iterative solver not currently supported"
    end if
#else
    allocate(allocator)
    allocator = allocator_t(globs%nx_loc, globs%ny_loc, globs%nz_loc, SZ)
    allocate(omp_backend_t :: backend)
    backend = omp_backend_t(globs, allocator)
#endif
    print *, "INIT BACKEND"
    call init_backend(backend, globs, nproc, irank)

    pressure => backend%allocator%get_block(DIR_X)
    f => backend%allocator%get_block(DIR_X)

  end subroutine initialise_test
  
  subroutine init_globs(globs, n, nproc, L)
    !! Initialisation for the globs object
    
    type(globs_t), intent(out) :: globs
    integer, dimension(3), intent(in) :: n  ! The grid sizes
    integer, intent(in) :: nproc            ! The number of processors
    real(dp), dimension(3), intent(in) :: L ! The domain dimensions

    globs%nx = n(1)
    globs%ny = n(2)
    globs%nz = n(3)

    globs%Lx = L(1)
    globs%Ly = L(2)
    globs%Lz = L(3)

    ! XXX: This requires the grid to be factorisable by nproc
    globs%nx_loc = globs%nx / nproc
    globs%ny_loc = globs%ny / nproc
    globs%nz_loc = globs%nz / nproc

    globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
    globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
    globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

    ! XXX: currently no way to select between periodic/non-periodic BCs
    globs%dx = globs%Lx / globs%nx
    globs%dy = globs%Ly / globs%ny
    globs%dz = globs%Lz / globs%nz

  end subroutine init_globs

  subroutine init_backend(backend, globs, nproc, irank)

    class(base_backend_t), intent(inout) :: backend
    type(globs_t), intent(in) :: globs
    integer, intent(in) :: nproc
    integer, intent(in) :: irank

    print *, "INIT DERPS"
    allocate(backend%xdirps)
    allocate(backend%ydirps)
    allocate(backend%zdirps)
    call init_dirps(backend%xdirps, backend%ydirps, backend%zdirps, globs, nproc, irank)

    print *, "ALLOCATE derivatives"
    call backend%alloc_tdsops(backend%xdirps%der2nd, backend%xdirps%n, &
                              backend%xdirps%d, "second-deriv", "compact6")
    call backend%alloc_tdsops(backend%ydirps%der2nd, backend%ydirps%n, &
                              backend%ydirps%d, "second-deriv", "compact6")
    call backend%alloc_tdsops(backend%zdirps%der2nd, backend%zdirps%n, &
                              backend%zdirps%d, "second-deriv", "compact6")
    print *, "DONE"
    
  end subroutine init_backend
  
  subroutine init_dirps(xdirps, ydirps, zdirps, globs, nproc, irank)

    use m_domain, only: domain_decomposition
    
    type(dirps_t), intent(out) :: xdirps, ydirps, zdirps
    type(globs_t), intent(in) :: globs
    integer, intent(in) :: nproc
    integer, intent(in) :: irank

    xdirps%nproc_dir = globs%nproc_x
    ydirps%nproc_dir = globs%nproc_y
    zdirps%nproc_dir = globs%nproc_z

    xdirps%n = globs%nx_loc
    ydirps%n = globs%ny_loc
    zdirps%n = globs%nz_loc

    xdirps%L = globs%Lx
    ydirps%L = globs%Ly
    zdirps%L = globs%Lz

    xdirps%d = globs%dx
    ydirps%d = globs%dy
    zdirps%d = globs%dz

    xdirps%n_blocks = globs%n_groups_x
    ydirps%n_blocks = globs%n_groups_y
    zdirps%n_blocks = globs%n_groups_z

    xdirps%dir = DIR_X
    ydirps%dir = DIR_Y
    zdirps%dir = DIR_Z

    print *, "DECOMP"
    call domain_decomposition(xdirps, ydirps, zdirps, &
                              irank, nproc)
      
  end subroutine init_dirps

  subroutine test_constant_field(pressure, lapl, f)

    class(field_t), intent(in) :: pressure
    type(laplace_operator_t), intent(in) :: lapl
    class(field_t), intent(inout) :: f

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
    call check_soln(f, expect)

  end subroutine test_constant_field

  subroutine test_variable_field(pressure, lapl, f)

    use m_common, only: pi
    use m_ordering, only: get_index_dir
    
    class(field_t), intent(in) :: pressure
    type(laplace_operator_t), intent(in) :: lapl
    class(field_t), intent(inout) :: f

    real(dp), dimension(:, :, :), allocatable :: expect

    real(dp) :: x, y, z
    integer :: i, j, k
    integer :: ii, jj, kk
    
    if (irank == 0) then
      print *, "Testing variable field"
    end if

    ! Set pressure field to some variable
    allocate(expect, mold = f%data)
    associate(xdirps => backend%xdirps, &
              ydirps => backend%ydirps, &
              zdirps => backend%zdirps, &
              dx => globs%dx, dy => globs%dy, dz => globs%dz, &
              Lx => globs%Lx, Ly => globs%Ly, Lz => globs%Lz)
      do k = 1, zdirps%n
        do j = 1, ydirps%n
          do i = 1, xdirps%n
            x = (xdirps%n_offset + (i - 1)) * dx
            y = (ydirps%n_offset + (j - 1)) * dy
            z = (zdirps%n_offset + (k - 1)) * dz

            ! Need to get Cartesian -> memory layout mapping
            call get_index_dir(ii, jj, kk, i, j, k, &
                               DIR_X, &
                               SZ, xdirps%n, ydirps%n, zdirps%n)

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
    call check_soln(f, expect)
    
  end subroutine test_variable_field

  subroutine check_soln(soln, expect, opttol)

    class(field_t), intent(in) :: soln
    real(dp), dimension(:, :, :), intent(in) :: expect
    real(dp), intent(in), optional :: opttol

    real(dp) :: rms
    real(dp) :: tol

    integer :: n
    
    if (present(opttol)) then
      tol = opttol
    else
      tol = 1.0e-8_dp
    end if

    n = backend%xdirps%n * backend%ydirps%n * backend%zdirps%n
    
    rms = sum((soln%data - expect)**2)
    call MPI_Allreduce(MPI_IN_PLACE, rms, 1, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, &
                       ierr)
    rms = sqrt(rms / n)
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

  end subroutine check_soln

end program test_poisson_cg_eval
