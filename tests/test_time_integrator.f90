program test_omp_adamsbashforth
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_common, only: dp, DIR_X, pi
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_time_integrator, only: time_intg_t
#ifdef CUDA
  use cudafor

  use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
  use m_cuda_backend, only: cuda_backend_t
#else
  use m_omp_backend, only: omp_backend_t
#endif

  implicit none

  logical :: allpass = .true.
  integer :: i, j, k, istartup, nrank, nproc, ierr
  integer :: nstep0 = 64, nstep, nrun = 4, norder = 4
  real(dp), allocatable, dimension(:) :: err
  real(dp), allocatable, dimension(:) :: norm
  real(dp) :: dt0 = 0.01_dp, dt, order
  real(dp) :: u0
  class(field_t), pointer :: u, v, w
  class(field_t), pointer :: du, dv, dw
  integer, dimension(3) :: dims_global, nproc_dir
  real(dp), dimension(3) :: L_global

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(mesh_t), allocatable :: mesh
#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
#endif
  class(time_intg_t), allocatable :: time_integrator

  ! initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
#endif

  ! Global number of cells in each direction
  dims_global = [1, 1, 1]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, 1]

  mesh = mesh_t(dims_global, nproc_dir, L_global)

  ! allocate object
#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh, 1)
  allocator => cuda_allocator
  if (nrank == 0) print *, 'CUDA allocator instantiated'

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  if (nrank == 0) print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(mesh, 1)
  allocator => omp_allocator
  if (nrank == 0) print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  if (nrank == 0) print *, 'OpenMP backend instantiated'
#endif

  time_integrator = time_intg_t(allocator=allocator, &
                                backend=backend, order=norder)
  print *, 'time integrator instantiated'

  ! allocate memory
  u => allocator%get_block(DIR_X)
  v => allocator%get_block(DIR_X)
  w => allocator%get_block(DIR_X)

  du => allocator%get_block(DIR_X)
  dv => allocator%get_block(DIR_X)
  dw => allocator%get_block(DIR_X)

  allocate (norm(nrun))

  ! compute l2 norm for various step sizes
  do k = 1, norder
    time_integrator%order = k
    dt = dt0
    nstep = nstep0
    do j = 1, nrun
      ! initial condition
      time_integrator%istep = 1

      ! compute l2 norm for a given step size
      allocate (err(nstep))
#ifdef CUDA
      select type (u)
      type is (cuda_field_t)
        u%data_d(1, 1, 1) = 1.0_dp
      end select
#else
      u%data(1, 1, 1) = 1.0_dp
#endif
      ! startup
      istartup = k - 1
      do i = 1, istartup
#ifdef CUDA
        select type (u)
        type is (cuda_field_t)
          u0 = u%data_d(1, 1, 1)
        end select
        select type (du)
        type is (cuda_field_t)
          du%data_d(1, 1, 1) = dahlquist_rhs(u0)
        end select
#else
        du%data(1, 1, 1) = dahlquist_rhs(u%data(1, 1, 1))
#endif
        call time_integrator%step(u, v, w, du, dv, dw, dt)
#ifdef CUDA
        select type (u)
        type is (cuda_field_t)
          u%data_d(1, 1, 1) = dahlquist_exact_sol(real(i, dp)*dt)
        end select
#else
        u%data(1, 1, 1) = dahlquist_exact_sol(real(i, dp)*dt)
#endif
      end do

      ! post-startup
      do i = 1, nstep
#ifdef CUDA
        select type (u)
        type is (cuda_field_t)
          u0 = u%data_d(1, 1, 1)
        end select
        select type (du)
        type is (cuda_field_t)
          du%data_d(1, 1, 1) = dahlquist_rhs(u0)
        end select
#else
        du%data(1, 1, 1) = dahlquist_rhs(u%data(1, 1, 1))
#endif
        call time_integrator%step(u, v, w, du, dv, dw, dt)
#ifdef CUDA
        select type (u)
        type is (cuda_field_t)
          u0 = u%data_d(1, 1, 1)
        end select
        err(i) = u0 - dahlquist_exact_sol( &
                 real(i + istartup, dp)*dt)
#else
        err(i) = u%data(1, 1, 1) - dahlquist_exact_sol( &
                 real(i + istartup, dp)*dt)
#endif
      end do

      ! compute l2 norms
      norm(j) = norm2(err)
      norm(j) = sqrt(norm(j)*norm(j)/real(nstep, dp))
      print *, err(nstep)
      deallocate (err)

      ! refine time stepping
      dt = dt/2.0_dp
      nstep = nstep*2
    end do

    ! check order of convergence
    order = log(norm(nrun - 1)/norm(nrun))/log(2.0_dp)
    print *, 'order', order
    if (abs(order - real(k, dp)) > 0.1_dp) then
      allpass = .false.
      write (stderr, '(a)') 'Check order... failed'
    else
      write (stderr, '(a)') 'Check order... passed'
    end if

  end do

  if (allpass) then
    write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  ! deallocate memory
  deallocate (norm)

  call allocator%release_block(du)
  call allocator%release_block(dv)
  call allocator%release_block(dw)

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)

  ! finalize MPI
  call MPI_Finalize(ierr)

contains
  function dahlquist_rhs(u)
    implicit none

    real(dp) :: dahlquist_rhs
    real(dp), intent(in) :: u

    real(dp) :: lambda = -1.0_dp

    dahlquist_rhs = lambda*u

  end function dahlquist_rhs

  function dahlquist_exact_sol(time)
    implicit none

    real(dp) :: dahlquist_exact_sol
    real(dp), intent(in) :: time

    real(dp) :: lambda = -1.0_dp

    dahlquist_exact_sol = exp(lambda*time)

  end function dahlquist_exact_sol

end program test_omp_adamsbashforth

