program test_omp_adamsbashforth
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_common, only: dp, globs_t, DIR_X
  use m_allocator, only: allocator_t, field_t
  use m_omp_backend, only: omp_backend_t, base_backend_t
  use m_time_integrator, only: time_intg_t

  implicit none

  logical :: allpass = .true.
  integer :: i, j, k, stage, istartup, ierr
  integer :: nstep0 = 64, nstep, nrun = 4, nmethod = 8
  real(dp), allocatable, dimension(:) :: err
  real(dp), allocatable, dimension(:) :: norm
  real(dp) :: dt0 = 0.01_dp, dt, order
  class(field_t), pointer :: u, v, w
  class(field_t), pointer :: du, dv, dw

  type(globs_t) :: globs
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(omp_backend_t), target :: omp_backend
  type(allocator_t), target :: omp_allocator
  class(time_intg_t), allocatable :: time_integrator

  ! initialize MPI
  call MPI_Init(ierr)

  ! set globs parameters
  globs%nx = 1
  globs%ny = 1
  globs%nz = 1

  globs%nx_loc = globs%nx
  globs%ny_loc = globs%ny
  globs%nz_loc = globs%nz

  globs%n_groups_x = globs%ny_loc*globs%nz_loc
  globs%n_groups_y = globs%nx_loc*globs%nz_loc
  globs%n_groups_z = globs%nx_loc*globs%ny_loc

  ! allocate object
  omp_allocator = allocator_t(globs%nx, globs%ny, globs%nz, 1)
  allocator => omp_allocator
  print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(globs, allocator)
  backend => omp_backend
  print *, 'OpenMP backend instantiated'

  ! allocate memory
  u => allocator%get_block(DIR_X)
  v => allocator%get_block(DIR_X)
  w => allocator%get_block(DIR_X)

  du => allocator%get_block(DIR_X)
  dv => allocator%get_block(DIR_X)
  dw => allocator%get_block(DIR_X)

  allocate (norm(nrun))

  ! compute l2 norm for various step sizes
  do k = 1, nmethod

    ! initialize time-integrator
    time_integrator = time_intg_t(allocator=allocator, &
                                  backend=backend, method=k)

    dt = dt0
    nstep = nstep0
    do j = 1, nrun
      ! initial condition
      time_integrator%istep = 1

      ! compute l2 norm for a given step size
      allocate (err(nstep))
      u%data(1, 1, 1) = 1.0_dp

      ! startup
      istartup = time_integrator%nstep - 1
      do i = 1, istartup
        du%data(1, 1, 1) = dahlquist_rhs(u%data(1, 1, 1))
        call time_integrator%step(u, v, w, du, dv, dw, dt)
        u%data(1, 1, 1) = dahlquist_exact_sol(real(i, dp)*dt)
      end do

      ! post-startup
      do i = 1, nstep
        do stage = 1, time_integrator%nstage
          du%data(1, 1, 1) = dahlquist_rhs(u%data(1, 1, 1))
          call time_integrator%step(u, v, w, du, dv, dw, dt)
        end do
        err(i) = u%data(1, 1, 1) - dahlquist_exact_sol( &
                 real(i + istartup, dp)*dt)
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
    if (abs(order - real(time_integrator%order, dp)) > 0.25_dp) then
      allpass = .false.
      write (stderr, '(a)') 'Check order... failed'
    else
      write (stderr, '(a)') 'Check order... passed'
    end if

    call time_integrator%finalize

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

