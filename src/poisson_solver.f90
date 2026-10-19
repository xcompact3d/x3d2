program x3d2_poisson_solver
  !! Standalone driver for the Poisson solver.
  !! It links only the modules that the Poisson solve needs, from the
  !! `x3d2_poisson` CMake target.
  !! It builds the mesh and the FFT Poisson solver for the grid and the
  !! boundary conditions on the command line.
  !! It fills the right-hand side and solves it the requested number of times.
  !! It then reports the solve time.
  !!
  !! The default right-hand side is `cos(2 pi x)`.
  !! This field has zero mean and zero normal derivative at a non-periodic
  !! wall.
  !! It is valid for every supported combination of boundary conditions.
  !! To use another right-hand side, edit `fill_rhs` and rebuild.
  !!
  !! `dirichlet` is the velocity boundary condition.
  !! The pressure gets homogeneous Neumann walls in that direction.
  !! The test `test_poisson_bc` checks the solve for all four boundary
  !! condition configurations.
  !!
  !! Usage: x3d2_poisson_solver nx ny nz [bc_x bc_y bc_z] [--repeat N]
  !!   bc_x, bc_y, bc_z: 'periodic' (default) or 'dirichlet'.
  !!                     The value applies to both ends of the direction.
  !!                     A Dirichlet direction needs an odd grid size.
  !!   --repeat N:       number of times to solve (default 1).

  use iso_fortran_env, only: stderr => error_unit, real64

  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size, &
                   MPI_Finalize, MPI_Init, MPI_Barrier, MPI_Wtime, &
                   MPI_Allreduce, MPI_IN_PLACE, MPI_SUM, MPI_MAX, &
                   MPI_DOUBLE_PRECISION

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, CELL, MPI_X3D2_DP
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t
  use m_tdsops_setup, only: allocate_tdsops

  use m_backend_setup, only: backend_setup_t, backend_name

  implicit none

  ! Wavenumber of the default right-hand side.
  integer, parameter :: N_WAVE = 2

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  type(dirps_t) :: xdirps, ydirps, zdirps

  type(backend_setup_t), target :: setup

  integer :: dims_global(3), nproc_dir(3), dims(3)
  integer :: nrank, nproc, ierr
  real(dp) :: L_global(3)
  character(len=9) :: bc_x_name, bc_y_name, bc_z_name
  character(len=9) :: BC_x(2), BC_y(2), BC_z(2)
  integer :: repeat_count, irep
  ! `t0` and `t1` hold raw `MPI_Wtime` readings in double precision.
  ! The code subtracts them in double precision.
  ! It then converts the difference to `dp`.
  double precision :: t0, t1
  real(dp) :: elapsed, min_time, sum_time
  real(dp) :: n_pi

  class(field_t), pointer :: f_device, temp
  class(field_t), pointer :: host_field

  ! ---- MPI init ----
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! ---- Argument parsing ----
  call parse_arguments()

  ! ---- Backend and allocator selection ----
  nproc_dir = [1, 1, nproc]
  L_global = [1.0_dp, 1.0_dp, 1.0_dp]

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                use_2decomp=(backend_name == 'OMP'))

  call setup%init(mesh)
  backend => setup%backend
  allocator => setup%allocator
  host_allocator => setup%host_allocator

  ! Get the local cell dimensions.
  ! The loops in `fill_rhs` use them.
  dims = mesh%get_dims(CELL)

  ! ---- tdsops and Poisson FFT setup ----
  xdirps%dir = DIR_X
  ydirps%dir = DIR_Y
  zdirps%dir = DIR_Z
  call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')
  call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')
  call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6')

  call backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)

  ! ---- Default right-hand side ----
  n_pi = real(N_WAVE, dp)*pi

  f_device => backend%allocator%get_block(DIR_C, CELL)
  temp => backend%allocator%get_block(DIR_C)
  host_field => host_allocator%get_block(DIR_C)

  call fill_rhs(host_field)
  call warn_if_rhs_mean_nonzero(host_field)

  min_time = huge(1.0_dp)
  sum_time = 0.0_dp

  do irep = 1, repeat_count
    call backend%set_field_data(f_device, host_field%data, DIR_C)
    call f_device%set_data_loc(CELL)

    call backend%sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t0 = MPI_Wtime()
    call backend%poisson_fft%solve_poisson(f_device, temp)
    call backend%sync()
    t1 = MPI_Wtime()

    elapsed = real(t1 - t0, dp)
    min_time = min(min_time, elapsed)
    sum_time = sum_time + elapsed
  end do

  call backend%allocator%release_block(temp)
  call host_allocator%release_block(host_field)
  call backend%allocator%release_block(f_device)

  ! ---- Report ----
  if (nrank == 0) then
    write (stderr, '(A)') ''
    write (stderr, '(A,I0,A,I0,A,I0)') &
      'Grid: ', dims_global(1), ' x ', dims_global(2), ' x ', dims_global(3)
    write (stderr, '(A,A,A,A,A,A)') &
      'BC: x=', trim(bc_x_name), ' y=', trim(bc_y_name), &
      ' z=', trim(bc_z_name)
    write (stderr, '(A,I0)') 'Ranks: ', nproc
    write (stderr, '(A,A)') 'Backend: ', trim(backend_name)
    write (stderr, '(A,ES12.4,A)') 'Solve time (min):  ', min_time, ' s'
    write (stderr, '(A,ES12.4,A)') 'Solve time (mean): ', &
      sum_time/real(repeat_count, dp), ' s'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine parse_arguments()
    !! Read `nx ny nz`, then the optional `bc_x bc_y bc_z`, then the optional
    !! `--repeat N`.
    !! The `bc` values are `periodic` or `dirichlet`.
    !! The default is `periodic`.
    !! The default for `--repeat` is 1.
    !! Any other input prints one usage line and stops.
    integer :: nargs, iarg, ios
    character(len=32) :: arg

    nargs = command_argument_count()

    bc_x_name = 'periodic'; bc_y_name = 'periodic'; bc_z_name = 'periodic'
    repeat_count = 1

    if (nargs < 3) call usage_error()

    do iarg = 1, 3
      call get_command_argument(iarg, arg)
      read (arg, *, iostat=ios) dims_global(iarg)
      if (ios /= 0 .or. dims_global(iarg) < 1) call usage_error()
    end do

    iarg = 4
    if (nargs >= iarg) then
      call get_command_argument(iarg, arg)
      if (trim(arg) /= '--repeat') then
        if (nargs < iarg + 2) call usage_error()
        call get_command_argument(iarg, bc_x_name)
        call get_command_argument(iarg + 1, bc_y_name)
        call get_command_argument(iarg + 2, bc_z_name)
        if (trim(bc_x_name) /= 'periodic' .and. &
            trim(bc_x_name) /= 'dirichlet') call usage_error()
        if (trim(bc_y_name) /= 'periodic' .and. &
            trim(bc_y_name) /= 'dirichlet') call usage_error()
        if (trim(bc_z_name) /= 'periodic' .and. &
            trim(bc_z_name) /= 'dirichlet') call usage_error()
        iarg = iarg + 3
      end if
    end if

    if (nargs >= iarg) then
      call get_command_argument(iarg, arg)
      if (trim(arg) /= '--repeat') call usage_error()
      if (nargs < iarg + 1) call usage_error()
      call get_command_argument(iarg + 1, arg)
      read (arg, *, iostat=ios) repeat_count
      if (ios /= 0 .or. repeat_count < 1) call usage_error()
      iarg = iarg + 2
    end if

    if (nargs >= iarg) call usage_error()

    if (trim(bc_x_name) == 'dirichlet' .and. mod(dims_global(1), 2) == 0) &
      call usage_error()
    if (trim(bc_y_name) == 'dirichlet' .and. mod(dims_global(2), 2) == 0) &
      call usage_error()
    if (trim(bc_z_name) == 'dirichlet' .and. mod(dims_global(3), 2) == 0) &
      call usage_error()

    BC_x = [bc_x_name, bc_x_name]
    BC_y = [bc_y_name, bc_y_name]
    BC_z = [bc_z_name, bc_z_name]

  end subroutine parse_arguments

  subroutine usage_error()
    if (nrank == 0) then
      write (stderr, '(A)') &
        'Usage: x3d2_poisson_solver nx ny nz [bc_x bc_y bc_z] &
        &[--repeat N]  (bc in {periodic,dirichlet}, dirichlet dims odd)'
    end if
    call MPI_Finalize(ierr)
    error stop 'invalid arguments'
  end subroutine usage_error

  subroutine fill_rhs(host_field)
    !! Fill `host_field` with the right-hand side of the Poisson equation.
    !! The routine uses the host-associated `mesh`, `dims` and `n_pi`.
    !! Edit this routine to supply another right-hand side, then rebuild.
    !! The coordinates are the cell centres from
    !! `mesh%get_coordinates(i, j, k, CELL)`.
    !! The right-hand side must have zero mean.
    !! In every supported combination of boundary conditions, the pressure
    !! has periodic or homogeneous Neumann conditions.
    !! The default is `cos(2 pi x)`.
    class(field_t), intent(inout) :: host_field

    integer :: i, j, k
    real(dp) :: coords(3)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = mesh%get_coordinates(i, j, k, CELL)
          host_field%data(i, j, k) = cos(n_pi*coords(1))
        end do
      end do
    end do
  end subroutine fill_rhs

  subroutine warn_if_rhs_mean_nonzero(host_field)
    !! Warn on rank 0 if the right-hand side does not have zero mean.
    !! The solver needs a zero-mean right-hand side for these boundary
    !! conditions.
    !! It discards the mean component.
    !! The solve runs in both cases.
    class(field_t), intent(in) :: host_field

    real(dp) :: local_max
    real(real64) :: local_sum, global_mean, n_cells_global
    integer :: global_dims(3), ierr

    ! The sum uses double precision for every `dp`.
    ! A `real32` sum of a valid zero-mean field can stall past 2^24 terms.
    ! This gives a false warning in a `SINGLE_PREC` build.
    local_sum = sum(real(host_field%data(1:dims(1), 1:dims(2), 1:dims(3)), &
                         real64))
    local_max = maxval(abs(host_field%data(1:dims(1), 1:dims(2), 1:dims(3))))

    call MPI_Allreduce(MPI_IN_PLACE, local_sum, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, local_max, 1, MPI_X3D2_DP, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)

    ! The cell count uses `real64`.
    ! The `product()` of a default integer array overflows above 2^31 cells.
    global_dims = mesh%get_global_dims(CELL)
    n_cells_global = real(global_dims(1), real64) &
                     *real(global_dims(2), real64) &
                     *real(global_dims(3), real64)
    global_mean = local_sum/n_cells_global

    ! The tolerance scales with the precision of the build through
    ! `epsilon(1.0_dp)`.
    ! A fixed 1e-10 is below the single-precision epsilon of about 1.19e-7.
    ! A `SINGLE_PREC` build rounds the discrete mean of `cos(2 pi x)` to about
    ! 1e-9 to 1e-8.
    ! The fixed tolerance gives a false warning for that rounding.
    ! The tolerance `100*epsilon(1.0_dp)` is far above that rounding.
    ! It is also far above the double-precision rounding of the same mean.
    if (abs(global_mean) > &
        100.0_dp*epsilon(1.0_dp)*max(local_max, tiny(1.0_dp))) then
      if (nrank == 0) then
        write (stderr, '(A,ES12.4,A)') &
          'Warning: right-hand side mean is ', global_mean, &
          '; the solver needs a zero-mean right-hand side for these &
          &boundary conditions, and the mean component is discarded.'
      end if
    end if
  end subroutine warn_if_rhs_mean_nonzero

end program x3d2_poisson_solver
