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

  use iso_fortran_env, only: stderr => error_unit

  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size, &
                   MPI_Finalize, MPI_Init, MPI_Barrier, MPI_Wtime

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, CELL, VERT
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t
  use m_tdsops_setup, only: allocate_tdsops

#ifdef CUDA
  use cudafor
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#elif defined(OMP_TGT)
  use omp_lib, only: omp_get_num_devices, omp_set_default_device, &
                     omp_get_default_device
  use m_omptgt_common, only: SZ
  use m_omptgt_allocator, only: omptgt_allocator_t
  use m_omptgt_backend, only: omptgt_backend_t
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  ! Wavenumber of the default right-hand side.
  integer, parameter :: N_WAVE = 2

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  type(dirps_t), pointer :: xdirps, ydirps, zdirps

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#elif defined(OMP_TGT)
  type(omptgt_backend_t), target :: omptgt_backend
  type(omptgt_allocator_t), target :: omptgt_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
#endif
  type(allocator_t), target :: omp_allocator

  character(len=8) :: backend_name
  integer :: dims_global(3), nproc_dir(3), vert_dims(3), dims(3)
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

#ifdef CUDA
  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                use_2decomp=.false.)
#else
  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                use_2decomp=.true.)
#endif

  ! Get the local vertex dimensions.
  vert_dims = mesh%get_dims(VERT)

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
  backend_name = "CUDA"

  cuda_allocator = cuda_allocator_t(vert_dims, SZ)
  allocator => cuda_allocator

  omp_allocator = allocator_t(vert_dims, SZ)
  host_allocator => omp_allocator

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#elif defined(OMP_TGT)
  ndevs = omp_get_num_devices()
  call omp_set_default_device(mod(nrank, ndevs)) ! round-robin
  devnum = omp_get_default_device()
  backend_name = "OMP_TGT"

  omptgt_allocator = omptgt_allocator_t(vert_dims, SZ)
  allocator => omptgt_allocator

  omp_allocator = allocator_t(vert_dims, SZ)
  host_allocator => omp_allocator

  omptgt_backend = omptgt_backend_t(mesh, allocator)
  backend => omptgt_backend
#else
  backend_name = "OMP"

  omp_allocator = allocator_t(vert_dims, SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  ! Get the local cell dimensions.
  ! The loops in `fill_rhs` use them.
  dims = mesh%get_dims(CELL)

  ! ---- tdsops and Poisson FFT setup ----
  allocate (xdirps, ydirps, zdirps)
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

  min_time = huge(1.0_dp)
  sum_time = 0.0_dp

  do irep = 1, repeat_count
    call fill_rhs(host_field)
    call backend%set_field_data(f_device, host_field%data, DIR_C)
    call f_device%set_data_loc(CELL)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t0 = MPI_Wtime()
    call backend%poisson_fft%solve_poisson(f_device, temp)
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

end program x3d2_poisson_solver
