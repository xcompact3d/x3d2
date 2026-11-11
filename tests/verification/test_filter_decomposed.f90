program test_filter_decomposed
  !! The low-pass filter must give the same result whether or not its
  !! direction is decomposed.
  !!
  !! x and z get the same boundary conditions and the same number of points,
  !! and only z is decomposed. The same profile is then filtered along x,
  !! where Thomas solves the filter exactly, and along z, where the
  !! distributed solver does with its reduced system solved exactly. The two
  !! results must agree point for point. Before the reduced system was solved
  !! exactly, they differed by about 1e-2 at Incompact3d's filter_alpha = 0.49
  !! and 16 points per rank, even for a constant, which grew the mean flow of
  !! a z-decomposed ABL run exponentially.
  !!
  !! Argument 1 selects the boundary conditions of x and z: 'periodic' (the
  !! ABL's spanwise direction) or 'walls' (Dirichlet start, Neumann end).
  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_size, MPI_Init, MPI_Finalize, &
                   MPI_Allreduce, MPI_IN_PLACE, MPI_MAX

  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, DIR_Z, DIR_C, VERT, RDR_X2Z, RDR_Z2X, &
                      MPI_X3D2_DP
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t, tdsops_t
  use m_tdsops_setup, only: allocate_tdsops

  implicit none

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  type(dirps_t), target :: xdirps, zdirps
  class(field_t), pointer :: f, filtered, f_z, filt_z

  real(dp), parameter :: filter_alpha = 0.49_dp
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  character(len=9), parameter :: bc_per(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_walls(2) = ['dirichlet', 'neumann  ']

  real(dp), allocatable :: d(:, :, :), ref(:)
  character(len=16) :: mode
  character(len=9) :: bc_xz(2)
  integer :: dims_global(3), dims(3), dpad(3), vert_global(3), nproc, ierr
  integer :: z_off, n
  logical :: all_pass

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  all_pass = .true.

  call get_command_argument(1, mode)
  select case (trim(mode))
  case ('periodic')
    bc_xz = bc_per
  case ('walls')
    bc_xz = bc_walls
  case default
    error stop 'argument 1 must be periodic or walls'
  end select

  ! 16 z points per rank, where the truncation was 1.6e-2
  dims_global = [16*nproc, 4, 16*nproc]
  mesh = mesh_t(dims_global, [1, 1, nproc], lengths, bc_xz, bc_per, bc_xz)
  dims = mesh%get_dims(VERT)
  vert_global = mesh%get_global_dims(VERT)
  n = vert_global(1)
  z_off = mesh%par%n_offset(3)
  if (vert_global(3) /= n) &
    error stop 'x and z must have the same number of vertices'

  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend

  xdirps%dir = DIR_X; zdirps%dir = DIR_Z
  call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6', filter_alpha=filter_alpha)
  call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6', filter_alpha=filter_alpha)

  f => allocator%get_block(DIR_X, VERT)
  filtered => allocator%get_block(DIR_X, VERT)
  f_z => allocator%get_block(DIR_Z)
  filt_z => allocator%get_block(DIR_Z)
  dpad = allocator%get_padded_dims(DIR_C)
  allocate (d(dpad(1), dpad(2), dpad(3)), ref(n))

  if (mesh%par%is_root()) print *, 'mode ', trim(mode), ', ranks along z ', &
    nproc
  call compare('constant, odd operator ', constant, xdirps%lowpass, &
               zdirps%lowpass)
  call compare('constant, even operator', constant, xdirps%lowpass_sym, &
               zdirps%lowpass_sym)
  call compare('profile, odd operator  ', profile, xdirps%lowpass, &
               zdirps%lowpass)
  call compare('profile, even operator ', profile, xdirps%lowpass_sym, &
               zdirps%lowpass_sym)

  call allocator%release_block(f)
  call allocator%release_block(filtered)
  call allocator%release_block(f_z)
  call allocator%release_block(filt_z)

  if (.not. all_pass) error stop 'FAIL'
  if (mesh%par%is_root()) print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  pure real(dp) function constant(m)
    integer, intent(in) :: m

    constant = 3.75_dp + 0._dp*m
  end function constant

  pure real(dp) function profile(m)
    !! Every scale from the mean to the 2*dx mode
    integer, intent(in) :: m

    profile = 3.75_dp + sin(2._dp*pi*3._dp*m/n) &
              + 0.3_dp*cos(2._dp*pi*7._dp*m/n) + 0.5_dp*(-1)**m &
              + 0.1_dp*modulo(37*m, 11)
  end function profile

  subroutine compare(label, g, op_x, op_z)
    !! Filter g along x (undecomposed) and along z (decomposed) and compare
    character(*), intent(in) :: label
    interface
      pure real(dp) function g(m)
        import :: dp
        integer, intent(in) :: m
      end function g
    end interface
    class(tdsops_t), intent(in) :: op_x, op_z

    real(dp) :: err, scale
    integer :: i, k

    ! Along x: every rank holds the whole line
    d = 0._dp
    do i = 1, dims(1)
      d(i, 1:dims(2), 1:dims(3)) = g(i)
    end do
    call backend%set_field_data(f, d)
    call backend%tds_solve(filtered, f, op_x)
    call backend%get_field_data(d, filtered)
    ref = d(1:n, 1, 1)

    ! Along z: each rank holds its block of the line
    d = 0._dp
    do k = 1, dims(3)
      d(1:dims(1), 1:dims(2), k) = g(z_off + k)
    end do
    call backend%set_field_data(f, d)
    call backend%reorder(f_z, f, RDR_X2Z)
    call backend%tds_solve(filt_z, f_z, op_z)
    call backend%reorder(filtered, filt_z, RDR_Z2X)
    call backend%get_field_data(d, filtered)

    err = 0._dp
    do k = 1, dims(3)
      err = max(err, maxval(abs(d(1:dims(1), 1:dims(2), k) - ref(z_off + k))))
    end do
    call MPI_Allreduce(MPI_IN_PLACE, err, 1, MPI_X3D2_DP, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    ! The system is solved exactly either way, so only round-off separates
    ! the two: a few epsilon of the O(5) field, amplified by the conditioning
    ! of the near-singular filter system
    scale = maxval(abs(ref))
    if (err > 1.e3_dp*epsilon(1._dp)*scale) then
      if (mesh%par%is_root()) print *, 'FAIL: ', label, ' max difference ', err
      all_pass = .false.
    else
      if (mesh%par%is_root()) print *, 'PASS: ', label, ' max difference ', err
    end if
  end subroutine compare

end program test_filter_decomposed
