program test_filter
  !! Verifies the tridiagonal low-pass filter against its defining properties.
  !!
  !! The filter exists to remove the 2*dx mode, which the compact schemes
  !! cannot dissipate and the staggered pressure projection cannot see. So the
  !! properties that matter are its transfer function at the two ends of the
  !! spectrum: it must leave smooth fields alone and annihilate the sawtooth.
  !!
  !! The zero-wavenumber check also pins down how large alpha may be. The
  !! filter's tridiagonal system has diagonal 1 and off-diagonals alpha, so it
  !! is only marginally diagonally dominant as alpha approaches 0.5, and the
  !! distributed solver's truncation then shows up as a DC error: at 32 points
  !! per rank it is 1e-9 for alpha=0.4 but 1e-3 for Incompact3d's 0.49.
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, &
                      RDR_X2Y, RDR_Y2X
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
  use m_tdsops, only: dirps_t

#ifdef CUDA
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  type(mesh_t), target :: mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  type(dirps_t), target :: xdirps, ydirps, zdirps
  class(field_t), pointer :: f, filtered

  integer, parameter :: dims_global(3) = [32, 33, 32]
  integer, parameter :: nproc_dir(3) = [1, 1, 1]
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  real(dp), parameter :: filter_alpha = 0.4_dp
  character(len=9), parameter :: bc_per(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_slip(2) = ['neumann  ', 'neumann  ']

  real(dp), allocatable :: d(:, :, :), before(:, :, :)
  real(dp) :: tolerance, dc_tolerance, err, amp_before, amp_after
  integer :: dims(3), dpad(3), i, j, k, ierr
  logical :: all_pass

  call MPI_Init(ierr)
  all_pass = .true.
  tolerance = 1000._dp*epsilon(1._dp)
  ! Zero wavenumber is where the distributed solver's truncation shows, so it
  ! cannot be held to machine precision. This bound still catches a wrong
  ! stencil, which would miss unity by order 1e-2 or more.
  dc_tolerance = 1e-8_dp

  mesh = mesh_t(dims_global, nproc_dir, lengths, bc_per, bc_slip, bc_per)
  dims = mesh%get_dims(VERT)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(dims, SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(dims, SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  xdirps%dir = DIR_X; ydirps%dir = DIR_Y; zdirps%dir = DIR_Z
  call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6', filter_alpha=filter_alpha)
  call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6', filter_alpha=filter_alpha)
  call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
                       'classic', 'compact6', filter_alpha=filter_alpha)

  f => allocator%get_block(DIR_X, VERT)
  filtered => allocator%get_block(DIR_X, VERT)
  dpad = allocator%get_padded_dims(DIR_C)
  allocate (d(dpad(1), dpad(2), dpad(3)), before(dpad(1), dpad(2), dpad(3)))

  ! --- A constant must pass through untouched -----------------------------
  ! The stencil weights sum to one, so the transfer function is exactly 1 at
  ! zero wavenumber. This also catches a sign slip in any boundary row.
  d = 3.75_dp
  call backend%set_field_data(f, d)
  call backend%tds_solve(filtered, f, xdirps%lowpass)
  call backend%get_field_data(d, filtered)
  err = maxval(abs(d(1:dims(1), 1:dims(2), 1:dims(3)) - 3.75_dp))
  call check('constant preserved in x', err, dc_tolerance, all_pass)

  ! In y the operator acts along the y-pencil, so the field is reordered
  ! first. The odd operator pins the boundary to zero, so only the even one
  ! can carry a constant through a free-slip boundary.
  block
    class(field_t), pointer :: f_y, filt_y
    f_y => allocator%get_block(DIR_Y)
    filt_y => allocator%get_block(DIR_Y)
    d = 3.75_dp
    call backend%set_field_data(f, d)
    call backend%reorder(f_y, f, RDR_X2Y)
    call backend%tds_solve(filt_y, f_y, ydirps%lowpass_sym)
    call backend%reorder(filtered, filt_y, RDR_Y2X)
    call backend%get_field_data(d, filtered)
    err = maxval(abs(d(1:dims(1), 1:dims(2), 1:dims(3)) - 3.75_dp))
    call check('constant preserved in y (even)', err, dc_tolerance, all_pass)
    call allocator%release_block(f_y)
    call allocator%release_block(filt_y)
  end block

  ! --- The 2*dx sawtooth must be annihilated ------------------------------
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        d(i, j, k) = real((-1)**i, dp)
      end do
    end do
  end do
  call backend%set_field_data(f, d)
  call backend%tds_solve(filtered, f, xdirps%lowpass)
  call backend%get_field_data(d, filtered)
  amp_after = maxval(abs(d(1:dims(1), 1:dims(2), 1:dims(3))))
  call check('2dx mode removed in x', amp_after, 1e-10_dp, all_pass)

  ! --- A well-resolved mode must survive ----------------------------------
  ! Filtering is only acceptable if it leaves the physics alone; an 8-point
  ! wave should lose no more than a few percent.
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        d(i, j, k) = sin(2._dp*pi*real(i - 1, dp)/8._dp)
      end do
    end do
  end do
  before = d
  amp_before = maxval(abs(before(1:dims(1), 1:dims(2), 1:dims(3))))
  call backend%set_field_data(f, d)
  call backend%tds_solve(filtered, f, xdirps%lowpass)
  call backend%get_field_data(d, filtered)
  amp_after = maxval(abs(d(1:dims(1), 1:dims(2), 1:dims(3))))
  if (abs(amp_after - amp_before)/amp_before > 0.05_dp) then
    print *, 'FAIL: 8dx mode damped by more than 5%, amplitude ', &
      amp_before, ' -> ', amp_after
    all_pass = .false.
  else
    print *, 'PASS: 8dx mode preserved, amplitude ratio ', &
      amp_after/amp_before
  end if

  call allocator%release_block(f)
  call allocator%release_block(filtered)
  deallocate (d, before)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  subroutine check(label, error, tol, pass)
    character(*), intent(in) :: label
    real(dp), intent(in) :: error, tol
    logical, intent(inout) :: pass

    if (error > tol) then
      print *, 'FAIL: ', label, ' error=', error, ' tolerance=', tol
      pass = .false.
    else
      print *, 'PASS: ', label, ' error=', error
    end if
  end subroutine check

end program test_filter
