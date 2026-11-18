program test_abl_shifted_periodic
  !! The ABL shifted periodic fringe must reproduce Incompact3d's
  !! fringe_region (ishiftedper = 1, iconcprec = 0) exactly. The reference
  !! below is the legacy loop transcribed on the global field: the ramp, the
  !! band sizes and the spanwise shift are all recomputed here rather than
  !! taken from m_abl. z is decomposed over every rank, so on several ranks
  !! part of the shift crosses rank boundaries.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Comm_size, MPI_Init, MPI_Finalize

  use m_abl, only: abl_t
  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, DIR_C, VERT
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  type(allocator_t), pointer :: host_allocator
  class(field_t), pointer :: u, v, w
  type(abl_config_t) :: cfg
  type(abl_t) :: abl

  real(dp), parameter :: lengths(3) = [1.44_dp, 0.46_dp, 0.72_dp]
  character(len=9), parameter :: bc_periodic(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_wall(2) = ['dirichlet', 'neumann  ']

  real(dp), allocatable :: host(:, :, :), reference(:, :, :)
  real(dp) :: max_error
  integer :: dims_global(3), nproc_dir(3), dims(3), dims_padded(3)
  integer :: nproc, ierr, comp
  logical :: all_pass

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  all_pass = .true.

  dims_global = [64, 9, 16*nproc]
  nproc_dir = [1, 1, nproc]
  mesh = mesh_t(dims_global, nproc_dir, lengths, &
                bc_periodic, bc_wall, bc_periodic)
  dims = mesh%get_dims(VERT)

  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend
  host_allocator => runtime%host_allocator

  dims_padded = allocator%get_padded_dims(DIR_C)
  allocate (host(dims_padded(1), dims_padded(2), dims_padded(3)))

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  call set_field(u, 1)
  call set_field(v, 2)
  call set_field(w, 3)

  call cfg%read(nml_string='&abl_nml u_star = 0.11, pressure_gradient = &
                &.true., shifted_periodic = .true. /')
  abl = abl_t(backend, mesh, host_allocator, cfg)
  call abl%configure_shifted_periodic()
  call abl%apply_shifted_periodic(u, v, w)

  do comp = 1, 3
    select case (comp)
    case (1)
      call backend%get_field_data(host, u)
    case (2)
      call backend%get_field_data(host, v)
    case (3)
      call backend%get_field_data(host, w)
    end select
    reference = legacy_fringe(comp)
    max_error = maxval(abs(host(1:dims(1), 1:dims(2), 1:dims(3)) - reference))
    print *, 'component', comp, 'max error', max_error
    ! The field is O(10), so a last-bit difference in the blend is ~2e-15
    if (max_error > 10._dp*epsilon(1._dp)*maxval(abs(reference))) &
      all_pass = .false.
  end do

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  pure real(dp) function initial(comp, i, j, k)
    !! A distinct value at every global vertex and component
    integer, intent(in) :: comp, i, j, k

    initial = real(comp, dp) + 0.1_dp*real(i, dp) + 0.01_dp*real(j, dp) &
              + sin(2._dp*pi*real(k, dp)/real(dims_global(3), dp))
  end function initial

  subroutine set_field(f, comp)
    class(field_t), intent(inout) :: f
    integer, intent(in) :: comp

    integer :: i, j, k

    host = 0._dp
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          host(i, j, k) = initial(comp, i, j, mesh%par%n_offset(3) + k)
        end do
      end do
    end do
    call backend%set_field_data(f, host)
  end subroutine set_field

  function legacy_fringe(comp) result(local)
    !! Incompact3d fringe_region on the global field, returning this rank's
    !! z block
    integer, intent(in) :: comp
    real(dp), allocatable :: local(:, :, :)

    real(dp), allocatable :: ux(:, :, :), ux_s(:, :, :)
    real(dp) :: xlx, zlz, dx, dz, dshift, x, edl, frl, fre, frs, frd, lambda
    integer :: nx, ny, nz, nshift, npe, nfe, i, j, k

    nx = dims_global(1)
    ny = dims_global(2)
    nz = dims_global(3)
    xlx = lengths(1)
    zlz = lengths(3)
    dx = xlx/real(nx, dp)
    dz = zlz/real(nz, dp)
    allocate (ux(nx, ny, nz), ux_s(nx, ny, nz))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ux(i, j, k) = initial(comp, i, j, k)
        end do
      end do
    end do

    dshift = zlz/8
    nshift = int(dshift/dz)
    do k = 1, nz - nshift
      ux_s(:, :, k + nshift) = ux(:, :, k)
    end do
    do k = 1, nshift
      ux_s(:, :, k) = ux(:, :, nz - nshift + k)
    end do

    edl = 2._dp/3._dp
    frl = 1._dp/6._dp
    fre = xlx*edl
    npe = nx
    frs = fre - fre*frl
    frd = fre - frs
    nfe = int(nx*fre/xlx)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nfe
          x = real(i - 1, dp)*dx
          if (x < frs) then
            lambda = 0._dp
          else if ((x >= frs) .and. (x < (fre - frd/4._dp))) then
            lambda = 0.5_dp*(1._dp - cos(4._dp*pi/3._dp*(x - frs)/frd))
          else if ((x >= (fre - frd/4._dp)) .and. (x < fre)) then
            lambda = 1._dp
          else
            lambda = 0._dp
          end if
          ux(i + npe - nfe, j, k) = lambda*ux_s(i, j, k) &
                                    + (1._dp - lambda)*ux(i + npe - nfe, j, k)
        end do
      end do
    end do

    local = ux(:, :, mesh%par%n_offset(3) + 1:mesh%par%n_offset(3) + dims(3))
  end function legacy_fringe

end program test_abl_shifted_periodic
