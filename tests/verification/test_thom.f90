program test_thom
  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda, backend_sz
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, BC_PERIODIC, BC_DIRICHLET, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test, checkerr

  implicit none

  ! The second-derivative roundoff floor is ~eps/dx^2: at n=1024 this is
  ! ~3e-3 in single precision (~6e-12 in double), so the tolerance must
  ! account for it.
#ifdef SINGLE_PREC
  real(dp), parameter :: residual_tol = 5.0e-2_dp
#else
  real(dp), parameter :: residual_tol = 1.0e-8_dp
#endif

  logical :: allpass = .true.
  integer :: n_glob, n_groups
  integer :: nrank, nproc

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(tdsops_t), allocatable :: periodic_tdsops, dirichlet_tdsops

  call initialise_mpi(nrank, nproc)

  n_glob = 1024
  if (backend_is_cuda) then
    n_groups = 128*128/backend_sz
  else
    n_groups = 64*64/backend_sz
  end if

  mesh = build_mesh(n_glob, n_groups)

  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator

  call backend%alloc_tdsops(periodic_tdsops, n_glob, 2*pi/n_glob, &
                            operation='second-deriv', scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  call backend%alloc_tdsops(dirichlet_tdsops, n_glob, 2*pi/(n_glob - 1), &
                            operation='second-deriv', scheme='compact6', &
                            bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

  if (nrank == 0) then
    print *, trim(runtime%backend_name), 'allocator instantiated'
    print *, trim(runtime%backend_name), 'backend instantiated'
  end if

  call run_case('=== Testing periodic case ===', 2*pi/n_glob, &
                periodic_tdsops, 'thom_periodic')
  call run_case('=== Testing Dirichlet case ===', 2*pi/(n_glob - 1), &
                dirichlet_tdsops, 'thom_dirichlet')

  call finalise_test(allpass, nrank)

contains

  function build_mesh(nx, nz) result(mesh_out)
    integer, intent(in) :: nx, nz
    type(mesh_t) :: mesh_out

    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)

    BC_x = ['periodic', 'periodic']
    BC_y = ['periodic', 'periodic']
    BC_z = ['periodic', 'periodic']

    mesh_out = mesh_t([nx, backend_sz, nz], [1, 1, 1], &
                      [2*pi, 1.0_dp, 1.0_dp], &
                      BC_x, BC_y, BC_z)
  end function build_mesh

  subroutine run_case(banner, delta, tdsops, label)
    character(len=*), intent(in) :: banner, label
    real(dp), intent(in) :: delta
    class(tdsops_t), intent(in) :: tdsops

    class(field_t), pointer :: u_field, du_field
    real(dp), allocatable :: u_data(:, :, :), du_data(:, :, :)
    integer :: i, j, k

    if (nrank == 0) print *, banner

    allocate (u_data(backend_sz, n_glob, n_groups))
    allocate (du_data(backend_sz, n_glob, n_groups))

    do k = 1, n_groups
      do j = 1, n_glob
        do i = 1, backend_sz
          u_data(i, j, k) = sin((j - 1)*delta)
        end do
      end do
    end do

    u_field => allocator%get_block(DIR_X, VERT)
    du_field => allocator%get_block(DIR_X, VERT)
    call backend%set_field_data(u_field, u_data, DIR_X)

    call backend%thom_solve(du_field, u_field, tdsops)
    call backend%get_field_data(du_data, du_field, DIR_X)

    call checkerr(u_data, du_data, residual_tol, label, allpass)

    call allocator%release_block(u_field)
    call allocator%release_block(du_field)

  end subroutine run_case

end program test_thom
