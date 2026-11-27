program test_curl_uniform
  !! The curl of a uniform flow must be zero at free-slip walls
  !!
  !! Sets u = 1, v = 0.5 and w = 0.25 on a small mesh and checks that, to
  !! round-off, max|omega| from solver_t%curl (behind the enstrophy in the
  !! monitoring output) is zero, and so is the vorticity magnitude from
  !! compute_derived_fields (the snapshot output). All three components are
  !! non-zero, so every curl derivative acts on a non-zero field.
  !!
  !! The Q-criterion is not checked: it also needs dv/dy and dw/dz, and a
  !! uniform v or w breaks the free-slip condition v = 0 or w = 0 at the
  !! wall, so those derivatives are not zero there.
  !!
  !! A uniform field has no vorticity, whatever the boundary conditions. The
  !! curl differentiates tangential components, which are even across a
  !! free-slip wall, so it must use the symmetric first derivative. With the
  !! plain first derivative the Neumann closure leaves a non-zero curl at the
  !! wall (issue #507).
  !!
  !! Configurations (x is an inflow/outflow direction, as in the cylinder):
  !!   Config 000 : all periodic
  !!   Config 110 : x dirichlet, y neumann, z periodic
  !!   Config 111 : x dirichlet, y neumann, z neumann
  !!
  !! Each run takes one configuration label as its only argument, because
  !! 2DECOMP&FFT allows one mesh per process. The test builds the solver_t by
  !! hand, because solver_t%init reads an input file.

  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda
  use m_common, only: dp, MPI_X3D2_DP, DIR_C, DIR_X, DIR_Y, DIR_Z, VERT
  use m_mesh, only: mesh_t
  use m_postprocess, only: compute_derived_fields
  use m_mpi, only: MPI_COMM_WORLD, MPI_Allreduce, MPI_MAX
  use m_solver, only: solver_t
  use m_tdsops_setup, only: allocate_tdsops
  use m_vector_calculus, only: vector_calculus_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  real(dp), parameter :: TOLERANCE = 1.0e-12_dp

  character(len=9), parameter :: PERIODIC(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: DIRICHLET(2) = ['dirichlet', 'dirichlet']
  character(len=9), parameter :: NEUMANN(2) = ['neumann  ', 'neumann  ']

  integer :: nrank, nproc
  character(len=3) :: label
  logical :: allpass

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  allpass = .true.

  if (command_argument_count() /= 1) &
    error stop 'Usage: test_curl_uniform <000|110|111>'
  call get_command_argument(1, label)

  select case (label)
  case ('000')
    call run_config(label, [32, 32, 32], PERIODIC, PERIODIC, PERIODIC, &
                    allpass)
  case ('110')
    call run_config(label, [33, 33, 32], DIRICHLET, NEUMANN, PERIODIC, &
                    allpass)
  case ('111')
    call run_config(label, [33, 33, 33], DIRICHLET, NEUMANN, NEUMANN, &
                    allpass)
  case default
    error stop 'Unknown configuration: use 000, 110 or 111'
  end select

  call finalise_test(allpass, nrank)

contains

  subroutine run_config(label, dims_global, BC_x, BC_y, BC_z, allpass)
    character(len=*), intent(in) :: label
    integer, intent(in) :: dims_global(3)
    character(len=*), intent(in) :: BC_x(2), BC_y(2), BC_z(2)
    logical, intent(inout) :: allpass

    type(backend_runtime_t), target :: runtime
    type(mesh_t), target :: mesh
    type(solver_t) :: solver
    class(field_t), pointer :: u, v, w, o_i, o_j, o_k

    real(dp), allocatable :: host(:, :, :)
    real(dp) :: local_max, global_max
    integer :: dims(3), dims_padded(3), ierr
    logical :: passed

    mesh = mesh_t(dims_global, [1, 1, nproc], [1.0_dp, 1.0_dp, 1.0_dp], &
                  BC_x, BC_y, BC_z, use_2decomp=.not. backend_is_cuda)

    call runtime%init(mesh)
    solver%backend => runtime%backend
    solver%mesh => mesh
    solver%host_allocator => runtime%host_allocator

    allocate (solver%xdirps, solver%ydirps, solver%zdirps)
    solver%xdirps%dir = DIR_X
    solver%ydirps%dir = DIR_Y
    solver%zdirps%dir = DIR_Z
    call allocate_tdsops(solver%xdirps, solver%backend, mesh, &
                         'compact6', 'compact6', 'classic', 'compact6')
    call allocate_tdsops(solver%ydirps, solver%backend, mesh, &
                         'compact6', 'compact6', 'classic', 'compact6')
    call allocate_tdsops(solver%zdirps, solver%backend, mesh, &
                         'compact6', 'compact6', 'classic', 'compact6')
    solver%vector_calculus = vector_calculus_t(solver%backend)

    u => solver%backend%allocator%get_block(DIR_X, VERT)
    v => solver%backend%allocator%get_block(DIR_X, VERT)
    w => solver%backend%allocator%get_block(DIR_X, VERT)
    o_i => solver%backend%allocator%get_block(DIR_X, VERT)
    o_j => solver%backend%allocator%get_block(DIR_X, VERT)
    o_k => solver%backend%allocator%get_block(DIR_X, VERT)
    call u%fill(1.0_dp)
    call v%fill(0.5_dp)
    call w%fill(0.25_dp)

    call solver%curl(o_i, o_j, o_k, u, v, w)

    dims = mesh%get_dims(VERT)
    dims_padded = solver%backend%allocator%get_padded_dims(DIR_C)
    allocate (host(dims_padded(1), dims_padded(2), dims_padded(3)))

    local_max = 0.0_dp
    call solver%backend%get_field_data(host, o_i)
    local_max = max(local_max, maxval(abs(host(1:dims(1), 1:dims(2), &
                                               1:dims(3)))))
    call solver%backend%get_field_data(host, o_j)
    local_max = max(local_max, maxval(abs(host(1:dims(1), 1:dims(2), &
                                               1:dims(3)))))
    call solver%backend%get_field_data(host, o_k)
    local_max = max(local_max, maxval(abs(host(1:dims(1), 1:dims(2), &
                                               1:dims(3)))))

    ! The snapshot vorticity magnitude reads solver%u, v, w
    solver%u => u
    solver%v => v
    solver%w => w
    call compute_derived_fields(solver, .true., .false.)
    call solver%backend%get_field_data(host, solver%vort)
    local_max = max(local_max, maxval(abs(host(1:dims(1), 1:dims(2), &
                                               1:dims(3)))))

    call MPI_Allreduce(local_max, global_max, 1, MPI_X3D2_DP, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)

    ! A NaN fails the comparison, so it counts as a failure
    passed = global_max < TOLERANCE
    if (nrank == 0) then
      write (stderr, '(A,A,A,ES12.4,A,ES9.2,2X,A)') &
        '  Config ', label, ': max|omega| = ', global_max, &
        '  tolerance = ', TOLERANCE, merge('PASS', 'FAIL', passed)
    end if
    if (.not. passed) allpass = .false.

    call solver%backend%allocator%release_block(u)
    call solver%backend%allocator%release_block(v)
    call solver%backend%allocator%release_block(w)
    call solver%backend%allocator%release_block(o_i)
    call solver%backend%allocator%release_block(o_j)
    call solver%backend%allocator%release_block(o_k)
    call solver%backend%allocator%release_block(solver%vort)
  end subroutine run_config

end program test_curl_uniform
