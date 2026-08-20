program perf_thom
  use mpi, only: MPI_Wtime

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda, backend_sz
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, BC_PERIODIC, BC_DIRICHLET, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test, &
                          write_perf_metric, write_perf_summary, &
                          write_device_bw_metric

  implicit none

  integer :: n_glob, n, n_groups
  integer :: n_iters, n_warmup
  integer :: ndof
  integer :: nrank, nproc
  integer :: mem_clock_rt, mem_bus_width
  logical :: has_device_bw_info
  real(dp) :: dx_per, dx
  real(dp) :: periodic_bw, dirichlet_bw
  character(len=4) :: backend_label

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(field_t), pointer :: u_field, du_field
  class(tdsops_t), allocatable :: periodic_tdsops, dirichlet_tdsops
  real(dp), allocatable :: u_data(:, :, :)

  call initialise_mpi(nrank, nproc)
  call configure_benchmark()

  mesh = build_mesh(n_glob, n_groups)
  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator

  u_field => allocator%get_block(DIR_X, VERT)
  du_field => allocator%get_block(DIR_X, VERT)
  allocate (u_data(backend_sz, n, n_groups))

  call backend%get_device_bw_info(mem_clock_rt, mem_bus_width, &
                                  has_device_bw_info)

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  call backend%alloc_tdsops(periodic_tdsops, n, dx_per, &
                            operation='second-deriv', scheme='compact6', &
                            bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  call backend%alloc_tdsops(dirichlet_tdsops, n, dx, &
                            operation='second-deriv', scheme='compact6', &
                            bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

  call run_case('periodic', dx_per, periodic_tdsops, periodic_bw)
  call run_case('dirichlet', dx, dirichlet_tdsops, dirichlet_bw)

  if (has_device_bw_info) then
    call write_device_bw_metric(mem_clock_rt, mem_bus_width)
  end if

  call allocator%release_block(u_field)
  call allocator%release_block(du_field)
  call finalise_test(.true., nrank)

contains

  subroutine configure_benchmark()
    n_glob = 512
    if (backend_is_cuda) then
      backend_label = 'cuda'
      periodic_bw = 6.0_dp
      dirichlet_bw = 4.0_dp
      n_groups = 512*512/backend_sz
      n_iters = 1000
    else
      backend_label = 'omp '
      periodic_bw = 3.0_dp
      dirichlet_bw = 3.0_dp
      n_groups = 128*128/backend_sz
      n_iters = 500
    end if
    n_warmup = 10
    n = n_glob
    ndof = n_glob*n_groups*backend_sz
  end subroutine configure_benchmark

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

  subroutine run_case(case_name, delta, tdsops, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: delta, consumed_bw
    class(tdsops_t), intent(in) :: tdsops

    integer :: iter
    real(dp) :: tstart, tend

    call initialise_input(delta)
    call backend%set_field_data(u_field, u_data, DIR_X)

    do iter = 1, n_warmup
      call backend%thom_solve(du_field, u_field, tdsops)
    end do
    call backend%sync()

    tstart = MPI_Wtime()
    do iter = 1, n_iters
      call backend%thom_solve(du_field, u_field, tdsops)
    end do
    call backend%sync()
    tend = MPI_Wtime()

    call write_perf_metric(trim(backend_label)//'_thom_'//trim(case_name), &
                           tend - tstart, n_iters, ndof, consumed_bw)
    if (has_device_bw_info) then
      call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw, &
                              mem_clock_rt, mem_bus_width)
    else
      call write_perf_summary(tend - tstart, n_iters, ndof, consumed_bw)
    end if
  end subroutine run_case

  subroutine initialise_input(delta)
    real(dp), intent(in) :: delta
    integer :: i, j, k

    do k = 1, n_groups
      do j = 1, n
        do i = 1, backend_sz
          u_data(i, j, k) = sin((j - 1)*delta)
        end do
      end do
    end do
  end subroutine initialise_input

end program perf_thom
