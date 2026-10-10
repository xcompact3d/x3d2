program perf_transeq
  !! Benchmark of a single distributed transport-equation component
  !! (backend%transeq_species with the velocity halo exchange), on a domain
  !! decomposed in x.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier, MPI_Wtime

  use m_allocator, only: allocator_t, field_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_device, &
                               backend_label, backend_sz
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_solver, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi, finalise_test, report_perf_minmax

  implicit none

  integer, parameter :: n_glob = 512, n_warmup = 10
  real(dp), parameter :: periodic_bw = 16.0_dp, nu = 1.0_dp
  integer :: n, n_groups, n_iters, ndof
  integer :: nrank, nproc, ierr
  integer :: mem_clock_rt, mem_bus_width
  logical :: has_device_bw_info

  type(backend_runtime_t), target :: runtime
  type(mesh_t), target :: mesh
  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(dirps_t) :: xdirps
  class(field_t), pointer :: u, spec, dspec

  call initialise_mpi(nrank, nproc)
  if (nrank == 0) print *, 'Performance benchmark with', nproc, 'ranks'

  if (backend_is_device) then
    n_groups = 512*512/backend_sz
    n_iters = 100
  else
    n_groups = 128*128/backend_sz
    n_iters = 100
  end if

  mesh = mesh_t([n_glob, backend_sz, n_groups], [nproc, 1, 1], &
                [2*pi, 1._dp, 1._dp], &
                ['periodic', 'periodic'], ['periodic', 'periodic'], &
                ['periodic', 'periodic'])
  call runtime%init(mesh)
  backend => runtime%backend
  allocator => runtime%allocator

  xdirps%dir = DIR_X
  call allocate_tdsops(xdirps, backend, mesh, &
                       'compact6', 'compact6', 'classic', 'compact6')

  n = mesh%get_n(DIR_X, VERT)
  ndof = n*n_groups*backend_sz

  u => allocator%get_block(DIR_X, VERT)
  spec => allocator%get_block(DIR_X, VERT)
  dspec => allocator%get_block(DIR_X, VERT)

  call backend%get_device_bw_info(mem_clock_rt, mem_bus_width, &
                                  has_device_bw_info)

  call run_case('periodic', periodic_bw)

  call allocator%release_block(u)
  call allocator%release_block(spec)
  call allocator%release_block(dspec)
  call finalise_test(.true., nrank)

contains

  subroutine run_case(case_name, consumed_bw)
    character(len=*), intent(in) :: case_name
    real(dp), intent(in) :: consumed_bw

    real(dp), allocatable :: data(:, :, :)
    real(dp) :: dx, tstart, tend
    integer :: iter, j, dims(3), offset
    character(len=:), allocatable :: label

    dims = u%get_shape()
    allocate (data(dims(1), dims(2), dims(3)))
    dx = mesh%geo%d(DIR_X)
    offset = mesh%par%nrank_dir(DIR_X)*n

    data = 0._dp
    do j = 1, n
      data(:, j, :) = sin((j - 1 + offset)*dx)
    end do
    call backend%set_field_data(u, data, DIR_X)
    do j = 1, n
      data(:, j, :) = cos((j - 1 + offset)*dx)
    end do
    call backend%set_field_data(spec, data, DIR_X)

    do iter = 1, n_warmup
      call backend%transeq_species(dspec, u, spec, nu, xdirps, .true.)
    end do
    call backend%sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    tstart = MPI_Wtime()
    do iter = 1, n_iters
      call backend%transeq_species(dspec, u, spec, nu, xdirps, .true.)
    end do
    call backend%sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    tend = MPI_Wtime()

    label = trim(backend_label)//'_transeq_'//case_name
    if (has_device_bw_info) then
      call report_perf_minmax(label, tend - tstart, n_iters, ndof, &
                              consumed_bw, mem_clock_rt, mem_bus_width)
    else
      call report_perf_minmax(label, tend - tstart, n_iters, ndof, &
                              consumed_bw)
    end if
  end subroutine run_case

end program perf_transeq
