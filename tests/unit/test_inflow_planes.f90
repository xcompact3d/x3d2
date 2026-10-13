program test_inflow_planes
  !! Precursor inflow planes: sample a known field, write blocks, read them
  !! back on a different rank layout and check every value and time stamp.
  !! The writer decomposes x, so only some ranks own the sampling plane. A
  !! run restarted inside a block must skip that block.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_C, VERT
  use m_config, only: inflow_config_t
  use m_field, only: field_t
  use m_inflow_planes, only: inflow_plane_writer_t, inflow_file_name, &
                             read_inflow_block
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test, global_all

#ifdef CUDA
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  type(mesh_t), target :: mesh, read_mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  type(allocator_t), target :: host_allocator
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(omp_backend_t), target :: omp_backend
#endif
  class(field_t), pointer :: u, v, w

  integer, parameter :: dims_global(3) = [16, 8, 8]
  real(dp), parameter :: lengths(3) = [4._dp, 1._dp, 2._dp]
  real(dp), parameter :: dt = 0.01_dp
  integer, parameter :: x_index = 5, ntimesteps = 3, sample_every = 2, &
                        start_iter = 1
  character(len=*), parameter :: prefix = 'test_inflow', &
                                 restart_prefix = 'test_inflow_restart'
  character(len=9), parameter :: bc_periodic(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_wall(2) = ['dirichlet', 'neumann  ']

  type(inflow_config_t) :: cfg
  type(inflow_plane_writer_t) :: writer
  integer :: dims(3), irank, nproc, ierr
  logical :: allpass, exists

  call initialise_mpi(irank, nproc)
  allpass = .true.
  if (nproc /= 4) error stop 'test_inflow_planes needs 4 ranks'

  mesh = mesh_t(dims_global, [2, 2, 1], lengths, &
                bc_periodic, bc_wall, bc_periodic)
  dims = mesh%get_dims(VERT)

  host_allocator = allocator_t(dims, SZ)
#ifdef CUDA
  cuda_allocator = cuda_allocator_t(dims, SZ)
  allocator => cuda_allocator
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  allocator => host_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)

  cfg%write_planes = .true.
  cfg%plane_x_index = x_index
  cfg%ntimesteps = ntimesteps
  cfg%sample_every = sample_every
  cfg%start_iter = start_iter

  ! Samples at iterations 3, 5, 7 (block 1) and 9, 11, 13 (block 2)
  cfg%inflow_prefix = prefix
  writer = inflow_plane_writer_t(backend, mesh, cfg, dt)
  call advance(writer, 1, 14)

  ! A restart at iteration 5 sees only part of block 1
  cfg%inflow_prefix = restart_prefix
  writer = inflow_plane_writer_t(backend, mesh, cfg, dt)
  call advance(writer, 5, 14)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! Read on a different decomposition: z split across all ranks
  read_mesh = mesh_t(dims_global, [1, 1, 4], lengths, &
                     bc_periodic, bc_wall, bc_periodic)
  call check_block(inflow_file_name(prefix, 1), 3)
  call check_block(inflow_file_name(prefix, 2), 9)
  call check_block(inflow_file_name(restart_prefix, 2), 9)

  inquire (file=inflow_file_name(restart_prefix, 1)//'/md.idx', exist=exists)
  if (exists) then
    allpass = .false.
    print *, 'A partly sampled block was written'
  end if

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)

  call global_all(allpass)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (irank == 0) call execute_command_line( &
    'rm -rf '//prefix//'_*.bp '//restart_prefix//'_*.bp')

  call finalise_test(allpass, irank)

contains

  pure real(dp) function expected(c, iter, i, j, k)
    !! Field value for component c at a global vertex and iteration
    integer, intent(in) :: c, iter, i, j, k

    expected = c*1e7_dp + iter*1e4_dp + i*1e2_dp + j*10._dp + k
  end function expected

  subroutine advance(plane_writer, first, last)
    !! Set the fields for each iteration and pass them to the writer
    type(inflow_plane_writer_t), intent(inout) :: plane_writer
    integer, intent(in) :: first, last

    real(dp), allocatable :: data(:, :, :)
    integer :: padded(3), off(3), iter, c, i, j, k

    padded = allocator%get_padded_dims(DIR_C)
    allocate (data(padded(1), padded(2), padded(3)), source=0._dp)
    off = mesh%par%n_offset

    do iter = first, last
      do c = 1, 3
        do k = 1, dims(3)
          do j = 1, dims(2)
            do i = 1, dims(1)
              data(i, j, k) = expected(c, iter, i + off(1), j + off(2), &
                                       k + off(3))
            end do
          end do
        end do
        select case (c)
        case (1)
          call backend%set_field_data(u, data)
        case (2)
          call backend%set_field_data(v, data)
        case (3)
          call backend%set_field_data(w, data)
        end select
      end do
      call plane_writer%sample(iter, iter*dt, u, v, w)
    end do
  end subroutine advance

  subroutine check_block(filename, first_iter)
    !! The block must hold the planes sampled from first_iter onward
    character(*), intent(in) :: filename
    integer, intent(in) :: first_iter

    real(dp), allocatable :: pu(:, :, :), pv(:, :, :), pw(:, :, :), times(:)
    integer :: rdims(3), off(3), s, j, k, iter
    real(dp) :: err, time_err

    call read_inflow_block(filename, read_mesh, pu, pv, pw, times)
    rdims = read_mesh%get_dims(VERT)
    off = read_mesh%par%n_offset

    if (size(pu, 1) /= ntimesteps .or. size(pu, 2) /= rdims(2) &
        .or. size(pu, 3) /= rdims(3)) then
      allpass = .false.
      print *, trim(filename), ': wrong plane array shape', shape(pu)
      return
    end if

    ! Plane values are copied, so they must match exactly; the time stamps
    ! may differ by rounding in iter*dt
    err = 0._dp
    time_err = 0._dp
    do s = 1, ntimesteps
      iter = first_iter + (s - 1)*sample_every
      time_err = max(time_err, abs(times(s) - iter*dt))
      do k = 1, rdims(3)
        do j = 1, rdims(2)
          err = max(err, abs(pu(s, j, k) &
                             - expected(1, iter, x_index, j + off(2), &
                                        k + off(3))))
          err = max(err, abs(pv(s, j, k) &
                             - expected(2, iter, x_index, j + off(2), &
                                        k + off(3))))
          err = max(err, abs(pw(s, j, k) &
                             - expected(3, iter, x_index, j + off(2), &
                                        k + off(3))))
        end do
      end do
    end do

    if (err > 0._dp .or. time_err > 1e-12_dp) then
      allpass = .false.
      print *, 'Rank', irank, trim(filename), ': max errors', err, time_err
    end if
  end subroutine check_block

end program test_inflow_planes
