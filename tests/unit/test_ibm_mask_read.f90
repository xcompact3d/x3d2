program test_ibm_mask_read
  !! Reads an IBM mask through ibm_t and checks it lands on the right points.
  !!
  !! The test writes the mask itself, each rank writing its own slab, with
  !! every point set to a value made from its global (i, j, k) index. A wrong
  !! offset, a swapped axis or a bad copy of the padded block then gives a
  !! wrong value. The file is named after the test executable, so tests that
  !! run at the same time do not share it. The grid has a different size in each direction and none
  !! of them is a multiple of SZ.
  !!
  !!   - no argument: the mask matches the grid and every rank checks its part
  !!     of ep1 against the index pattern;
  !!   - `mismatch`: the mask has one point more in x than the grid, so ibm_t
  !!     must stop the run. ctest looks for the error message in the output.
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, i8, DIR_C, VERT
  use m_field, only: field_t
  use m_ibm, only: ibm_t, iibm_basic
  use m_io_session, only: writer_session_t
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_Barrier
  use m_test_utils, only: initialise_mpi, finalise_test, global_all

  implicit none

  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend
  type(ibm_t) :: ibm
  type(writer_session_t) :: writer
  class(field_t), pointer :: host_ep1
  real(dp), allocatable :: mask(:, :, :)
  integer :: irank, nproc, ierr, i, j, k, nerr
  integer :: dims_global(3), dims(3), offset(3), file_extra(3)
  integer(i8) :: file_dims(3), start_dims(3), count_dims(3)
  character(len=32) :: mode
  character(len=256) :: exe_path, mask_file
  logical :: allpass

  call initialise_mpi(irank, nproc)

  mode = ''
  if (command_argument_count() >= 1) call get_command_argument(1, mode)
  call get_command_argument(0, exe_path)
  mask_file = exe_path(index(exe_path, '/', back=.true.) + 1:)
  mask_file = trim(mask_file)//'.bp'

  ! Split in z only, so each rank has a different z offset
  dims_global = [37, 35, 6*nproc]
  mesh = mesh_t(dims_global, [1, 1, nproc], [2.0_dp, 3.0_dp, 4.0_dp], &
                ['dirichlet', 'dirichlet'], &
                ['periodic ', 'periodic '], &
                ['periodic ', 'periodic '])
  call runtime%init(mesh)
  backend => runtime%backend

  dims = mesh%get_dims(VERT)
  offset = mesh%par%n_offset

  ! The mismatch case writes one point more in x on every rank
  file_extra = 0
  if (trim(mode) == 'mismatch') file_extra(1) = 1

  allocate (mask(dims(1) + file_extra(1), dims(2), dims(3)))
  do k = 1, size(mask, 3)
    do j = 1, size(mask, 2)
      do i = 1, size(mask, 1)
        mask(i, j, k) = expected(i, j, k)
      end do
    end do
  end do

  file_dims = int(mesh%get_global_dims(VERT) + file_extra, i8)
  start_dims = int(offset, i8)
  count_dims = int(shape(mask), i8)

  call writer%open(trim(mask_file), MPI_COMM_WORLD)
  if (irank == 0) call writer%write_data('iibm', int(iibm_basic, i8))
  call writer%write_data('ep1', mask, file_dims, start_dims, count_dims)
  call writer%close()
  deallocate (mask)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ibm = ibm_t(backend, mesh, runtime%host_allocator, trim(mask_file))

  if (trim(mode) == 'mismatch') then
    print *, 'ibm_t read a mask of the wrong shape: &
      &this line must never be reached'
    error stop 1
  end if

  allpass = ibm%iibm == iibm_basic .and. associated(ibm%ep1)
  if (.not. allpass) print *, 'Rank', irank, ': ibm_t did not read ep1'

  if (allpass) then
    host_ep1 => runtime%host_allocator%get_block(DIR_C)
    call backend%get_field_data(host_ep1%data, ibm%ep1)

    nerr = 0
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          if (host_ep1%data(i, j, k) /= expected(i, j, k)) then
            nerr = nerr + 1
            if (nerr <= 5) print *, 'Rank', irank, ': ep1 at', i, j, k, &
              'is', host_ep1%data(i, j, k), 'expected', expected(i, j, k)
          end if
        end do
      end do
    end do
    if (nerr > 0) then
      allpass = .false.
      print *, 'Rank', irank, ':', nerr, 'points of ep1 are wrong'
    end if

    call runtime%host_allocator%release_block(host_ep1)
    call backend%allocator%release_block(ibm%ep1)
  end if

  call global_all(allpass)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (irank == 0) call execute_command_line('rm -rf '//trim(mask_file))

  call finalise_test(allpass, irank)

contains

  pure real(dp) function expected(i, j, k)
    !! The mask value at local point (i, j, k): its global index, packed so
    !! that each axis has its own digits. Exact in single precision too.
    integer, intent(in) :: i, j, k

    expected = real((i + offset(1)) + 100*(j + offset(2)) &
                    + 10000*(k + offset(3)), dp)
  end function expected

end program test_ibm_mask_read
