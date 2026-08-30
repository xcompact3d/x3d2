program test_boundary_planes
  !! Face stamping must land on the intended y-planes and leave every other
  !! plane untouched.
  !!
  !! DIR_X fields are (SZ, nx, n_groups) with SZ consecutive y-vertices in the
  !! first index, but the two backends order the groups differently:
  !!
  !!   OMP:  dir_k = n_y_blocks*(z - 1) + y_block   (y-block fast)
  !!   CUDA: b     = nz*(y_block - 1) + z           (z fast)
  !!
  !! A routine written against the wrong convention still touches plausible
  !! groups, so it fails silently by scribbling on interior planes. The two
  !! orderings coincide when ny <= SZ, so ny here is deliberately larger than
  !! SZ on both backends (SZ is 32 on CUDA, 16 on OMP) and is not a multiple
  !! of either, which also exercises the partially-filled last y-block.
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_C, VERT, Y_FACE
  use m_field, only: field_t
  use m_mesh, only: mesh_t

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
  class(field_t), pointer :: f, f_start

  ! ny is not a multiple of either backend's SZ, so the last y-block is
  ! partially filled and the group orderings do not coincide.
  integer, parameter :: dims_global(3) = [16, 65, 8]
  integer, parameter :: nproc_dir(3) = [1, 1, 1]
  real(dp), parameter :: lengths(3) = [1._dp, 1._dp, 1._dp]
  character(len=9), parameter :: bc_periodic(2) = ['periodic ', 'periodic ']
  character(len=9), parameter :: bc_slip(2) = ['neumann  ', 'neumann  ']

  real(dp), parameter :: sentinel = 7._dp, bottom = -1._dp, top = -2._dp

  real(dp), allocatable :: data(:, :, :), expected(:, :, :)
  real(dp) :: tolerance
  integer :: dims(3), dims_padded(3), i, j, k, ierr
  logical :: all_pass

  call MPI_Init(ierr)
  all_pass = .true.
  tolerance = 2000._dp*epsilon(1._dp)

  mesh = mesh_t(dims_global, nproc_dir, lengths, &
                bc_periodic, bc_slip, bc_periodic)
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

  if (dims(2) <= SZ) error stop &
    'test_boundary_planes needs ny > SZ to distinguish the group orderings.'

  dims_padded = allocator%get_padded_dims(DIR_C)
  allocate (data(dims_padded(1), dims_padded(2), dims_padded(3)))
  allocate (expected(dims_padded(1), dims_padded(2), dims_padded(3)))

  f => allocator%get_block(DIR_X, VERT)
  f_start => allocator%get_block(DIR_X, VERT)

  ! --- field_set_face: scalar values on the two y-faces ---------------------
  call f%fill(sentinel)
  call backend%field_set_face(f, bottom, top, Y_FACE)
  call backend%get_field_data(data, f)

  expected = sentinel
  expected(:, 1, :) = bottom
  expected(:, dims(2), :) = top
  call check_planes('field_set_face', data, expected, dims, tolerance, &
                    all_pass)

  ! --- field_set_face_from_field: plane values copied from another field ----
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        ! Vary within each plane so a mis-indexed copy cannot coincide.
        data(i, j, k) = 100._dp*real(i, dp) + real(k, dp) + 0.5_dp*real(j, dp)
      end do
    end do
  end do
  call backend%set_field_data(f_start, data)

  expected = sentinel
  expected(:, 1, :) = data(:, 1, :)
  expected(:, dims(2), :) = data(:, dims(2), :)

  call f%fill(sentinel)
  call backend%field_set_face_from_field(f, f_start, 0._dp, Y_FACE)
  call backend%get_field_data(data, f)
  call check_planes('field_set_face_from_field', data, expected, dims, &
                    tolerance, all_pass)

  call allocator%release_block(f)
  call allocator%release_block(f_start)
  deallocate (data, expected)

  if (.not. all_pass) error stop 'FAIL'
  print *, 'PASS'
  call MPI_Finalize(ierr)

contains

  subroutine check_planes(label, actual, want, extents, tol, pass)
    !! Compare every y-plane separately so a failure names the stray plane.
    character(*), intent(in) :: label
    real(dp), intent(in) :: actual(:, :, :), want(:, :, :)
    integer, intent(in) :: extents(3)
    real(dp), intent(in) :: tol
    logical, intent(inout) :: pass

    real(dp) :: plane_error, worst_error
    integer :: plane, worst_plane

    worst_error = 0._dp
    worst_plane = 0
    do plane = 1, extents(2)
      plane_error = maxval(abs(actual(1:extents(1), plane, 1:extents(3)) - &
                               want(1:extents(1), plane, 1:extents(3))))
      if (plane_error > worst_error) then
        worst_error = plane_error
        worst_plane = plane
      end if
    end do

    if (worst_error > tol) then
      print *, 'FAIL: ', label, ' worst y-plane=', worst_plane, &
        ' error=', worst_error
      pass = .false.
    else
      print *, 'PASS: ', label
    end if
  end subroutine check_planes

end program test_boundary_planes
