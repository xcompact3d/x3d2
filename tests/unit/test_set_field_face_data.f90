program test_set_field_face_data
  !! set_field_face_data must write the given host planes into the face
  !! planes of a DIR_X field and leave every other point unchanged.

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t
  use m_common, only: dp, CELL, DIR_C, DIR_X, VERT, X_FACE, Y_FACE
  use m_mesh, only: mesh_t
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  type(backend_runtime_t), target :: runtime
  class(base_backend_t), pointer :: backend
  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime_c
  class(base_backend_t), pointer :: backend_c
  type(mesh_t), target :: mesh_c

  class(field_t), pointer :: fld, fld_c
  real(dp), dimension(:, :, :), allocatable :: ref, arr
  real(dp), dimension(:, :), allocatable :: p_start, p_end
  integer :: dims(3), shape_c(3), i, j, k
  integer :: dims_c(3)
  integer :: nrank, nproc
  logical :: allpass

  call initialise_mpi(nrank, nproc)

  ! ny = 37 is not a multiple of the pencil size, so the padded last
  ! y-block is exercised.
  mesh = mesh_t([16, 37, 24], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"], &
                ["periodic", "periodic"])

  call runtime%init(mesh)
  backend => runtime%backend

  dims = mesh%get_dims(VERT)
  ! The host arrays take the shape of a DIR_C block, which a device
  ! backend may pad. Only the first dims points are compared.
  fld_c => backend%allocator%get_block(DIR_C, VERT)
  shape_c = fld_c%get_shape()
  call backend%allocator%release_block(fld_c)
  allocate (ref(shape_c(1), shape_c(2), shape_c(3)))
  allocate (arr(shape_c(1), shape_c(2), shape_c(3)))
  ref = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        ref(i, j, k) = real(i + 100*j + 10000*k, dp)
      end do
    end do
  end do

  fld => backend%allocator%get_block(DIR_X, VERT)
  allpass = .true.

  ! Y_FACE with both planes
  call backend%set_field_data(fld, ref)
  allocate (p_start(dims(1), dims(3)), p_end(dims(1), dims(3)))
  do k = 1, dims(3)
    do i = 1, dims(1)
      p_start(i, k) = -real(i + 1000*k, dp)
      p_end(i, k) = -real(i + 1000*k, dp) - 0.5_dp
    end do
  end do
  call backend%set_field_face_data(fld, p_start, p_end, Y_FACE)
  arr = 0._dp
  call backend%get_field_data(arr, fld)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        if (j == 1) then
          if (arr(i, j, k) /= p_start(i, k)) allpass = .false.
        else if (j == dims(2)) then
          if (arr(i, j, k) /= p_end(i, k)) allpass = .false.
        else
          if (arr(i, j, k) /= ref(i, j, k)) allpass = .false.
        end if
      end do
    end do
  end do
  if (.not. allpass) print *, "Y_FACE with both planes failed"

  ! Y_FACE without data_end must leave the top plane alone
  call backend%set_field_data(fld, ref)
  call backend%set_field_face_data(fld, p_start, face=Y_FACE)
  arr = 0._dp
  call backend%get_field_data(arr, fld)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        if (j == 1) then
          if (arr(i, j, k) /= p_start(i, k)) allpass = .false.
        else
          if (arr(i, j, k) /= ref(i, j, k)) allpass = .false.
        end if
      end do
    end do
  end do
  if (.not. allpass) print *, "Y_FACE without data_end failed"
  deallocate (p_start, p_end)

  ! X_FACE with both planes
  call backend%set_field_data(fld, ref)
  allocate (p_start(dims(2), dims(3)), p_end(dims(2), dims(3)))
  do k = 1, dims(3)
    do j = 1, dims(2)
      p_start(j, k) = -real(j + 1000*k, dp)
      p_end(j, k) = -real(j + 1000*k, dp) - 0.5_dp
    end do
  end do
  call backend%set_field_face_data(fld, p_start, p_end, X_FACE)
  arr = 0._dp
  call backend%get_field_data(arr, fld)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        if (i == 1) then
          if (arr(i, j, k) /= p_start(j, k)) allpass = .false.
        else if (i == dims(1)) then
          if (arr(i, j, k) /= p_end(j, k)) allpass = .false.
        else
          if (arr(i, j, k) /= ref(i, j, k)) allpass = .false.
        end if
      end do
    end do
  end do
  if (.not. allpass) print *, "X_FACE with both planes failed"

  ! X_FACE without data_end must leave the outlet plane alone
  call backend%set_field_data(fld, ref)
  call backend%set_field_face_data(fld, p_start, face=X_FACE)
  arr = 0._dp
  call backend%get_field_data(arr, fld)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        if (i == 1) then
          if (arr(i, j, k) /= p_start(j, k)) allpass = .false.
        else
          if (arr(i, j, k) /= ref(i, j, k)) allpass = .false.
        end if
      end do
    end do
  end do
  if (.not. allpass) print *, "X_FACE without data_end failed"

  deallocate (p_start, p_end, ref, arr)
  call backend%allocator%release_block(fld)

  ! A CELL field on a non-periodic y. ny = 129 vertices (128 is a multiple
  ! of every SZ) gives one y-block more for the vertices than for the cells.
  ! The DIR_X group stride follows the vertices, so the planes of every
  ! z > 1 need the padded stride, not the CELL y-block count.
  mesh_c = mesh_t([16, 129, 24], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
                  ["periodic ", "periodic "], &
                  ["neumann  ", "neumann  "], &
                  ["periodic ", "periodic "])
  call runtime_c%init(mesh_c)
  backend_c => runtime_c%backend

  dims_c = mesh_c%get_dims(CELL)
  fld_c => backend_c%allocator%get_block(DIR_C, CELL)
  shape_c = fld_c%get_shape()
  call backend_c%allocator%release_block(fld_c)
  allocate (ref(shape_c(1), shape_c(2), shape_c(3)))
  allocate (arr(shape_c(1), shape_c(2), shape_c(3)))
  ref = 0._dp
  do k = 1, dims_c(3)
    do j = 1, dims_c(2)
      do i = 1, dims_c(1)
        ref(i, j, k) = real(i + 100*j + 10000*k, dp)
      end do
    end do
  end do

  fld => backend_c%allocator%get_block(DIR_X, CELL)
  call backend_c%set_field_data(fld, ref)
  allocate (p_start(dims_c(1), dims_c(3)), p_end(dims_c(1), dims_c(3)))
  do k = 1, dims_c(3)
    do i = 1, dims_c(1)
      p_start(i, k) = -real(i + 1000*k, dp)
      p_end(i, k) = -real(i + 1000*k, dp) - 0.5_dp
    end do
  end do
  call backend_c%set_field_face_data(fld, p_start, p_end, Y_FACE)
  arr = 0._dp
  call backend_c%get_field_data(arr, fld)
  do k = 1, dims_c(3)
    do j = 1, dims_c(2)
      do i = 1, dims_c(1)
        if (j == 1) then
          if (arr(i, j, k) /= p_start(i, k)) allpass = .false.
        else if (j == dims_c(2)) then
          if (arr(i, j, k) /= p_end(i, k)) allpass = .false.
        else
          if (arr(i, j, k) /= ref(i, j, k)) allpass = .false.
        end if
      end do
    end do
  end do
  if (.not. allpass) print *, "CELL Y_FACE on a non-periodic y failed"

  deallocate (p_start, p_end, ref, arr)
  call backend_c%allocator%release_block(fld)

  call finalise_test(allpass, nrank)

end program test_set_field_face_data
