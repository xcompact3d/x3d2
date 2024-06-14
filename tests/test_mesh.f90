program test_allocator
  use mpi
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_common, only: DIR_X, pi, dp, CELL, DIR_C, DIR_Y, DIR_Z, VERT, X_FACE

  implicit none

  logical :: allpass
  integer, dimension(3) :: nproc_dir, dims_global
  logical, dimension(3) :: periodic_BC
  real(dp), dimension(3) :: L_global 
  class(allocator_t), allocatable :: allocator
  class(mesh_t), allocatable :: mesh
  class(field_t), pointer :: ptr1, ptr2, ptr3
  integer, dimension(3) :: dims
  integer, dimension(3) :: dims_check
  integer :: ierr
  integer :: n_cell, n_vert, n_x_face
  real(dp), parameter :: eps = 1e-8

  call MPI_Init(ierr)
  allpass = .true.

  ! Global number of cells in each direction
  dims_global = [16, 4, 4]

  ! Global domain dimensions
  L_global = [1._dp, 1._dp, 1._dp]

  ! Domain decomposition in each direction
  nproc_dir = [4, 1, 1]

  periodic_BC = [.false. , .true., .false.]

  mesh = mesh_t(dims_global, nproc_dir, L_global, periodic_BC)

  if (mesh%par%nproc /= 4) then
    allpass = .false.
    print *, "To be run on 4 ranks"
  end if

  allocator = allocator_t(mesh, 8)

  ptr1 => allocator%get_block(DIR_X, CELL)
  ptr2 => allocator%get_block(DIR_X, VERT)
  ptr3 => allocator%get_block(DIR_Z, X_FACE)

  n_cell = mesh%get_n(ptr1)
  n_vert = mesh%get_n(ptr2)

  ! if last rank
  if (mesh%par%nrank == 3) then
    if (.not. (n_cell == 3 .and. n_vert == 4)) then
      allpass = .false.
      print *, "error in get_n and last rank, n_cell=", n_cell, "n_vert=", n_vert
    end if
  else
    if (.not. (n_cell == 4 .and. n_vert == 4)) then
      allpass = .false.
      print *, "error in get_n, n_cell=", n_cell, "n_vert=", n_vert
    end if
  end if

  n_x_face = mesh%get_n(ptr3)
  if (.not. n_x_face == 3) then
      allpass = .false.
      print *, "error in get_n for x_face, n_x_face=", n_x_face
  end if

  dims = mesh%get_padded_dims(DIR_C)
  dims_check = [8, 8, 4] ! No padding in Z

  if (.not. all(dims(:) == dims_check(:))) then
    allpass = .false.
    print *, "error with padded dimensions, dims_padded=", dims
  end if

  if (.not. abs(mesh%geo%d(DIR_X) - 1._dp/(16-1)) .lt. eps ) then
    allpass = .false.
    print *, "error with geo%d, non periodic BC"
  end if

  if (.not. abs(mesh%geo%d(DIR_Y) - 1._dp/4) .lt. eps ) then
    allpass = .false.
    print *, "error with geo%d, periodic BC"
  end if

  call allocator%destroy()

  if (allpass) then
    write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)
end program test_allocator
