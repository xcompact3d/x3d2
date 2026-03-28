program test_mesh
  use mpi
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_common, only: DIR_X, pi, dp, CELL, DIR_C, DIR_Y, DIR_Z, VERT, Z_FACE

  implicit none

  logical :: allpass
  integer, dimension(3) :: nproc_dir, dims_global
  real(dp), dimension(3) :: L_global
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  class(allocator_t), allocatable :: allocator
  class(mesh_t), allocatable :: mesh
  class(field_t), pointer :: ptr1, ptr2, ptr3
  integer, dimension(3) :: dims
  integer, dimension(3) :: dims_check
  integer :: ierr, nrank, nproc
  integer :: n_cell, n_vert, n_x_face
  real(dp), parameter :: eps = 1e-8

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  allpass = .true.

  if (nproc /= 4) then
    allpass = .false.
    print *, "To be run on 4 ranks"
  end if

  if (nrank == 0) then
    print *, "Generic decomposition"
  end if
  call run_test_mesh(.false., allpass)

#ifdef WITH_2DECOMPFFT
  if (nrank == 0) then
    print *, "2decomp decomposition"
  end if
  call run_test_mesh(.true., allpass)
#endif

  if (allpass) then
    write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

contains

  subroutine run_test_mesh(use_2decomp, allpass)
    logical, intent(in) :: use_2decomp
    logical, intent(inout) :: allpass
    integer, dimension(3) :: nproc_dir, dims_global, dims
    real(dp), dimension(3) :: L_global
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
    integer, dimension(4) :: n_vert_z

    ! Global number of cells in each direction
    dims_global = [4, 4, 16]

    ! Global domain dimensions
    L_global = [1._dp, 1._dp, 1._dp]

    ! Domain decomposition in each direction
    nproc_dir = [1, 1, 4]

    BC_x = [character(len=20) :: 'neumann', 'dirichlet']
    BC_y = [character(len=20) :: 'periodic', 'periodic']
    BC_z = [character(len=20) :: 'dirichlet', 'neumann']

    mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, &
                  use_2decomp=use_2decomp)

    ! Expected decomposition by 2decomp and generic
    if (use_2decomp) then
      n_vert_z = [3, 4, 4, 5]
    else
      n_vert_z = [4, 4, 4, 4]
    end if

    allocator = allocator_t(mesh%get_dims(VERT), 8)

    ptr1 => allocator%get_block(DIR_Z, CELL)
    ptr2 => allocator%get_block(DIR_Z, VERT)
    ptr3 => allocator%get_block(DIR_X, Z_FACE)

    n_cell = mesh%get_n(ptr1)
    n_vert = mesh%get_n(ptr2)

    ! if last rank
    if (mesh%par%nrank == 3) then
      if (.not. (n_cell == n_vert_z(mesh%par%nrank + 1) - 1 &
                 .and. n_vert == n_vert_z(mesh%par%nrank + 1))) then
        allpass = .false.
        print *, mesh%par%nrank, "error in get_n and last rank, n_cell=", &
          n_cell, "n_vert=", n_vert
      end if
    else
      if (.not. (n_cell == n_vert_z(mesh%par%nrank + 1) &
                 .and. n_vert == n_vert_z(mesh%par%nrank + 1))) then
        allpass = .false.
        print *, mesh%par%nrank, "error in get_n, n_cell=", &
          n_cell, "n_vert=", n_vert
      end if
    end if

    n_x_face = mesh%get_n(ptr3)
    if (.not. n_x_face == 3) then
      allpass = .false.
      print *, "error in get_n for x_face, n_x_face=", n_x_face
    end if

    dims = allocator%get_padded_dims(DIR_C)
    dims_check = [8, 8, n_vert_z(mesh%par%nrank + 1)] ! No padding in Z

    if (.not. all(dims(:) == dims_check(:))) then
      allpass = .false.
      print *, "error with padded dimensions, dims_padded=", dims
    end if

    if (.not. abs(mesh%geo%d(DIR_Z) - 1._dp/(16 - 1)) < eps) then
      allpass = .false.
      print *, "error with geo%d, non periodic BC"
    end if

    if (.not. abs(mesh%geo%d(DIR_Y) - 1._dp/4) < eps) then
      allpass = .false.
      print *, "error with geo%d, periodic BC"
    end if

    call allocator%destroy()

  end subroutine

end program test_mesh
