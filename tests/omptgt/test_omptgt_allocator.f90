program test_allocator_omptgt
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, pi, DIR_X
  use m_mesh, only: mesh_t

  use m_omptg_allocator, only: omptgt_allocator_t

  implicit none

  logical :: allpass
  integer, dimension(3) :: dims, nproc_dir
  real(dp) :: L_global(3)
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  class(allocator_t), allocatable :: allocator
  class(mesh_t), allocatable :: mesh
  class(field_t), pointer :: ptr1, ptr2, ptr3
  integer, allocatable :: l(:)
  integer :: ierr

  call MPI_Init(ierr)

  dims = [8, 8, 8]
  nproc_dir = [1, 1, 1]
  L_global = [2*pi, 2*pi, 2*pi]

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims, nproc_dir, L_global, BC_x, BC_y, BC_z)

  allocator = omptgt_allocator_t(mesh, 8)

  allpass = .true.

  ! Get the list of ids for free blocks.  By default there are none
  ! and returned list is [0].
  l = allocator%get_block_ids()
  if (size(l) /= 1 .or. l(1) /= 0) then
    allpass = .false.
    write (stderr, '(a)') 'Free list is initialised empty... failed'
  else
    write (stderr, '(a)') 'Free list is initialised empty... passed'
  end if

  ! Request two blocks and release them in reverse order.  List should
  ! contain two free blocks. (1 -> 2)
  ptr1 => allocator%get_block(DIR_X)
  ptr2 => allocator%get_block(DIR_X)
  call allocator%release_block(ptr2)
  call allocator%release_block(ptr1)

  if (.not. all(allocator%get_block_ids() == [1, 2])) then
    allpass = .false.
    write (stderr, '(a)') 'Blocks are released correctly... failed'
  else
    write (stderr, '(a)') 'Blocks are released correctly... passed'
  end if

  !! Destroy the free list and check that the list is empty again.
  call allocator%destroy()
  l = allocator%get_block_ids()
  if (size(l) /= 1 .or. l(1) /= 0 .or. allocator%next_id /= 0) then
    allpass = .false.
    write (stderr, '(a)') 'Free list is correctly destroyed... failed'
  else
    write (stderr, '(a)') 'Free list is correctly destroyed... passed'
  end if

  ! Request a block from a list of three.  This should grab the first
  ! block on top of the pile and reduce the free list to two blocks.
  ptr1 => allocator%get_block(DIR_X)
  ptr2 => allocator%get_block(DIR_X)
  ptr3 => allocator%get_block(DIR_X)
  call allocator%release_block(ptr3)
  call allocator%release_block(ptr2)
  call allocator%release_block(ptr1)
  ptr1 => allocator%get_block(DIR_X)

  if (.not. all(allocator%get_block_ids() == [2, 3])) then
    allpass = .false.
    write (stderr, '(a)') 'Block is correctly allocated... failed'
  else
    write (stderr, '(a)') 'Block is correctly allocated... passed'
  end if

  ! TODO: Check that data is resident on device
  ! TODO: Check that data is freed from device

  call allocator%destroy()

  call MPI_Finalize(ierr)
end program test_allocator_omptgt
