  program test_allocator_cuda
    use iso_fortran_env, only: stderr => error_unit

    use m_allocator, only: allocator_t, field_t
    use m_common, only: DIR_X
    use m_cuda_allocator, only: cuda_allocator_t

    implicit none

    logical :: allpass
    integer, parameter :: dims(3) = [8, 8, 8]
    class(allocator_t), allocatable :: allocator
    class(field_t), pointer :: ptr1, ptr2, ptr3
    integer, allocatable :: l(:)

    allocator = cuda_allocator_t(dims(1), dims(2), dims(3), 8)

    allpass = .true.

    ! Get the list of ids for free blocks.  By default there are none
    ! and returned list is [0].
    l = allocator%get_block_ids()
    if (size(l) /= 1 .or. l(1) /= 0) then
      allpass = .false.
      write(stderr, '(a)') 'Free list is initialised empty... failed'
    else
      write(stderr, '(a)') 'Free list is initialised empty... passed'
    end if

    ! Request two blocks and release them in reverse order.  List should
    ! contain two free blocks. (1 -> 2)
    ptr1 => allocator%get_block(DIR_X)
    ptr2 => allocator%get_block(DIR_X)
    call allocator%release_block(ptr2)
    call allocator%release_block(ptr1)

    if (.not. all(allocator%get_block_ids() .eq. [1, 2])) then
      allpass = .false.
      write(stderr, '(a)') 'Blocks are released correctly... failed'
    else
      write(stderr, '(a)') 'Blocks are released correctly... passed'
    end if

    !! Destroy the free list and check that the list is empty again.
    call allocator%destroy()
    l = allocator%get_block_ids()
    if (size(l) /= 1 .or. l(1) /= 0 .or. allocator%next_id /=0) then
      allpass = .false.
      write(stderr, '(a)') 'Free list is correctly destroyed... failed'
    else
      write(stderr, '(a)') 'Free list is correctly destroyed... passed'
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

    if (.not. all(allocator%get_block_ids() .eq. [2, 3])) then
      allpass = .false.
      write(stderr, '(a)') 'Block is correctly allocated... failed'
    else
      write(stderr, '(a)') 'Block is correctly allocated... passed'
    end if

    call allocator%destroy()
  end program test_allocator_cuda
