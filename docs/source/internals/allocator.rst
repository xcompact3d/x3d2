Memory management
=================

The parallel processing carried out in x3d2 requires that the
computational domain is sliced into slabs along a specific direction
of space.

Data slabs (instances of the derived type `slab_t`) hold a pointer to
an `allocator_t` instance that is responsible for the allocation and
dispatching of memory blocks to data slabs. See `(API)allocator_t
<https://xcompact3d.github.io/x3d2/api/type/allocator_t.html>`_. Whenever
data slabs are no longer required, the associated memory blocks can be
released to the allocator, ready to be dispatched to other slabs or
objects who might need them.

An instance of `allocator_t` maintains a linked list of memory blocks,
which each element in the list labelled with a unique identifier:

.. graphviz::

   digraph freelist {
       rank=LR;
       memblock1 [label="memory block\lid: 1"];
       memblock2 [label="memory block\lid: 2"];
       memblock3 [label="memory block\lid: 3"];

       memblock1 -> memblock2 -> memblock3;
  }

Each memory block is an instance of the `memblock_t` derived type
(`(API)memblock_t
<https://xcompact3d.github.io/x3d2/api/type/memblock_t.html>`_).  In
addition to its integer identifier, a instance of the `memblock_t`
type holds a 3D data array and a pointer to the next block in the
list.

All memory blocks in an allocator's list have the same size.  An
allocator can be created by passing in the dimension of the memory
blocks:

.. code:: fortran

	  allocator = allocator_t([32 16 32])

Given an `allocator_t` instance, memory blocks can requested using the
`get_block` function.

.. code:: fortran

    use m_allocator, only: memory_block_t, allocator_t
    ! ...
    type(memory_block_t), pointer :: memblock_ptr
    ! ..
    ptr => allocator%get_block()

Note that `get_block` returns a *pointer* to a memory block instead of
a copy.  In practice, the head of the block list is detached from the
list (pop operation):

.. graphviz::

    digraph allocation {

    memblock1 [label="memory block\lid: 1"];
    memblock2 [label="memory block\lid: 2"];
    memblock3 [label="memory block\lid: 3"];
    ptr [shape=box, label="memblock_ptr"]

    memblock2 -> memblock3;
    ptr -> memblock1;
   }

Dispatched memory blocks ban be released to the allocator by calling
the `release_block` subroutine:

.. code:: fortran

   ! Request two blocks and release the first one, which becomes the
   ! head of the free block list
   ptr1 => allocator%get_block()
   ptr2 => allocator%get_block()
   call allocator%release_block(ptr1)

Released blocks are pushed on top of the allocator's free blocks list:

.. graphviz::

    digraph release {

    memblock1 [label="memory block\lid: 1"];
    memblock2 [label="memory block\lid: 2"];
    memblock3 [label="memory block\lid: 3"];
    ptr1 [shape=box]
    ptr2 [shape=box]
    null [shape=box]

    memblock1 -> memblock3;
    ptr2 -> memblock2;
    ptr1 -> null
   }

Memory block allocation
-----------------------

If a request is made to an allocator whose free block list is empty, a
new block is allocated before it is immediately dispacted to the
requesting object.  As a consequence the memory footprint of a program
is expected to grow in early stages of a program's lifetime, before
reaching maximum where enough memory has been allocated and some can
be reused.
