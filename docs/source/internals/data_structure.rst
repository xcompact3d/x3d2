Data Structure
==============

x3d2 uses a specialised data structure internally in the backends for
improving the performance of compact scheme operators. Data structure
implements data layout optimisation based on the spatial direction of
the compact operators. We define the following options in `common.f90`
`DIR_X`, `DIR_Y`, `DIR_Z`, `DIR_C`.

`DIR_C` is a typical Cartesian ordering with x-dir first, and the
remaining are the orderings used internally in the backends. The shape
is in the form of `(SZ, n_x, n_y_padded*n_z_padded/SZ)`, where the first
dim is backend-specific constant SZ, the second is the number of points
in `DIR_?`, and the last dim is the number of groups in the domain. This
specific data layout helps with vectorisation on CPUs and thread level
parallelism on GPUs. As a consequence of this specialised layout, we
have a particular loop pattern in various kernels.

.. code:: fortran

  do k = 1, n_group
    do j = 1, n
      do i = 1, SZ
        f1(i, j, k) = f2(i, j, k)
      end do
    end do
  end do

The `k` loop is typically distributed along the available cores in a
CPU or along the thread blocks in a GPU. Then each core or thread block
carries out the explicit `j` loop while processing SZ many lines
concurrently.

Note: in GPU kernels the 'k' loop is managed by the blocks and 'i' loop
is managed by the threads. Therefore, a GPU kernel only contains the 'j'
loop.

.. code:: cuda

  ! kernel call from host side
  call gpu_kernel<<<dim3(n_group, 1, 1), dim3(SZ, 1, 1)>>>(...)

  ! the kernel implementation
  __global__ gpu_kernel(...)
    k = blockIdx.x
    i = threadIdx.x
    do i = 1, SZ
      f1(i, j, k) = f2(i, j, k)
    end do

