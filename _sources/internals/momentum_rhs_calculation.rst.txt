Computing the spatial contribution to the momentum equation
===========================================================

The first stage of solving the momentum equation at a given time-step
:math:`k` is to evaluate the spatial contribution

.. math::
   :label: rhs

   \mathbf{F}^k = - \frac{1}{2}[\nabla (\mathbf{u}^k \otimes \mathbf{u}^k)
   + (\mathbf{u}^k \cdot \nabla) \mathbf{u}^k] + \nu \nabla^2\mathbf{u}^k

This term involves spatial derivatives along the three directions
:math:`\hat{\mathbf{x}}`, :math:`\hat{\mathbf{y}}` and
:math:`\hat{\mathbf{z}}`.  In order to reduce the number of
communications between MPI processes, the computation is performed in
three stages corresponding to three different data configurations:

- In the x-configuration, the process has access to a data slab
  covering the x-y plane.  This a 1D domain decomposition along the
  :math:`\hat{\mathbf{z}}` direction.  Derivatives along
  :math:`\hat{\mathbf{x}}` can be computed without communications.
- In the z-configuration, a process has access to a data slab
  covering the y-z plane.  This a 1D domain decomposition along the
  :math:`\hat{\mathbf{x}}` direction.  Derivatives along
  :math:`\hat{\mathbf{z}}` can be computed without communications.
- In the y-configuration, a process has access to a data slab
  covering the x-y plane.  This a 1D domain decomposition along the
  :math:`\hat{\mathbf{z}}` direction.  Derivatives along
  :math:`\hat{\mathbf{y}}` can be computed without communications.

In the i-configuration, field data (e.g. the x-component of the
velocity field) is stored in a three-dimensional array of dimension
``dimension(SZ, N, M)`` where ``N`` is the size of the domain in the i
direction.  This data layout decomposes slab data into ``M * SZ`` ` 1D
arrays (called pencils) of size ``N``.  Pencils are processed in
batches of size ``SZ`` allowing for vectorised operations.

The algorithm for computing :math:numref:`rhs` then becomes:

.. code:: fortran

   ! Instances of type slab_t have access to slab data
   ! and can compute derivatives and contributions to
   ! momentum equation along a specific direction.
   type(slab_t) :: xslab, yslab, zslab

   ! ...

   ! Perform necesseary communications across processes
   ! asynchonously
   call async_rotate_x_to_z(xslab, z_comm_buffer)

   ! Meanwhile, compute the contribution in the x direction
   x_contrib = x_slab%transport()

   ! Assuming x to z transfer is done, arrange field data
   ! in the correct z-configuration layout and compute
   ! z contribution
   call from_comms_buffer(z_comm_buffer, z_slab)
   z_dir_contrib = z_slab%transport()

   ! Transfer the result back
   call async_rotate_z_to_x(z_dir_contrib, x_slab_z)

   ! Meanwhile, switch from x-config to y-config
   ! (no communications needed, only data reshuffling)
   call rotate_x_to_y(x_slab, y_slab)
   y_dir_contrib = yslab%transport()

   F = x_dir_contrib%data + y_dir_contrib%data + z_dir_contrib%data

The ``slab3d_t`` type
---------------------

Computation of :math:numref:`rhs` is supported by instances of the
``slab3d_t`` derived type.  These hold data for the three components
of a 3D vector in 3D space (thinks :math:`u_x`, :math:`u_y` and
:math:`u_z`) within a single data slab in the x, y or z-

The purpose of ``slab3d_t`` instances is to compute spatial
derivatives along a given direction referred to as the slab's primary
direction.  More precisely, user code is given access to type-bound
procedures such as ``transport`` or ``divergence`` which compute
contributions to physical quantities coming from spatial derivatives
in the slab's primary direction.

For more details, see the API docs.
