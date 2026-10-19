Building x3d2
=============

See :ref:`tooling` for details on the tool required to build x3d2.

Start by configuring the build directory, setting ``CMAKE_Fortran_COMPILER`` to
the MPI compiler wrapper to build with (either from the NVIDIA HPC SDK or from
Open MPI). If you give a bare name rather than an absolute path, it must be on
your ``PATH``.

.. code-block:: console

   $ cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Debug

The above is using ``Debug`` build, for release build use ``Release`` instead.

A build with ``-DWITH_MPI=OFF`` takes the compiler itself rather than a wrapper,
so name the one matching the backend (``gfortran``, ``flang``, ``nvfortran``
for ``ENABLE_BACKEND=CUDA``, ``ftn`` on Cray):

.. code-block:: console

   $ cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DWITH_MPI=OFF \
       -DENABLE_BACKEND=CUDA -DBACKEND_ARCH=cc80 -DCMAKE_BUILD_TYPE=Debug

.. note::

   Pass the compiler as a ``-D`` cache variable rather than through the ``FC``
   environment variable. CMake only reads ``FC`` on the first configure of a
   fresh build directory, so ``export FC=...`` has no effect when re-configuring
   an existing one.

Select a GPU backend with ``-DENABLE_BACKEND=CUDA`` or ``-DENABLE_BACKEND=OMP_TGT``;
the default ``OFF`` builds the CPU backend only. Either GPU backend also needs
the target architecture in ``-DBACKEND_ARCH`` (``cc80``, ``gfx942``, ...).
See :doc:`../user/advanced_build` for the full set of options.

Once the build directory is configured, you can build the executable and run the tests as follows:

.. code-block:: bash

   $ cd build
   $ make
   $ make test

Note that ``make test`` is only a launcher for the ``ctest`` executable. By default, ``ctest`` does not show the output of test executables on failure. If one or more tests fail, you probably want to run the tests with:

.. code-block:: bash

   $ ctest --output-on-failure

instead of ``make test``.

The main executable will be built in the ``build/bin`` directory as ``xcompact``.

Building the Standalone Poisson Solver
--------------------------------------

``x3d2_poisson_solver`` is a standalone driver for the FFT-based Poisson solver.
The build of this target links only the Poisson library ``x3d2_poisson``, the backend library and the driver.
It does not compile the full ``x3d2`` solver.

.. code-block:: bash

   $ cmake --build build --target x3d2_poisson_solver

The binary is ``build/bin/x3d2_poisson_solver``. Its usage line is:

.. code-block:: console

   x3d2_poisson_solver nx ny nz [bc_x bc_y bc_z] [--repeat N]

``nx``, ``ny`` and ``nz`` are the global grid dimensions.
Each of ``bc_x`` and ``bc_y`` is ``periodic`` (the default) or ``dirichlet``.
The argument ``bc_z`` must be ``periodic``.
A direction with ``dirichlet`` needs an odd grid size.
Here ``dirichlet`` is the velocity condition. The pressure gets homogeneous Neumann conditions in that direction.
``--repeat N`` sets the number of solves (default 1).

The driver supports four boundary condition cases.
Each case has a three digit name.
A digit is ``1`` if the direction is ``dirichlet`` and ``0`` if it is ``periodic``.
The cases are ``000`` (all periodic), ``100`` (``bc_x`` is ``dirichlet``), ``010`` (``bc_y`` is ``dirichlet``) and ``110`` (``bc_x`` and ``bc_y`` are ``dirichlet``).
The OpenMP and CUDA backends solve all four cases.
The ``010`` case runs on one MPI rank only.
On the CUDA backend, the ``100`` case on more than one rank needs the number of x cells to divide by the number of ranks.
On the CUDA backend, the ``110`` case on more than one rank needs the same for the number of x cells and the number of y cells.

The backend follows ``ENABLE_BACKEND``.
The ``OMP_TGT`` backend has no FFT Poisson solver, so the driver stops with that backend.
The OpenMP backend needs an MPI build.
The option ``WITH_MPI=OFF`` switches ``2decomp`` off.
The OpenMP backend then has no FFT Poisson solver.

The default right-hand side is ``cos(2 pi x)``.
To use another right-hand side, edit ``fill_rhs`` in ``src/poisson_solver.f90`` and rebuild.
The right-hand side must have a zero mean.
The driver prints a warning if the mean is not zero.

The driver reports the grid size, the boundary conditions, the number of ranks, the backend, and the minimum and mean solve time over the repeats.
It does not check the result.
The test ``tests/verification/test_poisson_bc.f90`` checks the solve against an analytical solution.
For example:

.. code-block:: console

   $ mpirun -np 2 build/bin/x3d2_poisson_solver 128 64 32 --repeat 5
