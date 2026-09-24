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
       -DENABLE_BACKEND=CUDA -DCMAKE_BUILD_TYPE=Debug

.. note::

   Pass the compiler as a ``-D`` cache variable rather than through the ``FC``
   environment variable. CMake only reads ``FC`` on the first configure of a
   fresh build directory, so ``export FC=...`` has no effect when re-configuring
   an existing one.

Select a GPU backend with ``-DENABLE_BACKEND=CUDA`` or ``-DENABLE_BACKEND=OMP_TGT``;
the default ``OFF`` builds the CPU backend only. See :doc:`../user/advanced_build`
for the full set of options.

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
---------------------------------------

``x3d2_poisson_solver`` is a standalone driver for the FFT-based Poisson
solver. Building it pulls in only the Poisson core library
(``x3d2_poisson``), the backend it runs against and the driver itself,
not the full x3d2 solver:

.. code-block:: bash

   $ cmake --build build --target x3d2_poisson_solver

The resulting binary is written to ``build/bin/x3d2_poisson_solver`` and
takes the following arguments:

.. code-block:: console

   x3d2_poisson_solver nx ny nz [bc_x bc_y bc_z] [--repeat N]

``nx``, ``ny`` and ``nz`` are the global grid dimensions. Each of
``bc_x``, ``bc_y`` and ``bc_z`` is ``periodic`` (the default) or
``dirichlet``; a direction using ``dirichlet`` requires an odd grid size.
``--repeat N`` sets how many times the solve is repeated (default 1). The
backend used (OpenMP, CUDA or OpenMP target offload) follows whichever
``ENABLE_BACKEND`` the build was configured with. The OpenMP backend
currently supports only the ``000`` and ``100`` boundary condition cases.
The ``010`` and ``110`` cases are not yet supported on OpenMP (the CUDA
backend supports all four).

The driver solves a manufactured cosine Poisson problem and reports the
grid size, boundary conditions, rank count, backend, and the minimum and
mean solve time over the requested number of repeats. It does not check
correctness: that is covered separately by
``tests/verification/test_poisson_bc.f90``, which verifies the Poisson
solve against an analytical solution across all four boundary condition
configurations. For example:

.. code-block:: console

   $ mpirun -np 2 build/bin/x3d2_poisson_solver 128 64 32 --repeat 5