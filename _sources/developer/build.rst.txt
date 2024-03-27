Building and running the unit tests
===================================

x3d2 is in early development phase, and there is no main executable
yet to be built.  However, currently implemented functionality is
covered by unit tests, which you can build and run on you development
machine.

To build x3d2, you will need git, a fortran compiler and CMake, see
:ref:`tooling`.

Start by configuring the build directory:

.. code-block:: console

   $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

To configure the build with the NVIDIA Fortran compiler, you can set
the `FC` environment variable to the compiler executable.  If you
specify an relative path, it must be present in your current `PATH`.

.. code-block:: console

   $ FC=nvfortran cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

Setting the Fortran compiler to the NVIDIA Fortran compiler will
automatically include the CUDA Fortran source files into the build
tree, which are ignored by default.

Once the build directory is configured, the tests can be built and run
as follows:

.. code-block:: console

   $ cd build
   $ make
   $ make test

Note that ``make test`` is only a launcher for the ``ctest``
executable.  By default ``ctest`` does not show the output of test
executables on failure.  If one of more tests fail, you probably want
to run the tests with:

.. code-block:: console

   $ ctest --output-on-failure

instead of ``make test``.
