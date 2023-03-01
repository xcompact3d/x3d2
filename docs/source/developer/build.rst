Building and running the unit tests
===================================

x3d2 is in early development phase, and there is no main executable
yet to be built.  However, currently implemented functionality is
covered by unit tests, which you can build and run on you development
machine.

To build x3d2, you will need git, a fortran compiler and CMake.

Start by configuring the build directory:

.. code-block:: console

   $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

Then, the tests can be built and run as follows:

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
