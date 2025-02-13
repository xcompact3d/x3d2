Building x3d2
=============

See :ref:`tooling` for details on the tool required to build x3d2.

Start by configuring the build directory:

.. code-block:: console

   $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

The above is using ``Debug`` build, for release build use ``Release`` instead.

Set the ``FC`` environment variable to the compiler executable (this can either be from NVIDIA HPC SDK or Open MPI). If you specify a relative path, it must be present in your current ``PATH``.

.. code-block:: bash

   $ FC=mpif90 cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

Once the build directory is configured, you can build the executable and run the tests as follows:

.. code-block:: bash

   $ cd build
   $ make
   $ make test

Note that ``make test`` is only a launcher for the ``ctest`` executable. By default, ``ctest`` does not show the output of test executables on failure. If one or more tests fail, you probably want to run the tests with:

.. code-block:: bash

   $ ctest --output-on-failure

instead of ``make test``.

The main executable will be built in the ``build/src`` directory as ``xcompact``.