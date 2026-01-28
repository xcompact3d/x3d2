Advanced Build Configuration
============================

Configuring ADIOS2 Support
--------------------------

x3d2 can leverage the ADIOS2 library for high-performance I/O operations used by the checkpoint and snapshot system. This section explains how to build x3d2 with ADIOS2 support.

Enabling ADIOS2 Support
~~~~~~~~~~~~~~~~~~~~~~~

To build x3d2 with ADIOS2 support, use the ``WITH_ADIOS2`` CMake option:

.. code-block:: bash

   -DWITH_ADIOS2=ON

ADIOS2 Installation Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

x3d2 provides two approaches for using ADIOS2:

1. **Built-in ADIOS2 (Default)**: By default, x3d2's build system automatically downloads and builds ADIOS2 for you.

2. **System ADIOS2**: Use an existing ADIOS2 installation on your system.

Using an Existing ADIOS2 Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, x3d2 will not use an ADIOS2 installation already present on your system. If you want to use an existing ADIOS2 installation, you need to set ``-DUSE_SYSTEM_ADIOS2=ON``.

If ADIOS2 is installed in a standard location, no additional configuration is needed:

.. code-block:: bash

   -DWITH_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON

For custom installation locations, provide the path to CMake:

.. code-block:: bash

   -DWITH_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON -DADIOS2_ROOT=/path/to/adios2/installation

When to Build a Custom ADIOS2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, you may need to build ADIOS2 specifically for use with x3d2. This is particularly important when:

- The system ADIOS2 was built with a different MPI implementation than what you're using for x3d2
- You don't have admin privileges to install ADIOS2 system-wide
- You need specific ADIOS2 features not available in your system's version

To use the built-in ADIOS2 (default behavior):

.. code-block:: bash

   -DWITH_ADIOS2=ON

Library Path Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

The project is configured to automatically download and build its own version of the ADIOS2 library. However, if you have another version of ADIOS2 already installed globally on your system, the runtime linker might mistakenly load the system's version.
This can lead to ``undefined symbol`` errors if the system's ADIOS2 was built with a different compiler than the one used for this project, or with an incompatible MPI implementation.

If you have ParaView installed, it often includes its own version of ADIOS2 which may conflict with the project's custom-built ADIOS2.

Prepending to LD_LIBRARY_PATH
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The safest way to prioritise your project's ADIOS2 while keeping access to other system libraries:

1. Navigate into your build directory:

   .. code-block:: bash

      cd <path-to-your-build-directory>

2. Run commands by prepending the project's ADIOS2 library path:

   .. code-block:: bash

      # For test suite
      LD_LIBRARY_PATH=./adios2-build/adios2-install/lib:$LD_LIBRARY_PATH make test

      # For running with mpirun
      LD_LIBRARY_PATH=./adios2-build/adios2-install/lib:$LD_LIBRARY_PATH mpirun -np 2 ./src/xcompact <input_file>

Full Replacement (Aggressive Isolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you continue to experience conflicts even with prepending, you can completely replace ``LD_LIBRARY_PATH``:

.. code-block:: bash

   # From build directory - replaces LD_LIBRARY_PATH entirely
   LD_LIBRARY_PATH=./adios2-build/adios2-install/lib mpirun -np 2 ./src/xcompact <input_file>

When to use this:

- When prepending doesn't resolve the conflict
- When you're certain your build uses RPATH for other dependencies (typical with modern CMake)

Caution:

- This removes all other paths from ``LD_LIBRARY_PATH``
- Only use if you understand your executable's dependency structure
- Your MPI, system libraries, etc. should be found via RPATH or system default paths

Verifying Your Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to verify your ADIOS2 installation is to run the test suite:

.. code-block:: bash

   make test

This will run a set of tests including ADIOS2 functionality tests. Look for passing tests related to checkpoint I/O and ADIOS2 operations.

You can also verify functionality by:

1. Creating a checkpoint namelist in your input file
2. Running a simulation with checkpointing enabled
3. Checking that checkpoint files are correctly generated in your output directory

If you encounter errors about missing libraries at runtime, check that the correct library path is set and that compatible MPI libraries are being used by both x3d2 and ADIOS2.

Troubleshooting
~~~~~~~~~~~~~~~

Checking Library Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter issues with ADIOS2 libraries, you can check which libraries x3d2 is actually using with the ``ldd`` command:

.. code-block:: bash

   ldd ./build/src/xcompact | grep adios2

This will show all the ADIOS2 libraries being loaded and their paths. Make sure they point to the expected location (either your system libraries or the custom-built ones).

Common issues include:

- Wrong ADIOS2 library is being loaded (system instead of custom or vice versa)
- MPI library mismatch between ADIOS2 and x3d2
- Missing libraries (shown as "not found")

Configuring Single Precision Mode
---------------------------------

x3d2 can be compiled to use single precision (32-bit) floating-point numbers as the default precision for all calculations, which can provide performance benefits and memory savings on some hardware.

Enabling Single Precision
~~~~~~~~~~~~~~~~~~~~~~~~~

To compile x3d2 in single precision mode, use the ``SINGLE_PREC`` CMake option:

.. code-block:: bash

   cmake -DSINGLE_PREC=ON ..

This will define the ``SINGLE_PREC`` preprocessor macro, causing the code to use single precision (``real32``) as the default floating-point type throughout the application.

Single Precision and Snapshot Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IO precision depends on two factors:

1. Compile-time precision (``-DSINGLE_PREC=ON``): Controls simulation precision. All I/O (checkpoints and snapshots) uses the simulation precision.

2. Runtime snapshot precision (``snapshot_sp=.true.`` in input file): Only available when compiled in double precision. Converts snapshots to single precision while keeping simulation and checkpoints in double precision.

Available combinations:

- Double precision simulation (default): checkpoints in double precision, snapshots configurable via ``snapshot_sp``
- Single precision simulation (``-DSINGLE_PREC=ON``): all I/O in single precision, ``snapshot_sp`` ignored
