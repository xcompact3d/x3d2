Advanced Build Configuration
============================

Configuring ADIOS2 Support
--------------------------

x3d2 can leverage the ADIOS2 library for high-performance I/O operations used by the checkpoint and snapshot system. This section explains how to build x3d2 with ADIOS2 support.

Enabling ADIOS2 Support
~~~~~~~~~~~~~~~~~~~~~~~

To build x3d2 with ADIOS2 support, use the ``ENABLE_ADIOS2`` CMake option:

.. code-block:: bash

   -DENABLE_ADIOS2=ON

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

   -DENABLE_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON

For custom installation locations, provide the path to CMake:

.. code-block:: bash

   -DENABLE_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON -DADIOS2_ROOT=/path/to/adios2/installation

When to Build a Custom ADIOS2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, you may need to build ADIOS2 specifically for use with x3d2. This is particularly important when:

- The system ADIOS2 was built with a different MPI implementation than what you're using for x3d2
- You don't have admin privileges to install ADIOS2 system-wide
- You need specific ADIOS2 features not available in your system's version

To use the built-in ADIOS2 (default behavior):

.. code-block:: bash

   -DENABLE_ADIOS2=ON

Library Path Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you already have an ADIOS2 installation but want to use the version built by x3d2, you may need to set the library path to ensure that the correct ADIOS2 libraries are used at runtime:

.. code-block:: bash

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/build/adios2/lib

This is particularly important when using a custom-built ADIOS2 with a different MPI implementation than what's available system-wide.

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

If the wrong libraries are being loaded, adjust your ``LD_LIBRARY_PATH`` environment variable:

.. code-block:: bash

   # To prioritise custom-built ADIOS2:
   export LD_LIBRARY_PATH=$(pwd)/build/adios2/lib:$LD_LIBRARY_PATH
   
   # Or to prioritise system ADIOS2 (if needed):
   export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib:$LD_LIBRARY_PATH
