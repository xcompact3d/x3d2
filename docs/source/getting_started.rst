Getting Started
===============

Installing on Linux
-------------------

Dependencies
~~~~~~~~~~~~

To build x3d2 on Linux, you will need the following dependencies:

Fortran Compiler
^^^^^^^^^^^^^^^^

A recent modern Fortran compiler is required to build x3d2. For example GNU Fortran compiler (``gfortran``) version 9 and above. You can install it using your package manager. For Ubuntu this can be obtained using:

.. code-block:: bash
   
   sudo apt install gfortran

OpenMP is typically included with most modern compilers, including GCC. Therefore, installing ``gfortran`` for example should also provide OpenMP support.

MPI
^^^

You can install Open MPI using your package manager. For Ubuntu, this can be obtained using

.. code-block:: bash

   sudo apt install openmpi-bin libopenmpi-dev

MPI is optional. Configuring with ``-DWITH_MPI=OFF`` builds a serial,
single-rank executable that needs no MPI installation at all; see
:doc:`user/advanced_build`.

CMake
^^^^^

To build x3d2, you will need CMake version 3.13 and above (3.20+ recommended). You can download the latest version from the `CMake website <https://cmake.org/download/>`_. Alternatively, you can install it using your package manager. For Ubuntu, this can be obtained using:

.. code-block:: bash

   sudo apt install cmake

NVIDIA HPC SDK
^^^^^^^^^^^^^^

If you want to install x3d2 for NVIDIA GPUs, please download and install the NVIDIA HPC SDK for your target platform from the `NVIDIA website <https://developer.nvidia.com/hpc-sdk-downloads>`_. To ensure that ``mpirun`` and ``mpif90`` pick up the version shipped with NVIDIA HPC SDK (instead of another installation) you will need to set up the environment variables (consider adding them to your ``~/.bashrc`` file if you want them to persist). For example:

.. code-block:: bash

   export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/bin:$PATH

   export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/lib:$LD_LIBRARY_PATH

.. note::

   By default, the SDK installs in ``/opt/nvidia/hpc_sdk/`` but the installation directory and version number may be different if you installed the NVIDIA HPC SDK in a custom location or used a different version. To ensure that the correct version of MPI is being used, you can check the paths of ``mpirun`` and ``mpif90`` by running:

   .. code-block:: bash

      which mpirun
      which mpif90
   
   Expected output:

   .. code-block:: bash

      /opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin/mpirun
      /opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin/mpif90
   
   If these do not match the expected paths, your system might be using a different version of MPI. You will need to update your ``PATH`` and ``LD_LIBRARY_PATH`` environment accordingly.

Compile from source
~~~~~~~~~~~~~~~~~~~~

To compile x3d2 from source, follow these steps:

1. Clone the repository and change into the x3d2 directory

.. code-block:: bash

   git clone https://github.com/xcompact3d/x3d2.git
   cd x3d2

2. Configure the build

.. code-block:: bash

   cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Release

``build`` is the directory the build configuration is written to, and
``CMAKE_BUILD_TYPE`` selects the build type (use ``Debug`` instead of
``Release`` for a debug build).

A default build uses MPI, so point ``CMAKE_Fortran_COMPILER`` at an MPI compiler
wrapper such as ``mpif90``. With ``-DWITH_MPI=OFF`` point it at the compiler
itself instead of a wrapper: ``gfortran`` for GNU, ``nvfortran`` for the NVIDIA
HPC SDK and therefore for ``ENABLE_BACKEND=CUDA``, or ``ftn`` on a Cray machine.
If you give a bare name rather than an absolute path, it has to be on your
``PATH``.

.. note::

   Pass the compiler as a ``-D`` cache variable rather than through the ``FC``
   environment variable. CMake only reads ``FC`` on the first configure of a
   fresh build directory, so ``export FC=...`` silently does nothing when
   re-configuring an existing one. ``-DCMAKE_Fortran_COMPILER=...`` is recorded
   in the cache and is what the third-party dependencies are forwarded.

3. Compile

.. code-block:: bash

   cmake --build build -j

This creates the ``xcompact`` binary in ``build/bin/``. Test executables are
placed in ``build/tests/bin/``.

4. Verify the installation by running the test suite

.. code-block:: bash

   cd build
   ctest

A successful installation should indicate 100% tests passed. ``ctest
--output-on-failure`` prints the output of any test that fails.


How the build is configured
---------------------------

The build is driven entirely by CMake cache variables passed at configure time.
There is a single configure step followed by a build step:

.. code-block:: bash

   cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpif90 <options>
   cmake --build build -j

Choosing a backend
~~~~~~~~~~~~~~~~~~

x3d2 has one CPU backend and two GPU backends, selected with ``ENABLE_BACKEND``.
The CPU (host OpenMP) backend is always built; ``ENABLE_BACKEND`` adds a GPU one
on top of it and decides which allocator, backend module and compiler flags are
used.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - ``ENABLE_BACKEND``
     - Compilers
     - Description
   * - ``OFF`` (default)
     - Cray, GNU, NVHPC, Flang
     - CPU only, using host OpenMP threading.
   * - ``CUDA``
     - NVHPC or PGI only
     - CUDA Fortran backend for NVIDIA GPUs. Configuring fails with any other
       compiler.
   * - ``OMP_TGT``
     - Cray, GNU, NVHPC, Flang
     - OpenMP target offload backend. Requires an OpenMP 4.5 or newer compiler.

"Flang" covers both upstream LLVM Flang and AMD's ROCm ``amdflang``; CMake
reports either as a ``CMAKE_Fortran_COMPILER_ID`` of ``LLVMFlang`` or
``Flang``, and the build treats the two identically. See :doc:`user/advanced_build`
for version requirements and Flang-specific notes.

For example, to build the OpenMP target offload backend:

.. code-block:: bash

   cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpif90 -DENABLE_BACKEND=OMP_TGT

Each supported compiler enables offload in its own way. The build detects the
compiler and supplies the flags it needs, so no manual offload flags are
required with any of them.

Both GPU backends need the target architecture, given with ``BACKEND_ARCH``.
Configuring a GPU build fails if it is unset. NVIDIA targets are named by their
compute capability as ``ccXX``, the convention nvfortran uses, and AMD targets
as ``gfx<model>``:

.. code-block:: bash

   -DENABLE_BACKEND=OMP_TGT -DBACKEND_ARCH=gfx942   # AMD MI300X
   -DENABLE_BACKEND=CUDA -DBACKEND_ARCH=cc80        # NVIDIA A100

The build passes the name to the selected compiler in the form it expects, for
example ``-gpu=cc80`` for NVHPC and ``--offload-arch=gfx942`` for GNU and
Flang. See :doc:`user/advanced_build` for the per-compiler details.

Build options
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 27 13 60

   * - Option
     - Default
     - Description
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - ``Release`` or ``Debug``.
   * - ``CMAKE_Fortran_COMPILER``
     - --
     - The MPI Fortran wrapper to build with, normally ``mpif90``. With
       ``WITH_MPI=OFF``, the compiler itself instead (``gfortran``,
       ``nvfortran``, ``ftn``, ``flang``) rather than a wrapper.
   * - ``ENABLE_BACKEND``
     - ``OFF``
     - GPU backend to build: ``OFF``, ``CUDA`` or ``OMP_TGT``.
   * - ``BACKEND_ARCH``
     - --
     - Target GPU architecture, required whenever ``ENABLE_BACKEND`` is not
       ``OFF``: ``cc80`` (A100), ``cc90`` (H100), ``gfx942`` (MI300X).
   * - ``SINGLE_PREC``
     - ``OFF``
     - Build in single precision.
   * - ``WITH_MPI``
     - ``ON``
     - Build against MPI. ``OFF`` builds a serial, single-rank executable.
   * - ``WITH_2DECOMPFFT``
     - ``ON``
     - Build the FFT-based Poisson solver against 2decomp-fft.
   * - ``WITH_ADIOS2``
     - ``OFF``
     - Enable ADIOS2 for checkpoint and snapshot I/O.

See :doc:`user/advanced_build` for the MPI, ADIOS2, 2decomp-fft, single
precision and CUDA debugging options in detail.

Third-party dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

2decomp-fft (when ``WITH_2DECOMPFFT=ON``) and ADIOS2 (when ``WITH_ADIOS2=ON``)
are downloaded and built into the build tree automatically. They are built
*during the build*, not during configuration, so configuring is quick and the
first ``cmake --build`` is what fetches and compiles them. Each is installed
into its own directory under ``build/`` keyed by version and configuration, and
is reused on later builds rather than rebuilt. Both can be pointed at an
existing installation instead; see :doc:`user/advanced_build`.

Build outputs
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Path
     - Contents
   * - ``build/bin/xcompact``
     - The solver executable.
   * - ``build/lib/``
     - The ``x3d2`` and ``x3d2_backends`` static libraries.
   * - ``build/tests/bin/``
     - Test executables, run via ``ctest``.


Installing on macOS
-------------------

Dependencies
~~~~~~~~~~~~

It is assumed that you have Xcode Command Line Tools and `Homebrew <https://brew.sh/>`_ installed. To build x3d2 on macOS, you will need the following dependencies:

Fortran Compiler
^^^^^^^^^^^^^^^^

A recent modern Fortran compiler is required to build x3d2. For example GNU Fortran compiler (``gfortran``) version 9 and above. You can install it using Homebrew:

.. code-block:: bash

   brew install gcc

OpenMP is typically included with most modern compilers, including GCC. Therefore, installing ``gcc`` for example should also provide OpenMP support. Next, identify your installed version (e.g., ``gcc-15``). You will need this for the next steps.

.. code-block:: bash

   ls $(brew --prefix)/bin/gcc-*

Open MPI
^^^^^^^^

You must build Open MPI from source to ensure it is compatible with GNU compilers. Replace ``15`` below with your specific GCC version.

.. code-block:: bash

   export HOMEBREW_CXX=g++-15
   export HOMEBREW_CC=gcc-15
   brew install open-mpi --build-from-source

CMake
^^^^^

To build x3d2, you will need CMake version 3.13 and above (3.20+ recommended). You can download the latest version from the `CMake website <https://cmake.org/download/>`_. Alternatively, you can install it using Homebrew:

.. code-block:: bash

   brew install cmake

Compile from source
~~~~~~~~~~~~~~~~~~~

To install x3d2 from source, follow these steps:

1. Clone the repository and change into the x3d2 directory

.. code-block:: bash

   git clone https://github.com/xcompact3d/x3d2.git

2. Configure the build

macOS users must explicitly select the GNU compilers for C and C++ so that
OpenMP support is detected correctly. Replace ``15`` with your installed GCC
version if it differs:

.. code-block:: bash

   cmake -S . -B build \
     -DCMAKE_Fortran_COMPILER=mpif90 \
     -DCMAKE_C_COMPILER=gcc-15 \
     -DCMAKE_CXX_COMPILER=g++-15 \
     -DCMAKE_BUILD_TYPE=Release

``build`` is the directory the build configuration is written to, and
``CMAKE_BUILD_TYPE`` selects the build type (use ``Debug`` for a debug build).

3. Compile

.. code-block:: bash

   cmake --build build -j

This creates the ``xcompact`` binary in ``build/bin/``.

4. Verify the installation by running the test suite

.. code-block:: bash

   cd build
   ctest

A successful installation should indicate 100% tests passed.
