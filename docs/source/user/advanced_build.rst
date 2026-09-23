Advanced Build Configuration
============================

This page covers the individual build options in detail. See
:doc:`../getting_started` for the overall CMake workflow and a summary table of
every option.

Selecting a Backend
-------------------

``ENABLE_BACKEND`` chooses which GPU backend is compiled alongside the CPU
(host OpenMP) backend, which is always built:

.. code-block:: bash

   -DENABLE_BACKEND=OFF       # CPU only (default)
   -DENABLE_BACKEND=CUDA      # CUDA Fortran, NVHPC/PGI only
   -DENABLE_BACKEND=OMP_TGT   # OpenMP target offload

Both GPU backends additionally require the target architecture in
``BACKEND_ARCH`` (see the per-backend sections below); only the default ``OFF``
builds without it.

The setting decides which allocator and backend sources are compiled, which
compiler and linker flags are applied, and which tests are registered. Building
more than one GPU backend at a time is not supported.

CUDA Backend
~~~~~~~~~~~~

Requires NVHPC or PGI; configuring fails with any other compiler. The compute
capability is required and comes from ``BACKEND_ARCH``, given in NVIDIA's
``sm_<cc>`` form; configuring fails if it is unset or does not start with
``sm_``. The build adds ``-cuda -gpu=cc<cc>`` when compiling and linking, and
links the Poisson solver against ``cuFFTMp``
(``-cudalib=cufftmp``) in an MPI build or against plain ``cuFFT``
(``-cudalib=cufft``) without one. cuFFTMp distributes a single transform across
ranks and its nvfortran wrapper calls ``MPI_Comm_f2c``, so it cannot be linked
into a serial build; a serial build is one rank, which plain cuFFT handles on
its own.

OpenMP Target Offload Backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires a compiler reporting OpenMP 4.5 or newer. Cray, GNU, NVHPC/PGI and
Flang are all supported. Each enables offload differently, and the build
detects the compiler and applies what it needs, so no manual offload flags are
required with any of them.

.. list-table::
   :header-rows: 1
   :widths: 16 44 40

   * - Compiler
     - Flags applied by the build
     - ``BACKEND_ARCH``
   * - Cray
     - ``-eF`` when compiling and ``-h omp`` when linking.
     - Must be set, but is not passed on: the offload target is selected by the
       Cray programming environment rather than by CMake.
   * - GNU
     - ``-fopenmp`` and ``--offload-arch=<arch>``.
     - ``gfx<model>``. Configuring fails on any other form.
   * - NVHPC, PGI
     - ``-mp=gpu -gpu=cc<cc>``, replacing the host-only ``-mp`` that CMake's
       ``FindOpenMP`` supplies. Without that substitution the target regions
       compile but run on the host against device pointers.
     - ``sm_<cc>``, translated to nvfortran's ``cc<cc>``. Configuring fails on
       any other form.
   * - Flang
     - ``-fopenmp --offload-arch=<arch>``, plus ``-fopenmp-version=50`` (see
       `Building with Flang`_ below).
     - ``gfx<model>``. Configuring fails on any other form.

``BACKEND_ARCH`` is required for every GPU build, whatever the compiler, and is
set at configure time:

.. code-block:: bash

   -DENABLE_BACKEND=OMP_TGT -DBACKEND_ARCH=gfx942

Vendor selection is automatic: the build defines ``OMP_TGT_NVIDIA`` for
NVHPC/PGI and ``OMP_TGT_AMD`` for Cray, GNU and Flang, which selects the
vendor-appropriate ``SZ`` parameter. Note that the GNU and Flang paths assume
an AMD target.

Building with Flang
--------------------

x3d2 builds with both upstream LLVM Flang and AMD's ROCm ``amdflang``; CMake
reports either as a ``CMAKE_Fortran_COMPILER_ID`` of ``LLVMFlang`` or
``Flang``, and the build treats the two identically.

Minimum version
~~~~~~~~~~~~~~~

Compiling the OpenMP target offload backend (``ENABLE_BACKEND=OMP_TGT``)
requires **flang-22 or newer**. Earlier releases cannot compile the OpenMP
constructs used by the offload sources, even with the version flag described
below. A CPU-only build (``ENABLE_BACKEND=OFF``) has no such requirement
beyond a Flang release recent enough to support ``-std=f2018``.

OpenMP version flag
~~~~~~~~~~~~~~~~~~~~

Flang defaults to parsing OpenMP 3.1, which rejects both ``declare target`` in
a pure procedure (added in OpenMP 4.5) and the ``loop`` construct (added in
OpenMP 5.0) that the sources use. The build always adds
``-fopenmp-version=50`` for Flang, for both the CPU and ``OMP_TGT`` backends,
so this does not need to be set manually.

Locating the OpenMP runtime (``libomp.so``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Flang links against ``libomp.so`` out of its own toolchain installation, which
is usually not on the loader's default search path. Neither the Flang driver
nor CMake's ``FindOpenMP`` (which reports only ``-fopenmp`` for Flang) leaves
a library path behind for CMake to turn into an RPATH, so binaries would link
successfully but fail to start with an error such as ``libomp.so: cannot open
shared object file``.

To avoid depending on ``LD_LIBRARY_PATH`` at runtime, the build locates
``libomp.so`` itself and adds its directory to the RPATH of every target that
links OpenMP. It first asks the compiler driver
(``flang -print-file-name=libomp.so``), then falls back to searching
``lib``/``lib64``/``lib/<triple>`` next to the compiler's install prefix,
which covers both a plain LLVM install and the ROCm layout.

If detection fails, configuring emits a warning and the resulting binaries
will only run with ``libomp.so``'s directory on ``LD_LIBRARY_PATH``. Point the
build at it directly instead:

.. code-block:: bash

   -DFLANG_OPENMP_RUNTIME_DIR=/path/to/directory/containing/libomp.so

Building without MPI
--------------------

MPI is enabled by default. Turning it off builds a serial executable that runs
as a single rank and needs no MPI installation or launcher:

.. code-block:: bash

   -DWITH_MPI=OFF

``CMAKE_Fortran_COMPILER`` is then the compiler itself rather than an MPI
wrapper. Do not pass ``mpif90`` or ``mpifort`` here: the wrapper puts the MPI
library back on the compile and link lines, which defeats the point of the
option. Name the compiler that matches the backend, for example ``nvfortran``
for ``ENABLE_BACKEND=CUDA``, ``ftn`` on a Cray machine, or ``gfortran`` for a
CPU build with GNU:

.. code-block:: bash

   cmake -S . -B build -DWITH_MPI=OFF \
     -DCMAKE_Fortran_COMPILER=nvfortran -DENABLE_BACKEND=CUDA -DBACKEND_ARCH=sm_80

The code still calls MPI unconditionally; ``src/mpi.f90`` supplies serial
stand-ins for the MPI entities x3d2 uses, under which every collective is the
identity, the rank is 0 and the communicator size is 1. Add new stubs there
rather than guarding call sites.

Three things follow from the single rank:

* ``WITH_2DECOMPFFT`` is forced off, because 2decomp-fft is an MPI library. The
  FFT-based Poisson solver is therefore unavailable in a serial build.
* ``ENABLE_BACKEND=CUDA`` links plain ``cuFFT`` rather than ``cuFFTMp``, as
  described under `CUDA Backend`_ above. This is transparent: the CUDA Poisson
  solver already falls back to plain cuFFT at runtime wherever cuFFTMp is
  unavailable.
* Only the single-rank tests are registered, and they run the executable
  directly instead of through ``mpirun``.

``WITH_ADIOS2`` still works: ADIOS2 is built without MPI to match, and x3d2
uses its serial bindings. See :ref:`adios2-and-mpi` below.

Configuring 2decomp-fft Support
-------------------------------

The FFT-based Poisson solver is built against `2decomp-fft
<https://github.com/xcompact3d/2decomp-fft>`_ and is enabled by default. Disable
it if you do not need the FFT Poisson solver, which also avoids the download:

.. code-block:: bash

   -DWITH_2DECOMPFFT=OFF

By default the library is downloaded and built into the build tree, into a
directory keyed by version and precision (``decomp2d-opt-dp-<version>``, or
``-sp-`` for a ``SINGLE_PREC=ON`` build) so that a double and a single precision
build do not collide. To use a pre-built installation instead, point
``decomp2d_install_dir`` at its prefix:

.. code-block:: bash

   -Ddecomp2d_install_dir=/path/to/2decomp-fft/install

The library must match the precision x3d2 is built with.

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

The built-in ADIOS2 is fetched and compiled *during the build*, not during
configuration, so ``cmake -S . -B build`` returns quickly and the first
``cmake --build build`` is what downloads and builds it. It is reused on
subsequent builds; to deliberately re-fetch the pinned tag, build the
``adios2-<version>-update`` target.

Using an Existing ADIOS2 Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, x3d2 will not use an ADIOS2 installation already present on your system. If you want to use an existing ADIOS2 installation, you need to set ``-DUSE_SYSTEM_ADIOS2=ON``. Configuring fails if ``USE_SYSTEM_ADIOS2=ON`` but no installation is found.

If ADIOS2 is installed in a standard location, no additional configuration is needed:

.. code-block:: bash

   -DWITH_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON

For custom installation locations, provide the path to CMake:

.. code-block:: bash

   -DWITH_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON -DADIOS2_ROOT_DIR=/path/to/adios2/installation

.. _adios2-and-mpi:

ADIOS2 and MPI
~~~~~~~~~~~~~~

ADIOS2 follows ``WITH_MPI``. The built-in ADIOS2 is configured with
``ADIOS2_USE_MPI`` to match, and the build links the Fortran bindings ADIOS2
exports for that configuration: ``adios2::fortran_mpi`` with MPI,
``adios2::fortran`` without it.

The two must agree, because ADIOS2's ``adios2_init`` and ``adios2_open`` take a
communicator only in an MPI build; without one they take no communicator at
all, which is a difference in the argument list rather than in a value that
could be passed through. x3d2 keeps the communicator flowing through its own
I/O layer regardless and drops it in the two wrappers at the top of
``src/io/adios2/io.f90``.

With ``USE_SYSTEM_ADIOS2=ON`` in an MPI build, the MPI component is required,
so configuring fails against a serial installation rather than failing to
compile later.

When to Build a Custom ADIOS2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, you may need to build ADIOS2 specifically for use with x3d2. This is particularly important when:

- The system ADIOS2 was built with a different MPI implementation than what you're using for x3d2
- You don't have admin privileges to install ADIOS2 system-wide
- You need specific ADIOS2 features not available in your system's version

To use the built-in ADIOS2 (default behaviour):

.. code-block:: bash

   -DWITH_ADIOS2=ON

CUDA Architecture for the ADIOS2 Build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When x3d2 is built with ``ENABLE_BACKEND=CUDA``, the built-in ADIOS2 is compiled
with CUDA support. It reuses ``BACKEND_ARCH``, stripped to the bare compute
capability that ``CMAKE_CUDA_ARCHITECTURES`` expects, so ``-DBACKEND_ARCH=sm_80``
builds ADIOS2's CUDA sources for ``80``. There is no separate option to set, and
nothing is probed from the build machine, so a GPU-less login or build node
configures the same way as a compute node.

Library Path Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

The project is configured to automatically download and build its own version of the ADIOS2 library. However, if you have another version of ADIOS2 already installed globally on your system, the runtime linker might mistakenly load the system's version.
This can lead to ``undefined symbol`` errors if the system's ADIOS2 was built with a different compiler than the one used for this project, or with an incompatible MPI implementation.

If you have ParaView installed, it often includes its own version of ADIOS2 which may conflict with the project's custom-built ADIOS2.

.. note::

   The built-in ADIOS2 is installed into ``adios2-<config>-<version>`` inside
   the build directory, where ``<config>`` is ``cuda`` for a
   ``ENABLE_BACKEND=CUDA`` build and ``cpu`` otherwise. A CUDA-enabled ADIOS2
   cannot be reused by a CPU build, which is why the two are kept apart. Adjust
   the paths below to match your build.

Prepending to LD_LIBRARY_PATH
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The safest way to prioritise your project's ADIOS2 while keeping access to other system libraries:

1. Navigate into your build directory:

   .. code-block:: bash

      cd <path-to-your-build-directory>

2. Run commands by prepending the project's ADIOS2 library path:

   .. code-block:: bash

      # For test suite
      LD_LIBRARY_PATH=./adios2-cpu-v2.12.1/lib:$LD_LIBRARY_PATH ctest

      # For running with mpirun
      LD_LIBRARY_PATH=./adios2-cpu-v2.12.1/lib:$LD_LIBRARY_PATH mpirun -np 2 ./bin/xcompact <input_file>

Full Replacement (Aggressive Isolation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you continue to experience conflicts even with prepending, you can completely replace ``LD_LIBRARY_PATH``:

.. code-block:: bash

   # From build directory - replaces LD_LIBRARY_PATH entirely
   LD_LIBRARY_PATH=./adios2-cpu-v2.12.1/lib mpirun -np 2 ./bin/xcompact <input_file>

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

   ctest

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

   ldd ./build/bin/xcompact | grep adios2

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

   -DSINGLE_PREC=ON

This will define the ``SINGLE_PREC`` preprocessor macro, causing the code to use single precision (``real32``) as the default floating-point type throughout the application.

Single Precision and Snapshot Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IO precision depends on two factors:

1. Compile-time precision (``-DSINGLE_PREC=ON``): Controls simulation precision. All I/O (checkpoints and snapshots) uses the simulation precision.

2. Runtime snapshot precision (``snapshot_sp=.true.`` in input file): Only available when compiled in double precision. Converts snapshots to single precision while keeping simulation and checkpoints in double precision.

Available combinations:

- Double precision simulation (default): checkpoints in double precision, snapshots configurable via ``snapshot_sp``
- Single precision simulation (``-DSINGLE_PREC=ON``): all I/O in single precision, ``snapshot_sp`` ignored

Debugging the CUDA Backend
--------------------------

FP Exception Trapping
~~~~~~~~~~~~~~~~~~~~~~

The NVHPC debug build omits ``-Ktrap=fp`` by default. Because it traps floating-point exceptions 
process-wide, benign FP operations inside linked libraries (for example cuFFTMp's PTX-JIT during 
Poisson plan creation, or HPC-X/UCX in ``ucp_worker_iface_open`` during ``MPI_Init``) can raise 
a fatal ``SIGFPE`` that is unrelated to x3d2. To enable trapping when hunting a genuine FP bug 
in x3d2's own kernels, use the ``CUDA_DEBUG_FP_TRAP`` CMake option:

.. code-block:: bash

   -DCUDA_DEBUG_FP_TRAP=ON

