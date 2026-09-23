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

The setting decides which allocator and backend sources are compiled, which
compiler and linker flags are applied, and which tests are registered. Building
more than one GPU backend at a time is not supported.

CUDA Backend
~~~~~~~~~~~~

Requires NVHPC or PGI; configuring fails with any other compiler. The build adds
``-cuda``, and links the Poisson solver against ``cuFFTMp``
(``-cudalib=cufftmp``) in an MPI build or against plain ``cuFFT``
(``-cudalib=cufft``) without one. cuFFTMp distributes a single transform across
ranks and its nvfortran wrapper calls ``MPI_Comm_f2c``, so it cannot be linked
into a serial build; a serial build is one rank, which plain cuFFT handles on
its own.

OpenMP Target Offload Backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires a compiler reporting OpenMP 4.5 or newer. Cray, GNU and NVHPC/PGI are
all supported. Each enables offload differently, and the build detects the
compiler and applies what it needs, so no manual offload flags are required with
any of them.

.. list-table::
   :header-rows: 1
   :widths: 16 44 40

   * - Compiler
     - Flags applied by the build
     - ``OMP_TGT_ARCH``
   * - Cray
     - ``-eF`` when compiling and ``-h omp`` when linking.
     - Not used. The offload target is selected by the Cray programming
       environment rather than by CMake.
   * - GNU
     - ``-fopenmp`` and ``--offload-arch=<arch>``.
     - Required. Configuring fails if it is unset.
   * - NVHPC, PGI
     - ``-mp=gpu``, replacing the host-only ``-mp`` that CMake's ``FindOpenMP``
       supplies. Without that substitution the target regions compile but run
       on the host against device pointers.
     - Optional. When unset, the compiler targets the GPU of the build machine.

Where the compiler needs the architecture, set it at configure time:

.. code-block:: bash

   -DENABLE_BACKEND=OMP_TGT -DOMP_TGT_ARCH=gfx942

Vendor selection is automatic: the build defines ``OMP_TGT_NVIDIA`` for
NVHPC/PGI and ``OMP_TGT_AMD`` for Cray and GNU, which selects the
vendor-appropriate ``SZ`` parameter. Note that the GNU path assumes an AMD
target.

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
     -DCMAKE_Fortran_COMPILER=nvfortran -DENABLE_BACKEND=CUDA

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When x3d2 is built with ``ENABLE_BACKEND=CUDA``, the built-in ADIOS2 is compiled
with CUDA support. By default its CUDA sources are built for the architecture of
the GPU present at configure time (``native`` detection). If you are building on
a node without a visible GPU (for example a GPU-less login or build node),
``native`` cannot probe the hardware, so set the target architecture explicitly
with the ``CUDA_ARCH`` CMake option:

.. code-block:: bash

   -DCUDA_ARCH=80

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

Enabling GPU-Aware ADIOS2 I/O
-----------------------------

When running on NVIDIA GPUs, x3d2 can hand fields to ADIOS2 directly from GPU memory for checkpoints and snapshots. Each field is packed into its output layout by a GPU kernel, replacing the reordering and repacking the host-staged path does on the CPU.

This does not remove the device-to-host transfer. Output files live on the host file system, so ADIOS2's BP5 engine copies each packed GPU buffer into its host serialisation buffer during ``Put``, and the same volume of data crosses PCIe as with host-staged output. The saving is the CPU work and host memory traffic around that copy.

Requirements
~~~~~~~~~~~~

GPU-aware ADIOS2 I/O requires:

- The NVHPC (or PGI) Fortran compiler
- ADIOS2 built with CUDA support (``-DADIOS2_USE_CUDA=ON``)

When using the built-in ADIOS2 (default), the build system builds ADIOS2 with CUDA support when ``ENABLE_BACKEND=CUDA`` is selected.

Build Configuration
~~~~~~~~~~~~~~~~~~~

To enable GPU-aware I/O, select the CUDA backend and enable both ADIOS2 options:

.. code-block:: bash

   cmake .. -DENABLE_BACKEND=CUDA -DWITH_ADIOS2=ON -DWITH_ADIOS2_GPU_AWARE=ON

Using a System ADIOS2 with CUDA Support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a system ADIOS2 installation that was built with ``-DADIOS2_USE_CUDA=ON``, you can use it directly:

.. code-block:: bash

   cmake .. -DENABLE_BACKEND=CUDA -DWITH_ADIOS2=ON -DUSE_SYSTEM_ADIOS2=ON -DWITH_ADIOS2_GPU_AWARE=ON

The build system will verify that the ``adios2::core_cuda`` target is available and report an error if it is not.

How Fields Are Written
~~~~~~~~~~~~~~~~~~~~~~

Snapshots and checkpoints use the GPU-aware path when every field to write lives in device memory and no output stride is set. Strided snapshots, and any field the backend cannot pack on the device, go through the host-staged path instead.

For each field, one kernel packs the solver's padded, direction-ordered storage into a contiguous buffer of the field's true extent. The same kernel converts to single precision when ``snapshot_sp`` is set. The buffer belongs to the open output file and is reused, so no device memory is allocated per write.

ADIOS2 2.12's BP5 engine copies GPU buffers to host during each ``Put``, even when a deferred put is requested. The default is therefore one staging buffer with synchronous puts. ``X3D2_ADIOS2_GPU_BATCH_FIELDS`` switches to deferred puts that keep several buffers alive until ``EndStep``. That only helps with an ADIOS2 engine that honours deferred GPU puts, and it costs one field of device memory per buffer.

The staging buffer adds device memory equal to one output field of the local subdomain, for example 128 MiB for a 256³ double-precision field on one rank. Check that this fits alongside the solver before enabling larger batches.

At startup, rank 0 prints ``ADIOS2 GPU write mode: ...``, which confirms whether the GPU-aware path is active.

Expected Performance
~~~~~~~~~~~~~~~~~~~~

The gain grows with output frequency and with the cost of CPU-side copies on the host. Strided snapshots take the host-staged path and see no gain. Checkpoints are never strided and always benefit.

As an example, a 256³ Taylor–Green vortex with IBM in double precision, run for 20 steps on one RTX A4000 with a snapshot every step at ``output_stride = 1, 1, 1`` and a checkpoint every 5 steps:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Run
     - Total (s)
     - I/O (s)
     - Peak GPU memory (MiB)
   * - Solver only, no output
     - 29.2
     - n/a
     - 3621
   * - Host-staged output
     - 76.6
     - 47.3
     - 3621
   * - GPU-aware output
     - 52.6
     - 23.4
     - 3879

I/O time roughly halves, and one snapshot step goes from about 1.5 s to 0.58 s. The time ADIOS2 spends writing the file is the same in both runs. With ``output_stride = 2, 2, 2`` in the same case, the gain is about 5%, all from the checkpoints.

Runtime Options
~~~~~~~~~~~~~~~

The ADIOS2 backend reads these environment variables when the first writer starts, and rank 0 prints the selected options. ``X3D2_NVTX`` is read separately, the first time an NVTX range is emitted.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Variable
     - Default
     - Effect
   * - ``X3D2_ADIOS2_GPU_WRITE_MODE``
     - ``auto``
     - ``auto`` (or ``gpu``) writes from device memory when possible. ``host`` always copies fields to host first, for example to compare the two paths.
   * - ``X3D2_ADIOS2_GPU_BATCH_FIELDS``
     - ``1``
     - Device staging buffers in flight. ``1`` uses one buffer with synchronous puts. ``N > 1`` uses deferred puts and calls ``PerformPuts`` whenever ``N`` buffers are pending. ``0`` keeps one buffer per field until ``EndStep``.
   * - ``X3D2_ADIOS2_IO_BENCH``
     - ``0``
     - ``1`` times ``Put`` and ``EndStep`` for every output step (maximum over ranks) and prints a summary with throughput when the file closes.
   * - ``X3D2_ADIOS2_IO_BENCH_WARMUP``
     - ``2``
     - Number of initial steps excluded from the benchmark summary.
   * - ``X3D2_ADIOS2_IO_BENCH_VERBOSE``
     - ``1``
     - ``1`` also prints the timings of every step, not only the summary.
   * - ``X3D2_NVTX``
     - ``1``
     - Emits NVTX ranges for Nsight Systems: ``ADIOS2_Put``, ``ADIOS2_EndStep`` and ``ADIOS2_DevicePack`` on the GPU-aware path, ``IO_HostStage`` and ``IO_HostPack`` on the host-staged path. CUDA builds only, with or without ADIOS2.

Boolean options accept ``1/0``, ``true/false``, ``yes/no`` and ``on/off``. For example, to benchmark snapshot output without per-step lines:

.. code-block:: bash

   X3D2_ADIOS2_IO_BENCH=1 X3D2_ADIOS2_IO_BENCH_VERBOSE=0 mpirun -np 4 ./bin/xcompact <input_file>

Profiling Output with Nsight Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CUDA builds label each output step with NVTX ranges, so the two write paths can be compared on a timeline:

.. code-block:: bash

   nsys profile -t cuda,nvtx -o gpu_io mpirun -np 1 ./bin/xcompact <input_file>
   nsys-ui gpu_io.nsys-rep

Each output step appears as ``ADIOS2_I/O_Step``. On the GPU-aware path it contains a short ``ADIOS2_DevicePack`` kernel and an ``ADIOS2_Put`` per field. The Put is mostly ADIOS2's device-to-host copy. On the host-staged path it contains ``IO_HostStage`` (device-to-host copy and reordering on the CPU) and ``IO_HostPack`` (copy into the write buffer) per field, then ``ADIOS2_Put``. Both end with ``ADIOS2_EndStep``, where ADIOS2 writes the file.

To compare the paths with a GPU-aware build only, profile a second run with ``X3D2_ADIOS2_GPU_WRITE_MODE=host``.

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
