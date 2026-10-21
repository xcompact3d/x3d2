GPU memory check
================

``x3d2-memcheck`` estimates the GPU memory that an x3d2 case needs. It reads
an ``.x3d`` input file. It reports whether the case fits on the card, for one
GPU and for larger GPU counts. Run it before you submit a job.

The tool exists only in the CUDA backend. It needs cuFFTMp. Tested on an
A100 40 GB with NVHPC 25.3 and CUDA 12.8.

Building
--------

Configure the build with ``-DENABLE_BACKEND=CUDA`` and build the target
``x3d2-memcheck``:

.. code-block:: bash

   cmake -S . -B build -DENABLE_BACKEND=CUDA -DBACKEND_ARCH=cc80
   cmake --build build --target x3d2-memcheck

``x3d2-memcheck`` needs an MPI build (``WITH_MPI=ON``) with cuFFTMp. A
``WITH_MPI=OFF`` build does not make the target. CMake prints a message and
skips it.

The binary is ``build/bin/x3d2-memcheck``. See :doc:`advanced_build` for the
CUDA backend options.

Running
-------

Use one MPI rank. The input file must be the first argument:

.. code-block:: bash

   mpirun -n 1 bin/x3d2-memcheck <input.x3d> [--static | --build [--extensive <n>]]

``--static``
   Run the formula estimate and a small calibration probe.

``--build``
   Run the default mode, then a real case build.

``--extensive <n>``
   Run ``<n>`` substeps in the real build. ``<n>`` must be a positive
   integer. The option also selects ``--build``.

The tool stops with an error when you give ``--static`` together with
``--build`` or ``--extensive``. It also stops with an error for a missing
input file argument, an unknown flag, a run with more than one MPI rank, and
a run without a usable CUDA device.

Modes
-----

The tool has three modes.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Mode
     - Description
   * - ``--static``
     - The formula estimate and one calibration probe on a 32^3 grid. The
       tool builds no full-grid plans.
   * - default
     - The static mode, then a query-only cuFFTMp plan on the full grid. The
       query allocates nothing and gives the plan size. The tool then runs
       the full-grid probe if the free memory holds the plan size plus a
       margin.
   * - ``--build``
     - The default mode, then a real case build and one substep on one GPU
       (``<n>`` substeps with ``--extensive``). The tool measures the memory
       with ``cudaMemGetInfo``. It prints a ``CHECK`` line that compares the
       measurement with the estimate for ``num_gpus=1``.

The tool adds the GPU-aware ADIOS2 staging buffer by formula in every mode. The
real build does not write a snapshot, so the ``CHECK`` line excludes this term.

The ``CHECK`` line shows ``OK`` when the estimate and the measurement differ
by 5 percent or less. Otherwise it shows ``MISMATCH``.

The tool skips the real build, and prints no ``CHECK`` line, in these cases:

* The estimate for ``num_gpus=1`` is ``DOES_NOT_FIT``.
* The input names a ``flow_case`` that the tool does not know. The tool logs
  that the flow case has no dispatch.
* The input sets ``ibm_on=T`` and the ``ibm_<BC-suffix>.bp`` mask file is
  not in the working directory.

The estimate then stands.

The real build writes no case output files. It switches off ``monitoring.csv``,
the turbine ``disc<t>.adm`` files and the ABL profile file.

Reading the table
-----------------

The report has one row for each GPU count: 1, 2, 4 and 8. The ``basis``
column tells how the tool obtained the FFT overhead of the row:

``measured``
   The row has ``num_gpus=1``. The value comes from the probe, or from the real
   build.

``extrapolated``
   The row has ``num_gpus>1`` and a full-grid probe or measurement exists. The
   tool scales the one-GPU result by 1/``num_gpus``. Without one, the row is
   ``formula-only`` or ``predicted``.

``predicted``
   The full-grid probe did not fit on the card. The tool uses the heap size
   that the plan query predicts. The report says
   ``no full-grid probe: ...`` with the reason.

``formula-only``
   The tool uses the formula alone. No full-grid probe or measurement backs
   the row. The 32^3 calibration probe always runs. A row is also
   ``formula-only`` when its formula floor is above 95 percent of the card
   and the tool has no predicted heap size. With a predicted heap size, such
   a row is ``predicted``.

Verdicts
--------

The verdict depends on the memory on each GPU as a fraction of the card
memory:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Verdict
     - Fraction
     - Meaning
   * - ``FITS``
     - below 80%
     - The case fits with margin.
   * - ``BORDERLINE``
     - 80% to 95%
     - The case can fit. Run ``--build`` to measure it.
   * - ``DOES_NOT_FIT``
     - above 95%
     - The case does not fit.

The verdict ``UNSUPPORTED`` means that x3d2 does not support the GPU count or
the decomposition. A GPU count that is not supported has no table row. The
tool prints ``<num_gpus> (skipped: <reason>)`` instead. For the layout in the
input file, the tool prints ``Requested nproc_dir [..]: not supported
(<reason>)``. The reasons are:

* ``num_gpus`` is below 1.
* ``num_gpus`` is above 1 and ``BC_y`` is not periodic.
* ``gdims(3)`` does not divide by ``num_gpus``.
* The y cell count does not divide by ``num_gpus`` (BC 000 and 010).
* The x cell count does not divide by ``num_gpus`` (BC 100).
* ``nproc_dir`` is not ``[1,1,num_gpus]``.

Exit codes
----------

The exit code follows the verdict of the ``nproc_dir`` in the input file:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Code
     - Verdict
   * - 0
     - ``FITS``
   * - 1
     - ``BORDERLINE``
   * - 2
     - ``DOES_NOT_FIT`` or ``UNSUPPORTED``

A run that aborts ends with a non-zero exit code from ``abort_run``. This
happens for a usage error, for a runtime error and for a plain cuFFT fallback.

Limitations
-----------

.. warning::

   The estimate is approximate. A job can still run out of memory after the
   check gives ``FITS``.

* The rows for ``num_gpus>1`` are extrapolated when a full-grid probe or
  measurement exists. The tool scales it by 1/``num_gpus``. It does not
  measure them.
* The full-grid probe guard is approximate. The tool compares the free memory
  with the query-only plan size plus a margin of 1 GiB. The sizes come from
  measurements on an A100 40 GB, and other cards can differ. The guard can
  under-estimate the memory. The probe can then abort for lack of memory.
* The tool needs cuFFTMp. It aborts in every mode when an FFT plan asks for
  cuFFTMp and gets plain cuFFT.
* ``x3d2-memcheck`` needs an MPI build (``WITH_MPI=ON``). A ``WITH_MPI=OFF``
  build does not make the target.
