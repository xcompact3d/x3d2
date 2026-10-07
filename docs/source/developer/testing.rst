Testing
=======

Test categories
---------------

Tests are organised into three directories by purpose:

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Category
     - Directory
     - When to use
   * - Unit
     - ``tests/unit/``
     - Fast checks of a single module (allocator, reordering, statistics)
   * - Verification
     - ``tests/verification/``
     - Comparing numerical results against analytical solutions
   * - Performance
     - ``tests/performance/``
     - Benchmarking throughput — large problems, many iterations, no correctness checks

Each test is registered with a CTest label matching its category (``unit``, ``verification``, or ``performance``) and its backend (``omp`` or ``cuda``).

Writing a test
--------------

Create a Fortran source file in the appropriate directory. Use ``#ifdef CUDA`` guards so the same file compiles for both backends.

**Unit test example** (``tests/unit/test_example.f90``):

.. code-block:: fortran

   program test_example
     use MPI
     use m_common, only: dp
     implicit none

     integer :: ierr
     logical :: allpass

     call MPI_Init(ierr)
     allpass = .true.

     ! --- your checks here ---
     if (1 + 1 /= 2) then
       print *, 'FAIL: arithmetic'
       allpass = .false.
     end if

     call MPI_Finalize(ierr)
     if (.not. allpass) error stop 1
   end program test_example

**Verification test** — use ``check_norm`` from ``m_test_utils``:

.. code-block:: fortran

   use m_test_utils, only: check_norm
   ...
   call check_norm(error_norm, 1e-8_dp, 'my_operator_periodic', allpass)

``check_norm`` prints a standardised ``PASSED``/``FAILED`` line with the norm value and tolerance.

**Performance test** — use ``report_perf`` from ``m_test_utils``:

.. code-block:: fortran

   use m_test_utils, only: report_perf
   ...
   call report_perf('my_kernel', elapsed_time, n_iters, ndof, bytes_per_dof)

This emits machine-parseable output: ``PERF_METRIC: <label> time=<X>s bw=<Y> GiB/s``

Registering a test
------------------

Add a line to ``tests/CMakeLists.txt`` using the function matching your category:

.. code-block:: cmake

   # Unit test, 1 MPI rank, OMP backend
   define_test(unit/test_example.f90 1 omp)

   # Verification test, 4 MPI ranks, OMP backend
   define_verification_test(verification/test_example.f90 4 omp)

   # Performance test, 1 MPI rank, CUDA backend
   # (inside the ENABLE_BACKEND STREQUAL "CUDA" block)
   define_performance_test(performance/perf_example.f90 1 cuda)

The backend argument is one of ``omp`` (CPU), ``cuda`` or ``omp_tgt``. Tests for
a GPU backend must be registered inside the matching ``ENABLE_BACKEND`` block at
the end of the file, since those sources are only compiled when that backend is
selected.

For dual-backend tests, register twice — once in the OMP section and once inside
the relevant ``ENABLE_BACKEND`` block:

.. code-block:: cmake

   # OMP section
   define_test(unit/test_example.f90 1 omp)

   # Inside the ENABLE_BACKEND STREQUAL "CUDA" block
   define_test(unit/test_example.f90 1 cuda)

   # Inside the ENABLE_BACKEND STREQUAL "OMP_TGT" block
   define_test(unit/test_example.f90 1 omp_tgt)

Running tests
-------------

.. code-block:: bash

   $ cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Debug
   $ cd build
   $ make
   $ ctest -R test_example --output-on-failure

Common CTest commands:

.. code-block:: bash

   # All unit + verification tests (same as CI)
   $ ctest -L "unit|verification" --output-on-failure

   # By category
   $ ctest -L unit
   $ ctest -L verification
   $ ctest -L performance

   # By backend
   $ ctest -L omp
   $ ctest -L cuda
   $ ctest -L omp_tgt

   # List all tests without running
   $ ctest -N

Conventions
-----------

- File naming: ``test_*.f90`` for unit/verification, ``perf_*.f90`` for performance.
- Exit code: Call ``error stop 1`` on failure so CTest detects it.
- Backend guards: Use ``#ifdef CUDA`` / ``#elif defined(OMP_TGT)`` / ``#else`` / ``#endif`` for backend-specific code.
- MPI: All tests call ``MPI_Init``/``MPI_Finalize``, even single-rank tests. Tests are
  launched with ``mpirun --oversubscribe -np <N>``, or with ``srun -n <N>`` when building
  with the Cray compiler, since Cray machines are Slurm-driven and have no ``mpirun``.
- Shared utilities: All tests automatically link against ``x3d2_test_utils`` (provides ``check_norm`` and ``report_perf``).
