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

Each test is registered with a CTest label matching its category (``unit``, ``verification``, or ``performance``) and its backend (``omp``, ``cuda`` or ``omp_tgt``).

Writing a test
--------------

Create a Fortran source file in the appropriate directory. Tests are backend
agnostic: the same source is compiled once per backend, so write it against
``base_backend_t`` and let ``m_backend_runtime`` (``tests/common/backend_runtime.f90``)
build the backend and allocator for the current build:

.. code-block:: fortran

   use m_backend_runtime, only: backend_runtime_t, backend_sz
   ...
   type(backend_runtime_t), target :: runtime
   class(base_backend_t), pointer :: backend
   ...
   call runtime%init(mesh)
   backend => runtime%backend

- Move data between host arrays and fields with ``backend%set_field_data`` and
  ``backend%get_field_data`` (or ``field%fill``). Do not read or write
  ``field%data`` on a backend allocator's field: on a GPU backend it is not where
  the data lives.
- Fill data in Cartesian order (the default orientation of ``set_field_data``)
  rather than computing a directional layout by hand; the layouts in memory differ
  between backends.
- Build operators with ``backend%alloc_tdsops`` rather than ``tdsops_init``.
- ``backend_sz`` and ``backend_is_cuda`` are there for sizing a problem to the
  backend; keep backend-specific code out of the test itself. If a test needs a
  kernel that ``base_backend_t`` does not expose yet, add a small helper to
  ``m_backend_runtime`` (see ``penta_solve``) instead of ``#ifdef``-ing the test.

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

Add a line to the ``define_backend_tests`` function in ``tests/CMakeLists.txt``,
using the function matching your category:

.. code-block:: cmake

   # Unit test, 1 MPI rank
   define_test(unit/test_example.f90 1 ${backend})

   # Verification test, CMAKE_CTEST_NPROCS MPI ranks
   define_verification_test(verification/test_example.f90 ${np} ${backend})

   # Performance test, 1 MPI rank
   define_performance_test(performance/perf_example.f90 1 ${backend})

``define_backend_tests`` is called once for every backend in the build: always
``omp`` (CPU), plus ``cuda`` or ``omp_tgt`` when ``ENABLE_BACKEND`` selects one.
A test that cannot run on some backend because of a known feature gap goes
inside an ``if`` on ``${backend}`` with a comment saying why.

Tests of host-side code that never touches a backend (mesh, statistics, ...) are
registered once, with ``omp``, above that function.

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
- Backend independence: no ``#ifdef CUDA`` / ``OMP_TGT`` in tests; backend-specific
  code lives in ``tests/common/backend_runtime.f90``. The only exceptions are tests
  of a feature that exists on one backend alone, such as GPU-aware I/O
  (``test_cuda_gpu_aware_io*``), which are registered for that backend only.
- MPI: All tests call ``MPI_Init``/``MPI_Finalize``, even single-rank tests. Tests are
  launched with ``mpirun --oversubscribe -np <N>``, or with ``srun -n <N>`` when building
  with the Cray compiler, since Cray machines are Slurm-driven and have no ``mpirun``.
- Shared utilities: All tests automatically link against ``x3d2_test_utils`` (provides ``check_norm`` and ``report_perf``).
