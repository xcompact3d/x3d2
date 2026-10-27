set(CMAKE_CTEST_ARGUMENTS "--output-on-failure" CACHE STRING "comma separated list of arguments to pass to ctest")
mark_as_advanced(CMAKE_CTEST_ARGUMENTS)
set(CMAKE_CTEST_NPROCS "4" CACHE STRING "number of process used in tests")
mark_as_advanced(CMAKE_CTEST_NPROCS)
if(NOT WITH_MPI)
  # A build without MPI is a single rank, so the parametrised tests run at one
  # process and the explicitly multi-rank ones below are not registered at all.
  set(CMAKE_CTEST_NPROCS 1)
endif()
option(WITH_NIGHTLY_REGRESSION
  "Register the heavier nightly regression cases (excluded from PR CI)." OFF)
mark_as_advanced(WITH_NIGHTLY_REGRESSION)

# Shared test utilities library
add_library(x3d2_test_utils STATIC common/test_utils.f90)
target_link_libraries(x3d2_test_utils PRIVATE x3d2)
target_link_libraries(x3d2_test_utils PRIVATE OpenMP::OpenMP_Fortran)
# test_utils.f90 reaches MPI through m_mpi, so with MPI enabled it needs the
# MPI module path itself: x3d2 links MPI privately and does not pass it on.
if(WITH_MPI)
  target_link_libraries(x3d2_test_utils PRIVATE MPI::MPI_Fortran)
endif()
set_target_properties(x3d2_test_utils PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_test_utils")
target_include_directories(x3d2_test_utils INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_test_utils")
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  target_compile_options(x3d2_test_utils PRIVATE "-eF")
else()
  target_compile_options(x3d2_test_utils PRIVATE "-cpp")
endif()

# Build the backend-dependent test runtime once per backend. A single build can
# contain CPU, CUDA, and OpenMP-target tests, so each variant needs its own 
# preprocessor definitions and module directory.
function(define_test_backend_runtime backend)
  set(runtime_target x3d2_test_runtime_${backend})
  if(TARGET ${runtime_target})
    return()
  endif()

  add_library(${runtime_target} STATIC common/backend_runtime.f90)
  set_target_properties(${runtime_target} PROPERTIES
    Fortran_MODULE_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/modules/${runtime_target}")
  target_include_directories(${runtime_target} INTERFACE
    "${CMAKE_CURRENT_BINARY_DIR}/modules/${runtime_target}")

  if(${backend} STREQUAL cuda)
    target_compile_options(${runtime_target} PRIVATE "-cpp" "-cuda" "-DCUDA")
    target_link_options(${runtime_target} INTERFACE "-cuda")
  else()
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
      target_compile_options(${runtime_target} PRIVATE "-eF")
    else()
      target_compile_options(${runtime_target} PRIVATE "-cpp")
    endif()
    if(${backend} STREQUAL omp_tgt)
      target_compile_options(${runtime_target} PRIVATE "-DOMP_TGT")
      if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" OR
         ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
        target_compile_options(${runtime_target} PRIVATE "-mp=gpu")
        target_link_options(${runtime_target} INTERFACE "-mp=gpu")
      else()
        target_compile_options(${runtime_target} PRIVATE "-fopenmp")
        target_link_options(${runtime_target} INTERFACE "-fopenmp")
      endif()
      # ifx takes its offload flags from OpenMP::OpenMP_Fortran instead
      if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
        target_compile_options(${runtime_target} PRIVATE
          "--offload-arch=${BACKEND_ARCH}")
        target_link_options(${runtime_target} INTERFACE
	  "--offload-arch=${BACKEND_ARCH}")
      endif()
    endif()
  endif()

  target_link_libraries(${runtime_target} PRIVATE x3d2_backends x3d2)
  # Match the OpenMP flags of the executables using this module (see
  # x3d2_test_utils above).  The CUDA tests build without OpenMP, so this
  # runtime does too.
  if(NOT ${backend} STREQUAL cuda)
    target_link_libraries(${runtime_target} PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  # backend_runtime.f90 reaches MPI through m_mpi (see x3d2_test_utils above).
  if(WITH_MPI)
    target_link_libraries(${runtime_target} PRIVATE MPI::MPI_Fortran)
  endif()
endfunction()

# Base test definition function
function(define_test testfile np backend)
  # Nothing can launch more than one rank without MPI.
  if(NOT WITH_MPI AND NOT ${np} EQUAL 1)
    return()
  endif()

  # Capture optional arguments passed to the function (e.g. input file or flags)
  set(extra_args "${ARGN}")
  string(REPLACE ";" " " extra_args_str "${extra_args}")

  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np}
    "${X3D2_TEST_NAME_SUFFIX}")

  define_test_backend_runtime(${backend})
  add_executable(${test_name} ${testfile})
  set_target_properties(${test_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/${test_name}")

  if(${backend} STREQUAL cuda)
    target_compile_options(${test_name} PRIVATE "-cpp")
    target_compile_options(${test_name} PRIVATE "-cuda")
    target_compile_options(${test_name} PRIVATE "-DCUDA")
    target_link_options(${test_name} INTERFACE "-cuda")
  else()
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
      target_compile_options(${test_name} PRIVATE "-eF")
      target_link_options(${test_name} PRIVATE -h omp)
      # Some CCE version define its internal $host_init$<proc>$<module>_ record both in
      # the object of the module owning the procedure and in every object that
      # uses it, and ld rejects the pair.  The records are zero-initialised .bss
      # data, so letting ld keep one is benign.
      target_link_options(${test_name} PRIVATE
        "-Wl,--allow-multiple-definition")
    else()
      target_compile_options(${test_name} PRIVATE "-cpp")
      if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        target_compile_options(${test_name} PRIVATE "-fallow-argument-mismatch")
      endif()
      if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
          ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wno-ambiguous-structure-constructor -Werror")
      endif()
    endif()
    # OpenMP is resolved once at the top level (see the root CMakeLists.txt).
    target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)

    if(${backend} STREQUAL omp_tgt)
      target_compile_options(${test_name} PRIVATE "-DOMP_TGT")
      if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
        target_compile_options(${test_name} PRIVATE "-fopenmp")
        target_link_options(${test_name} PRIVATE "-fopenmp")
      endif()
      if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" AND
         NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
        target_compile_options(${test_name} PRIVATE "--offload-arch=${BACKEND_ARCH}")
        target_link_options(${test_name} PRIVATE "--offload-arch=${BACKEND_ARCH}")
      endif()
    endif()
  endif()
  target_link_libraries(${test_name} PRIVATE x3d2_test_runtime_${backend})
  target_link_libraries(${test_name} PRIVATE x3d2_backends)
  target_link_libraries(${test_name} PRIVATE x3d2_test_utils)
  target_link_libraries(${test_name} PRIVATE x3d2)

  if(WITH_2DECOMPFFT)
    target_link_libraries(${test_name} PRIVATE decomp2d)
  endif()

  if(WITH_MPI)
    target_link_libraries(${test_name} PRIVATE MPI::MPI_Fortran)
  endif()

  if(WITH_ADIOS2)
    target_link_libraries(${test_name} PRIVATE ${ADIOS2_FORTRAN_TARGET})
  endif()

  if(NOT WITH_MPI)
    add_test(NAME ${test_name} COMMAND sh -c "bin/${test_name} ${extra_args_str}")
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    # A quoted COMMAND is one argv[0], so ctest looks for a program literally
    # named "srun bin/<test> ..." and reports every test as Not Run.  Go through
    # a shell, as the other two branches do, and pass the rank count on: Cray
    # machines are Slurm-driven and srun otherwise takes it from the allocation.
    add_test(NAME ${test_name}
      COMMAND sh -c "srun -n ${np} bin/${test_name} ${extra_args_str}")
  else()
    add_test(NAME ${test_name} COMMAND sh -c "mpirun ${MPIEXEC_OVERSUBSCRIBE_FLAG} -np ${np} bin/${test_name} ${extra_args_str}")
  endif()

  # default label includes backend
  set_tests_properties(${test_name} PROPERTIES LABELS "unit;${backend}")

  # The omp_tgt backend hands omp_target_alloc buffers to MPI. UCX tells
  # device from host buffers with a memory type cache fed by intercepting
  # CUDA allocations, which misses the OpenMP runtime's, so it reads them as
  # host memory and crashes in its CMA transport. Without the cache it asks
  # the driver about every buffer. MPIs not built on UCX ignore the variable.
  if(WITH_MPI AND ${backend} STREQUAL omp_tgt)
    set_property(TEST ${test_name} APPEND PROPERTY
      ENVIRONMENT "UCX_MEMTYPE_CACHE=n")
    # Intel MPI reads every buffer as host memory unless told to look for
    # device ones, and then faults in MPI_Isend on the first halo exchange.
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
      set_property(TEST ${test_name} APPEND PROPERTY
        ENVIRONMENT "I_MPI_OFFLOAD=1")
    endif()
  endif()
endfunction()

# Specialised test definition functions with category labels
function(define_verification_test testfile np backend)
  # define_test() registers nothing for a multi-rank test in a serial build.
  if(NOT WITH_MPI AND NOT ${np} EQUAL 1)
    return()
  endif()

  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np}
    "${X3D2_TEST_NAME_SUFFIX}")
  set_tests_properties(${test_name} PROPERTIES LABELS "verification;${backend}")
endfunction()

# Register a test for an executable that another target builds.
function(define_command_test test_name np labels pass_regex exe)
  # Nothing can launch more than one rank without MPI.
  if(NOT WITH_MPI AND NOT ${np} EQUAL 1)
    return()
  endif()

  set(extra_args "${ARGN}")
  string(REPLACE ";" " " extra_args_str "${extra_args}")

  # A test with a pass_regex needs a clean exit and the fixed string in the
  # output.  PASS_REGULAR_EXPRESSION alone accepts a failed run.
  set(wrapper "")
  if(pass_regex)
    set(wrapper "${CMAKE_CURRENT_SOURCE_DIR}/run_expect_output.sh '${pass_regex}' ")
  endif()

  if(NOT WITH_MPI)
    add_test(NAME ${test_name} COMMAND sh -c "${wrapper}${exe} ${extra_args_str}")
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    add_test(NAME ${test_name}
      COMMAND sh -c "${wrapper}srun -n ${np} ${exe} ${extra_args_str}")
  else()
    add_test(NAME ${test_name} COMMAND sh -c "${wrapper}mpirun --oversubscribe -np ${np} ${exe} ${extra_args_str}")
  endif()

  set_tests_properties(${test_name} PROPERTIES LABELS "${labels}")
endfunction()

function(define_performance_test testfile np backend)
  # define_test() registers nothing for a multi-rank test in a serial build.
  if(NOT WITH_MPI AND NOT ${np} EQUAL 1)
    return()
  endif()

  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "performance;${backend}")
endfunction()

function(define_regression_test testfile np backend)
  # define_test() registers nothing for a multi-rank test in a serial build.
  if(NOT WITH_MPI AND NOT ${np} EQUAL 1)
    return()
  endif()

  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "regression;${backend}")
endfunction()

# Full-case regression: run the xcompact solver on a small input file and diff
# the resulting monitoring.csv against a committed golden file. To (re)generate
# a golden after an intentional numerics change, run the case once and pass
# --bless to compare_monitoring.py (see that script's docstring). Golden files
# are only reproducible for a fixed compiler/precision/rank-count, so bless them
# in the same environment CI uses.
find_package(Python3 COMPONENTS Interpreter REQUIRED)

function(define_regression_case_test case_name input_file golden_file np backend)
  # Trailing arguments are passed through to the comparator (e.g. --ceiling).
  set(extra_args "${ARGN}")
  string(REPLACE ";" " " extra_args_str "${extra_args}")

  set(test_name regression_${case_name}_${backend}_${np})
  set(work_dir ${CMAKE_CURRENT_BINARY_DIR}/regression/${test_name})
  file(MAKE_DIRECTORY ${work_dir})
  set(comparator ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_monitoring.py)
  set(input ${CMAKE_CURRENT_SOURCE_DIR}/regression/${input_file})
  set(golden ${CMAKE_CURRENT_SOURCE_DIR}/regression/${golden_file})
  # Match the launcher define_test() uses: Cray machines are Slurm-driven and
  # have no mpirun, so an mpirun-based test cannot start there at all.
  if(NOT WITH_MPI)
    if(NOT ${np} EQUAL 1)
      return()
    endif()
    set(launcher "")
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    set(launcher "srun -n ${np}")
  else()
    set(launcher "mpirun ${MPIEXEC_OVERSUBSCRIBE_FLAG} -np ${np}")
  endif()
  add_test(
    NAME ${test_name}
    COMMAND sh -c
      "${launcher} $<TARGET_FILE:xcompact> ${input} && \
       ${Python3_EXECUTABLE} ${comparator} golden \
       --produced monitoring.csv --golden ${golden} ${extra_args_str}"
    WORKING_DIRECTORY ${work_dir})
  set_tests_properties(${test_name} PROPERTIES LABELS "regression;${backend}")
endfunction()

# Rank counts. Tests that decompose the domain over every rank run on np
# (CMAKE_CTEST_NPROCS) and must pass on any count. Tests built on a single-rank
# ([1, 1, 1]) decomposition, or on host code with no decomposition at all, run
# serially: extra ranks would only repeat the same work.
set(np ${CMAKE_CTEST_NPROCS})
set(serial 1)

# Tests of host-side code that never touches a backend: built once, under the
# CPU OpenMP backend runtime.
define_test(unit/test_allocator.f90 ${serial} omp)
define_test(unit/test_reorder_map.f90 ${serial} omp)
define_test(unit/test_poisson_bc_code.f90 ${serial} omp)
# Checks the exact per-rank split of a fixed grid, written for 4 ranks.
define_test(unit/test_mesh.f90 4 omp)
define_test(unit/test_statistics.f90 ${serial} omp)
define_test(unit/test_output_scalars.f90 ${serial} omp)
define_test(unit/test_les.f90 ${serial} omp)
define_test(unit/test_abl_diagnostics.f90 ${serial} omp)
define_test(unit/test_memcheck_stretched_y.f90 ${serial} omp)
define_test(unit/test_gpu_io_staging.f90 ${serial} omp)
define_test(unit/test_logger.f90 ${serial} omp)
# A serial build sets np to 1, which would register the same test twice.
if(NOT np EQUAL serial)
  define_test(unit/test_logger.f90 ${np} omp)
endif()
define_test(unit/test_relative_l2_error.f90 ${serial} omp)
if(NOT np EQUAL serial)
  define_test(unit/test_relative_l2_error.f90 ${np} omp)
endif()
# abort_run must stop every rank and report the error, whether one rank fails
# or all of them do. ctest checks the output rather than the exit code: a hang
# times out, and a run that returns from abort_run prints the FAIL line.
foreach(mode one_rank collective)
  set(X3D2_TEST_NAME_SUFFIX _${mode})
  define_test(unit/test_logger_abort.f90 ${np} omp ${mode})
  set_tests_properties(test_logger_abort_omp_${np}_${mode} PROPERTIES
    PASS_REGULAR_EXPRESSION "ERROR: deliberate"
    FAIL_REGULAR_EXPRESSION "must never be reached"
    TIMEOUT 60)
endforeach()
unset(X3D2_TEST_NAME_SUFFIX)
# A Poisson case type that leaves the periodicity reorder to the default of
# poisson_fft_t must stop, not fall through. Host code, so built once.
foreach(mode enforce undo)
  set(X3D2_TEST_NAME_SUFFIX _${mode})
  define_test(unit/test_poisson_fft_abort.f90 ${serial} omp ${mode})
  set_tests_properties(test_poisson_fft_abort_omp_${serial}_${mode} PROPERTIES
    PASS_REGULAR_EXPRESSION "ERROR: ${mode}_periodicity is not used"
    FAIL_REGULAR_EXPRESSION "must never be reached"
    TIMEOUT 60)
endforeach()
unset(X3D2_TEST_NAME_SUFFIX)
# relative_l2_error must stop on a zero reference or a non-finite input. ctest
# checks the output, not the exit code: MPI_Abort's depends on the launcher.
foreach(mode zero_ref nonfinite)
  set(X3D2_TEST_NAME_SUFFIX _${mode})
  define_test(unit/test_relative_l2_error.f90 ${np} omp ${mode})
  set_tests_properties(test_relative_l2_error_omp_${np}_${mode} PROPERTIES
    PASS_REGULAR_EXPRESSION "ERROR: relative_l2_error"
    FAIL_REGULAR_EXPRESSION "must never be reached"
    TIMEOUT 60)
endforeach()
unset(X3D2_TEST_NAME_SUFFIX)
if(WITH_ADIOS2)
  define_test(unit/test_adios2_mixed_precision.f90 ${serial} omp)
  define_test(unit/test_adm_checkpoint.f90 ${serial} omp)
  define_test(unit/test_snapshot_append.f90 ${np} omp)
endif()

# Backend tests. Every test below is written against base_backend_t through
# the m_backend_runtime helper (tests/common/backend_runtime.f90), so the same
# source is registered once per backend in the build: always the CPU OpenMP
# backend, plus CUDA or OpenMP-target when ENABLE_BACKEND selects one.
set(X3D2_TEST_BACKENDS omp)
if(${ENABLE_BACKEND} STREQUAL "CUDA")
  list(APPEND X3D2_TEST_BACKENDS cuda)
elseif(${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  list(APPEND X3D2_TEST_BACKENDS omp_tgt)
endif()

function(define_backend_tests backend)
  set(is_omp FALSE)
  set(is_device TRUE)
  if(${backend} STREQUAL omp)
    set(is_omp TRUE)
    set(is_device FALSE)
  endif()

  # Unit tests (fast, isolated component tests)
  define_test(unit/test_reordering.f90 ${serial} ${backend})
  define_test(unit/test_vecadd.f90 ${serial} ${backend})
  define_test(unit/test_setget_field.f90 ${serial} ${backend})
  define_test(unit/test_kinetic_energy.f90 ${np} ${backend})
  define_test(unit/test_snapshot_species_fields.f90 ${serial} ${backend})
  define_test(unit/test_scalar_product.f90 ${np} ${backend})
  define_test(unit/test_sum_intox.f90 ${np} ${backend})
  define_test(unit/test_boundary_planes.f90 ${serial} ${backend})
  define_test(unit/test_abl_seed.f90 ${serial} ${backend})
  define_test(unit/test_precorrect_walls.f90 ${serial} ${backend})
  define_test(unit/test_adm.f90 ${serial} ${backend})
  if(is_omp AND NOT np EQUAL serial)
    # Reordering is checked both without and with a decomposition.
    define_test(unit/test_reordering.f90 ${np} ${backend})
  endif()
  if(is_device)
    # Checks the allocator on a device selected per rank.
    define_test(unit/test_device_allocator.f90 ${serial} ${backend})
  endif()

  # Verification tests (numerical correctness vs analytical solutions)
  define_verification_test(verification/test_time_integrator.f90 ${serial} ${backend})
  define_verification_test(verification/test_tridiag.f90 ${np} ${backend})
  # tol_hyper must reject zero and sign-flipped output. The expected failure
  # gives a non-zero exit, so ctest checks the output, not the exit code.
  foreach(mode zero_output sign_flip)
    set(X3D2_TEST_NAME_SUFFIX _${mode})
    define_verification_test(verification/test_tridiag.f90 ${np} ${backend} ${mode})
    set_tests_properties(test_tridiag_${backend}_${np}_${mode} PROPERTIES
      PASS_REGULAR_EXPRESSION "Check 2nd ders, hyperviscous, dir-neu\\.\\.\\. failed"
      FAIL_REGULAR_EXPRESSION "Check 2nd ders, hyperviscous, dir-neu\\.\\.\\. passed"
      TIMEOUT 60)
  endforeach()
  unset(X3D2_TEST_NAME_SUFFIX)
  define_verification_test(verification/test_penta.f90 ${np} ${backend})
  define_verification_test(verification/test_transeq.f90 ${np} ${backend})
  define_verification_test(verification/test_les_field.f90 ${serial} ${backend})
  define_verification_test(verification/test_thom.f90 ${serial} ${backend})
  define_verification_test(verification/test_filter.f90 ${serial} ${backend})

  # The FFT Poisson solvers (2DECOMP&FFT on CPU, cuFFTMp on CUDA) are not at
  # feature parity, so the CUDA FFT test and the Poisson BC tests keep the rank
  # counts and BC configurations each backend has been validated with: several
  # BC combinations are not yet supported, or not yet accurate, on other counts.
  if(is_omp AND WITH_2DECOMPFFT)
    define_verification_test(verification/test_fft.f90 ${np} ${backend})
  elseif(${backend} STREQUAL cuda)
    define_verification_test(verification/test_fft.f90 1 ${backend})
  endif()
  if((is_omp AND WITH_2DECOMPFFT) OR ${backend} STREQUAL cuda)
    foreach(bc ppp pdp dpp ddp)
      set(X3D2_TEST_NAME_SUFFIX _${bc})
      define_verification_test(verification/test_poisson_bc.f90 1 ${backend} --bc ${bc})
    endforeach()
    # Validated 2 rank runs of the same test: the suffix keeps their names apart
    set(X3D2_TEST_NAME_SUFFIX _dpp)
    if(is_omp)
      define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc dpp --nproc_dir 1,2,1)
      set(X3D2_TEST_NAME_SUFFIX _pdp)
      define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc pdp)
    else()
      define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc dpp)
    endif()
    set(X3D2_TEST_NAME_SUFFIX _ddp)
    define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc ddp)
    # The backend's factory must stop on a BC code it has no case type for
    set(X3D2_TEST_NAME_SUFFIX _unsupported_bc)
    define_test(unit/test_poisson_fft_abort.f90 ${serial} ${backend} unsupported_bc)
    set_tests_properties(
      test_poisson_fft_abort_${backend}_${serial}_unsupported_bc PROPERTIES
      PASS_REGULAR_EXPRESSION "ERROR: Requested BCs are not supported"
      FAIL_REGULAR_EXPRESSION "must never be reached"
      TIMEOUT 60)
    unset(X3D2_TEST_NAME_SUFFIX)
    # The driver has its backend fixed at compile time, so register it for
    # the CUDA build only here; the CPU OMP registration is further down.
    if(${backend} STREQUAL cuda)
      define_command_test(poisson_solver_standalone_000_cuda 2 "unit;standalone" "Solve time (min):" ../bin/x3d2_poisson_solver 128 64 32)
    endif()
  endif()

  # I/O tests. A device backend needs an ADIOS2 that can stage device data.
  if(WITH_ADIOS2 AND (is_omp OR ADIOS2_WITH_CUDA))
    define_test(unit/test_adios2_read_write.f90 ${np} ${backend})
    define_test(unit/test_io_session.f90 ${np} ${backend})
    define_test(unit/test_ab_checkpoint.f90 ${serial} ${backend})
  endif()
  if(WITH_ADIOS2 AND (is_omp OR ${backend} STREQUAL cuda))
    # Writes on one rank layout and reads back on another: needs exactly 4.
    define_test(unit/test_inflow_planes.f90 4 ${backend})
  endif()
  if(${backend} STREQUAL cuda AND ADIOS2_WITH_CUDA AND WITH_ADIOS2_GPU_AWARE)
    define_test(unit/test_cuda_gpu_aware_io.f90 ${np} ${backend})
    define_test(unit/test_cuda_gpu_aware_io_padded.f90 ${np} ${backend})
  endif()

  # Performance tests (benchmarks, excluded from CI)
  define_performance_test(performance/perf_reorder.f90 ${serial} ${backend})
  define_performance_test(performance/perf_thom.f90 ${serial} ${backend})
  define_performance_test(performance/perf_tridiag.f90 ${np} ${backend})
  define_performance_test(performance/perf_penta.f90 ${np} ${backend})
  define_performance_test(performance/perf_transeq.f90 ${np} ${backend})
endfunction()

foreach(backend IN LISTS X3D2_TEST_BACKENDS)
  define_backend_tests(${backend})
endforeach()

# CUDA-only x3d2-memcheck tests, registered after the per-backend tests.
# x3d2-memcheck is not built without MPI (see src/CMakeLists.txt).
if(${ENABLE_BACKEND} STREQUAL "CUDA" AND WITH_MPI)
  # x3d2-memcheck regression tests: run the already-built binary against a
  # few example inputs and check its printed output, rather than compiling
  # a bespoke Fortran test program. "static" checks the tool runs and
  # prints a verdict with no full-grid cuFFTMp probe; "build" additionally
  # checks its real --build measurement completes and the static
  # estimate is close to it (the CHECK line) with no calibration drift (no
  # WARNING line) - this is the "is the guess close to the actual run"
  # check. PASS_REGULAR_EXPRESSION makes ctest ignore the process exit
  # code (x3d2-memcheck's own exit code is the verdict: 0/1/2 for
  # FITS/BORDERLINE/DOES_NOT_FIT, none of which should fail a correctly-
  # working test run). TGV/generic/channel are the three example inputs
  # that build unattended - cylinder needs an external IBM mask file not
  # in the repo, foil is not tracked here.
  function(define_memcheck_test case mode)
    set(name memcheck_${mode}_${case})
    set(input ${CMAKE_SOURCE_DIR}/examples/${case}/input.x3d)
    add_test(NAME ${name} COMMAND sh -c
      "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${input} --${mode}")
    set_tests_properties(${name} PROPERTIES LABELS "verification;cuda")
    if(${mode} STREQUAL build)
      set_tests_properties(${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "CHECK num_gpus=1: .* - OK"
        FAIL_REGULAR_EXPRESSION "WARNING: peak_fields_lookup|MISMATCH")
    else()
      set_tests_properties(${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "Requested nproc_dir gives num_gpus=")
    endif()
  endfunction()
  foreach(memcheck_case TGV generic channel)
    define_memcheck_test(${memcheck_case} static)
    define_memcheck_test(${memcheck_case} build)
  endforeach()

  # Every mode measures the CUDA context and the NVSHMEM heap on the GPU in
  # use, with the small-grid probe in the main process, and prints both on
  # one line. --static needs nothing else from the GPU.
  add_test(NAME memcheck_static_TGV_calibration COMMAND sh -c
    "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d --static")
  set_tests_properties(memcheck_static_TGV_calibration PROPERTIES
    LABELS "verification;cuda"
    PASS_REGULAR_EXPRESSION "Calibration on this GPU: CUDA context [0-9.]+ GiB, NVSHMEM heap [0-9.]+ GiB\\.")

  # Hiding every GPU from the process (CUDA_VISIBLE_DEVICES=) must be caught
  # cleanly by the device checks in memcheck.f90's main program, not crash.
  add_test(NAME memcheck_static_no_gpu COMMAND sh -c
    "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d --static")
  set_tests_properties(memcheck_static_no_gpu PROPERTIES
    LABELS "verification;cuda"
    ENVIRONMENT "CUDA_VISIBLE_DEVICES="
    PASS_REGULAR_EXPRESSION "no usable CUDA device")

  # An nproc_dir that does not divide the grid must be reported as
  # UNSUPPORTED (exit code 2), not silently estimated: TGV's z=256 does not
  # divide by 3.
  add_test(NAME memcheck_static_unsupported_ng COMMAND sh -c
    "d=$(mktemp -d); \
     sed 's/^nproc_dir *=.*/nproc_dir = 1, 1, 3/' \
       ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d > $d/input.x3d; \
     mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck \
       $d/input.x3d --static > $d/out.log 2>&1; \
     rc=$?; cat $d/out.log; \
     if [ $rc -eq 2 ] && \
        grep -q 'not supported (z=256 not divisible' $d/out.log; then \
       echo UNSUPPORTED_OK; \
     fi; \
     rm -rf $d")
  set_tests_properties(memcheck_static_unsupported_ng PROPERTIES
    LABELS "verification;cuda"
    PASS_REGULAR_EXPRESSION "UNSUPPORTED_OK")

  # peak_fields_lookup must follow the time integrator and the species
  # count: a --build of the TGV input with the edit applied must show no
  # calibration drift (no WARNING line) and a CHECK line that is OK. The
  # sed edit goes to a copy of the input in a mktemp directory. A sed
  # expression that leaves the input unchanged fails the test (SED_NOOP).
  # With the optional third argument the test passes on the token
  # <name>_OK printed by that shell snippet instead of on the CHECK line.
  # The snippet runs after the --build, with $d the directory, $rc the
  # exit code and $d/build.log the output of the --build.
  # With the optional fourth argument the run uses that option instead of
  # --build; an empty string runs the default mode.
  function(define_memcheck_build_edit_test name sed_expr)
    set(base ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d)
    set(pass_regex "CHECK num_gpus=1: .* - OK")
    set(snippet "true")
    set(mode_flag "--build")
    if(ARGC GREATER 3)
      set(mode_flag "${ARGV3}")
    endif()
    if(ARGC GREATER 2)
      set(pass_regex "${name}_OK")
      set(snippet "${ARGV2}")
    endif()
    add_test(NAME ${name} COMMAND sh -c
      "d=$(mktemp -d); \
       sed '${sed_expr}' ${base} > $d/input.x3d; \
       if cmp -s ${base} $d/input.x3d; then echo SED_NOOP; fi; \
       mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck \
         $d/input.x3d ${mode_flag} > $d/build.log 2>&1; \
       rc=$?; cat $d/build.log; \
       ${snippet}; \
       rm -rf $d")
    set_tests_properties(${name} PROPERTIES
      LABELS "verification;cuda"
      PASS_REGULAR_EXPRESSION "${pass_regex}"
      FAIL_REGULAR_EXPRESSION "SED_NOOP|WARNING: peak_fields_lookup|MISMATCH")
  endfunction()
  define_memcheck_build_edit_test(memcheck_build_TGV_rk4
    "s/^time_intg *=.*/time_intg = \"RK4\"/")
  define_memcheck_build_edit_test(memcheck_build_TGV_ab2
    "s/^time_intg *=.*/time_intg = \"AB2\"/")
  define_memcheck_build_edit_test(memcheck_build_TGV_species2
    "s/^time_intg *=.*/&\\nn_species = 2/")

  # A requested nproc_dir outside the GPU count list [1,2,4,8] must get the
  # same verdict in the initial report, in the measured table (its own row)
  # and in --static, with the same exit code. TGV's z=256 divides by 16.
  define_memcheck_build_edit_test(memcheck_build_ng16_matches_static
    "s/^nproc_dir *=.*/nproc_dir = 1, 1, 16/"
    "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck $d/input.x3d --static > $d/static.log 2>&1; rc_s=$?; v_s=$(sed -n 's/^Requested nproc_dir gives num_gpus=16: .* - //p' $d/static.log); v_i=$(sed -n 's/^Requested nproc_dir gives num_gpus=16: .* - //p' $d/build.log); v_r=$(sed -n 's/^ *16 |.*| //p' $d/build.log); if test x$v_s != x && test x$v_s = x$v_i && test x$v_s = x$v_r && test $rc_s -eq $rc; then echo memcheck_build_ng16_matches_static_OK; fi")

  # A requested nproc_dir that splits x must stay UNSUPPORTED after the real
  # build: exit code 2, never the FITS of the 1 GPU measurement.
  define_memcheck_build_edit_test(memcheck_build_x_split_exit2
    "s/^nproc_dir *=.*/nproc_dir = 2, 1, 1/"
    "if test $rc -eq 2 && grep -qF 'only nproc_dir = [1,1,num_gpus] is supported by the CUDA backend' $d/build.log; then echo memcheck_build_x_split_exit2_OK; fi")

  # A 2048^3 grid split over 64 GPUs fits the card per GPU, but the full-grid
  # probe on one GPU needs more than the free memory, so the tool runs no
  # full-grid probe: the plan would abort in NVSHMEM and end the tool. The
  # default mode must print the "no full-grid probe" note, mark the table
  # rows predicted, print the num_gpus=64 line, and exit with a verdict code
  # (0 to 2, exit 0 is legitimate as the grid fits per GPU) without an error
  # stop or an NVSHMEM abort. The calibration line is the only line that may
  # name NVSHMEM. On a card where every 2048^3/64 row exceeds 95% of the
  # memory (for example 16 GB) the note is suppressed, and the correct outcome
  # is exit 2 with DOES_NOT_FIT on the num_gpus=64 line, which the test also
  # accepts.
  define_memcheck_build_edit_test(memcheck_default_TGV_probe_skipped
    "s/^dims_global *=.*/dims_global = 2048, 2048, 2048/;s/^nproc_dir *=.*/nproc_dir = 1, 1, 64/"
    "if { test $rc -le 2 && grep -qF 'no full-grid probe: the plans need ' $d/build.log && grep -qF 'predicted rows follow' $d/build.log && grep -q '| predicted *| ' $d/build.log && grep -q 'Requested nproc_dir gives num_gpus=64' $d/build.log && ! grep -qi 'error stop' $d/build.log && ! grep -i nvshmem $d/build.log | grep -qv 'Calibration on this GPU'; } || { test $rc -eq 2 && grep -q 'Requested nproc_dir gives num_gpus=64: .* - DOES_NOT_FIT' $d/build.log && ! grep -qi 'error stop' $d/build.log && ! grep -i nvshmem $d/build.log | grep -qv 'Calibration on this GPU'; }; then echo memcheck_default_TGV_probe_skipped_OK; fi"
    "")

  # The report prints the CUDA context line once, and in the default mode it
  # marks each table row with how its overhead was obtained: the num_gpus=1 row
  # measured, the num_gpus>1 rows extrapolated by 1/num_gpus, with the note printed. The
  # sed edit only changes a comment line, so the input stays the TGV as is.
  define_memcheck_build_edit_test(memcheck_default_TGV_context_line_once
    "s/^! Flow case$/! Flow case, test copy/"
    "if test $(grep -c 'CUDA context (used at start)' $d/build.log) -eq 1; then echo memcheck_default_TGV_context_line_once_OK; fi"
    "")
  define_memcheck_build_edit_test(memcheck_default_TGV_extrapolation_note
    "s/^! Flow case$/! Flow case, test copy/"
    "if grep -qF '(num_gpus>1 rows extrapolate the single-GPU probe by 1/num_gpus.)' $d/build.log && grep -q '^ *1 | .*| measured *| ' $d/build.log && grep -q '^ *2 | .*| extrapolated *| ' $d/build.log; then echo memcheck_default_TGV_extrapolation_note_OK; fi"
    "")

  # A relative input path works: the real build re-reads argument 1 as the
  # input file, and the CHECK line stays OK.
  add_test(NAME memcheck_build_relative_input COMMAND sh -c
    "tmp=$(mktemp -d) && \
     cp ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d \"$tmp/input.x3d\" && \
     cd \"$tmp\" && \
     mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck \
       input.x3d --build; \
     rm -rf \"$tmp\"")
  set_tests_properties(memcheck_build_relative_input PROPERTIES
    LABELS "verification;cuda"
    PASS_REGULAR_EXPRESSION "CHECK num_gpus=1: .* - OK"
    FAIL_REGULAR_EXPRESSION "MISMATCH")

  # A real --build run writes no file in the working directory: the case
  # build writes no monitoring.csv. The directory holds only the input and
  # the log afterwards.
  add_test(NAME memcheck_build_leaves_no_files COMMAND sh -c
    "tmp=$(mktemp -d) && \
     cp ${CMAKE_SOURCE_DIR}/examples/TGV/input.x3d \"$tmp/input.x3d\" && \
     cd \"$tmp\" && \
     mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck \
       input.x3d --build > build.log 2>&1; \
     cat build.log; \
     if [ \"$(ls -A | sort | tr '\\n' ' ')\" = 'build.log input.x3d ' ]; then \
       echo NO_FILES_OK; else echo MISMATCH; ls -A; fi; \
     rm -rf \"$tmp\"")
  set_tests_properties(memcheck_build_leaves_no_files PROPERTIES
    LABELS "verification;cuda"
    PASS_REGULAR_EXPRESSION "NO_FILES_OK"
    FAIL_REGULAR_EXPRESSION "MISMATCH")

  # The query-only mode of the memcheck FFT workspace query reports the plan
  # memory and the workspace size without allocating the workspace.
  define_verification_test(verification/test_fft_query_only.f90 1 cuda)
  target_link_libraries(test_fft_query_only_cuda_1 PRIVATE x3d2_memcheck)

  # GPU-aware ADIOS2 device writes hold one extra unpadded local field on the
  # device while a snapshot/checkpoint field is packed (src/io/adios2/io.f90);
  # x3d2-memcheck adds that analytically. Only meaningful when the build has
  # the device path compiled in.
  if(WITH_ADIOS2_GPU_AWARE)
    set(memcheck_gpuio_input ${CMAKE_SOURCE_DIR}/tests/memcheck/tgv_snapshot_stride1.x3d)
    # One unpadded 256^3 field is 128 MiB in double precision, 64 MiB in
    # single precision.
    if(SINGLE_PREC)
      set(_gpuio_mib "64.0")
    else()
      set(_gpuio_mib "128.0")
    endif()
    add_test(NAME memcheck_gpuio_static_TGV COMMAND sh -c
      "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${memcheck_gpuio_input} --static")
    set_tests_properties(memcheck_gpuio_static_TGV PROPERTIES LABELS "verification;cuda"
      PASS_REGULAR_EXPRESSION "GPU-aware IO staging: ${_gpuio_mib} MiB/GPU")
    add_test(NAME memcheck_gpuio_static_TGV_hostmode COMMAND sh -c
      "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${memcheck_gpuio_input} --static")
    set_tests_properties(memcheck_gpuio_static_TGV_hostmode PROPERTIES LABELS "verification;cuda"
      ENVIRONMENT "X3D2_ADIOS2_GPU_WRITE_MODE=host"
      PASS_REGULAR_EXPRESSION "GPU-aware IO staging: none \\(write mode host\\)")
    # An unknown write mode must resolve as in the solver (auto), not print
    # the raw text.
    add_test(NAME memcheck_gpuio_static_TGV_unknown_mode COMMAND sh -c
      "mpirun --oversubscribe -np 1 ${PROJECT_BINARY_DIR}/bin/x3d2-memcheck ${memcheck_gpuio_input} --static")
    set_tests_properties(memcheck_gpuio_static_TGV_unknown_mode PROPERTIES LABELS "verification;cuda"
      ENVIRONMENT "X3D2_ADIOS2_GPU_WRITE_MODE=bogus"
      PASS_REGULAR_EXPRESSION "MiB/GPU at num_gpus=1 \\(write mode auto, "
      FAIL_REGULAR_EXPRESSION "write mode bogus")
  endif()
endif()

# Regression smoke test: fast full-solver TGV vs golden monitoring.csv.
# These reuse the main xcompact binary, whose backend follows the compiler
# (NVHPC/PGI -> CUDA). The committed goldens are CPU/OMP output, so only
# register them for a CPU OMP build; a CUDA golden would drive a separate case.
set(_cpu_omp_backend TRUE)
if(ENABLE_BACKEND)
  set(_cpu_omp_backend FALSE)
endif()

# div_u is a "stays at round-off" assertion, so it is bounded absolutely rather
# than diffed (see the comparator docstring). Round-off tracks the working
# precision, and the case, so the ceilings do too; rtol tracks precision
# because it sets what the suite can detect. The golden is double precision for
# both. The single-precision tolerances leave roughly an order of magnitude or
# more above the largest observed cross-build differences while retaining useful
# regression sensitivity. If a CI compiler drifts past these, the failure prints
# the observed |diff| so the tolerance can be revisited empirically.
if(SINGLE_PREC)
  set(_smoke_cmp_args --rtol 1e-5
      --ceiling div_u_max=1e-4 --ceiling div_u_mean=1e-5)
  set(_nightly_cmp_args --rtol 1e-4
      --ceiling div_u_max=1e-4 --ceiling div_u_mean=1e-5)
  # The channel's enstrophy drifts up to 3.6e-5 from the double precision
  # golden over its 50 noisy steps, past the smoke rtol.
  set(_channel_cmp_args --rtol 1e-4
      --ceiling div_u_max=1e-4 --ceiling div_u_mean=1e-5)
else()
  set(_smoke_cmp_args --rtol 1e-10
      --ceiling div_u_max=1e-8 --ceiling div_u_mean=1e-9)
  set(_nightly_cmp_args --rtol 1e-10
      --ceiling div_u_max=1e-12 --ceiling div_u_mean=1e-13)
  set(_channel_cmp_args ${_smoke_cmp_args})
endif()

if(WITH_2DECOMPFFT AND _cpu_omp_backend)
  define_regression_case_test(tgv_smoke tgv_smoke.x3d tgv_smoke_golden_omp.csv 2 omp
    ${_smoke_cmp_args})
endif()

# CUDA regression cases. The channel cases run the 010 FFT Poisson solver with
# a stretched y grid, one case per stretching operator. Only the CUDA backend
# supports y stretching in the FFT solver, and only on one rank. Their goldens
# are blessed from a double precision NVHPC build, and a single precision
# build checks against them too.
if(${ENABLE_BACKEND} STREQUAL "CUDA")
  # The CUDA backend matches the OMP golden to the smoke rtol, so the TGV case
  # shares it, which also checks that the two backends agree.
  define_regression_case_test(tgv_smoke tgv_smoke.x3d tgv_smoke_golden_omp.csv 2 cuda
    ${_smoke_cmp_args})
  define_regression_case_test(channel_smoke channel_smoke.x3d
    channel_smoke_golden_cuda.csv 1 cuda ${_channel_cmp_args})
  define_regression_case_test(channel_bottom_smoke channel_bottom_smoke.x3d
    channel_bottom_smoke_golden_cuda.csv 1 cuda ${_channel_cmp_args})
endif()

# The driver has its backend fixed at compile time, so register it as OMP only
# for a CPU OMP build.
if(WITH_2DECOMPFFT AND _cpu_omp_backend)
  define_command_test(poisson_solver_standalone_000_omp 2 "unit;standalone" "Solve time (min):" ../bin/x3d2_poisson_solver 128 64 32)
endif()

# Heavier nightly regression case (opt-in; kept out of PR CI to stay fast).
if(WITH_2DECOMPFFT AND _cpu_omp_backend AND WITH_NIGHTLY_REGRESSION)
  define_regression_case_test(tgv_nightly tgv_nightly.x3d tgv_nightly_golden_omp.csv 2 omp
    ${_nightly_cmp_args})
endif()
unset(_smoke_cmp_args)
unset(_nightly_cmp_args)
unset(_channel_cmp_args)
unset(_cpu_omp_backend)

set(EXECUTABLE_OUTPUT_PATH "bin")
