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
      if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" AND
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
    add_test(NAME ${test_name} COMMAND sh -c "mpirun --oversubscribe -np ${np} bin/${test_name} ${extra_args_str}")
  endif()

  # default label includes backend
  set_tests_properties(${test_name} PROPERTIES LABELS "unit;${backend}")
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
    set(launcher "mpirun --oversubscribe -np ${np}")
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
# Checks the exact per-rank split of a fixed grid, written for 4 ranks.
define_test(unit/test_mesh.f90 4 omp)
define_test(unit/test_statistics.f90 ${serial} omp)
define_test(unit/test_output_scalars.f90 ${serial} omp)
define_test(unit/test_les.f90 ${serial} omp)
if(WITH_ADIOS2)
  define_test(unit/test_adios2_mixed_precision.f90 ${serial} omp)
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
  set(is_full_backend TRUE)
  if(${backend} STREQUAL omp)
    set(is_omp TRUE)
    set(is_device FALSE)
  elseif(${backend} STREQUAL omp_tgt)
    # The OpenMP-target backend implements only field copies, reorder,
    # vecadd/veccopy and vector_norm_squared so far; every other operation
    # stops with "not implemented". Tests needing the rest are registered
    # under is_full_backend below.
    set(is_full_backend FALSE)
  endif()

  # Unit tests (fast, isolated component tests)
  define_test(unit/test_reordering.f90 ${serial} ${backend})
  define_test(unit/test_vecadd.f90 ${serial} ${backend})
  define_test(unit/test_setget_field.f90 ${serial} ${backend})
  define_test(unit/test_kinetic_energy.f90 ${np} ${backend})
  define_test(unit/test_snapshot_species_fields.f90 ${serial} ${backend})
  if(is_full_backend)
    define_test(unit/test_scalar_product.f90 ${np} ${backend})
    define_test(unit/test_sum_intox.f90 ${np} ${backend})
  endif()
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
  if(is_full_backend)
    define_verification_test(verification/test_tridiag.f90 ${np} ${backend})
    define_verification_test(verification/test_penta.f90 ${np} ${backend})
    define_verification_test(verification/test_transeq.f90 ${np} ${backend})
    define_verification_test(verification/test_les_field.f90 ${serial} ${backend})
    define_verification_test(verification/test_thom.f90 ${serial} ${backend})
  endif()

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
      set(X3D2_TEST_NAME_SUFFIX _ddp)
      define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc ddp)
    else()
      define_verification_test(verification/test_poisson_bc.f90 2 ${backend} --bc dpp)
    endif()
    unset(X3D2_TEST_NAME_SUFFIX)
  endif()

  # I/O tests. A device backend needs an ADIOS2 that can stage device data.
  if(WITH_ADIOS2 AND (is_omp OR ADIOS2_WITH_CUDA))
    define_test(unit/test_adios2_read_write.f90 ${np} ${backend})
    define_test(unit/test_io_session.f90 ${np} ${backend})
    define_test(unit/test_ab_checkpoint.f90 ${serial} ${backend})
  endif()
  if(${backend} STREQUAL cuda AND ADIOS2_WITH_CUDA AND WITH_ADIOS2_GPU_AWARE)
    define_test(unit/test_cuda_gpu_aware_io.f90 ${np} ${backend})
    define_test(unit/test_cuda_gpu_aware_io_padded.f90 ${np} ${backend})
  endif()

  # Performance tests (benchmarks, excluded from CI)
  define_performance_test(performance/perf_reorder.f90 ${serial} ${backend})
  if(is_full_backend)
    define_performance_test(performance/perf_thom.f90 ${serial} ${backend})
    define_performance_test(performance/perf_tridiag.f90 ${np} ${backend})
    define_performance_test(performance/perf_penta.f90 ${np} ${backend})
    define_performance_test(performance/perf_transeq.f90 ${np} ${backend})
  endif()
endfunction()

foreach(backend IN LISTS X3D2_TEST_BACKENDS)
  define_backend_tests(${backend})
endforeach()

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
else()
  set(_smoke_cmp_args --rtol 1e-10
      --ceiling div_u_max=1e-8 --ceiling div_u_mean=1e-9)
  set(_nightly_cmp_args --rtol 1e-10
      --ceiling div_u_max=1e-12 --ceiling div_u_mean=1e-13)
endif()

if(WITH_2DECOMPFFT AND _cpu_omp_backend)
  define_regression_case_test(tgv_smoke tgv_smoke.x3d tgv_smoke_golden_omp.csv 2 omp
    ${_smoke_cmp_args})
endif()

# Heavier nightly regression case (opt-in; kept out of PR CI to stay fast).
if(WITH_2DECOMPFFT AND _cpu_omp_backend AND WITH_NIGHTLY_REGRESSION)
  define_regression_case_test(tgv_nightly tgv_nightly.x3d tgv_nightly_golden_omp.csv 2 omp
    ${_nightly_cmp_args})
endif()
unset(_smoke_cmp_args)
unset(_nightly_cmp_args)
unset(_cpu_omp_backend)

set(EXECUTABLE_OUTPUT_PATH "bin")
