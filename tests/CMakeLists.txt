
set(CMAKE_CTEST_ARGUMENTS "--output-on-failure" CACHE STRING "comma separated list of arguments to pass to ctest")
set(CMAKE_CTEST_NPROCS "4" CACHE STRING "number of process used in tests")

option(WITH_NIGHTLY_REGRESSION
  "Register the heavier nightly regression cases (excluded from PR CI)." OFF)

# Shared test utilities library
add_library(x3d2_test_utils STATIC common/test_utils.f90)
target_link_libraries(x3d2_test_utils PRIVATE x3d2)
set_target_properties(x3d2_test_utils PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_test_utils")
target_include_directories(x3d2_test_utils INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_test_utils")
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  target_compile_options(x3d2_test_utils PRIVATE "-eF")
else()
  target_compile_options(x3d2_test_utils PRIVATE "-cpp")
endif()

# Base test definition function
function(define_test testfile np backend)
  # Capture optional arguments passed to the function (e.g. input file or flags)
  set(extra_args "${ARGN}")
  string(REPLACE ";" " " extra_args_str "${extra_args}")

  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})

  add_executable(${test_name} ${testfile}
                 ${CMAKE_SOURCE_DIR}/src/backend/backend_runtime.f90)
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
      target_link_options(${test_name} INTERFACE -h omp)
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
      # Note this is somewhat of a hack - hardcoded to build against MI300A
      target_compile_options(${test_name} PRIVATE "-fopenmp")
      target_compile_options(${test_name} PRIVATE "-DOMP_TGT")
      target_link_options(${test_name} PRIVATE "-fopenmp")
      if(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
        target_compile_options(${test_name} PRIVATE "--offload-arch=${OMP_TGT_ARCH}")
        target_link_options(${test_name} PRIVATE "--offload-arch=${OMP_TGT_ARCH}")
      endif()
    endif()
  endif()
  target_link_libraries(${test_name} PRIVATE x3d2_backends)
  target_link_libraries(${test_name} PRIVATE x3d2_test_utils)
  target_link_libraries(${test_name} PRIVATE x3d2)

  if(WITH_2DECOMPFFT)
    target_link_libraries(${test_name} PRIVATE decomp2d)
  endif()

  find_package(MPI REQUIRED)
  target_link_libraries(${test_name} PRIVATE MPI::MPI_Fortran)

  if(WITH_ADIOS2)
    target_link_libraries(${test_name} PRIVATE adios2::fortran_mpi)
  endif()

  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    add_test(NAME ${test_name} COMMAND "srun bin/${test_name} ${extra_args_str}")
  else()
    add_test(NAME ${test_name} COMMAND sh -c "mpirun --oversubscribe -np ${np} bin/${test_name} ${extra_args_str}")
  endif()

  # default label includes backend
  set_tests_properties(${test_name} PROPERTIES LABELS "unit;${backend}")
endfunction()

# Specialised test definition functions with category labels
function(define_verification_test testfile np backend)
  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "verification;${backend}")
endfunction()

function(define_performance_test testfile np backend)
  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "performance;${backend}")
endfunction()

function(define_regression_test testfile np backend)
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
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
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

# Unit tests (fast, isolated component tests)
define_test(unit/test_allocator.f90 1 omp)
define_test(unit/test_reordering.f90 1 omp)
define_test(unit/test_reordering.f90 2 omp)
define_test(unit/test_reordering.f90 4 omp)
define_test(unit/test_mesh.f90 4 omp)
define_test(unit/test_reorder_map.f90 1 omp)
define_test(unit/test_vecadd.f90 1 omp)
define_test(unit/test_scalar_product.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_kinetic_energy.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_sum_intox.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_setget_field.f90 1 omp)
define_test(unit/test_snapshot_species_fields.f90 1 omp)
define_test(unit/test_statistics.f90 1 omp)
define_test(unit/test_output_scalars.f90 1 omp)
define_test(unit/test_les.f90 1 omp)

if(WITH_ADIOS2)
  define_test(unit/test_adios2_read_write.f90 4 omp)
  define_test(unit/test_io_session.f90 4 omp)
  define_test(unit/test_ab_checkpoint.f90 1 omp)
endif()

# Verification tests (numerical correctness vs analytical solutions)
define_verification_test(verification/test_thom.f90 1 omp)
define_verification_test(verification/test_omp_tridiag.f90 ${CMAKE_CTEST_NPROCS} omp)
define_verification_test(verification/test_omp_penta.f90 ${CMAKE_CTEST_NPROCS} omp)
define_verification_test(verification/test_omp_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)
define_verification_test(verification/test_omp_transeq_species.f90 ${CMAKE_CTEST_NPROCS} omp)
define_verification_test(verification/test_omp_dist_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)
define_verification_test(verification/test_time_integrator.f90 1 omp)
define_verification_test(verification/test_les_field.f90 1 omp)
if(WITH_2DECOMPFFT)
  define_verification_test(verification/test_fft.f90 4 omp)
endif()

# Regression smoke test: fast full-solver TGV vs golden monitoring.csv.
# These reuse the main xcompact binary, whose backend follows the compiler
# (NVHPC/PGI -> CUDA). The committed goldens are CPU/OMP output, so only
# register them for a CPU OMP build; a CUDA golden would drive a separate case.
set(_cpu_omp_backend TRUE)
if((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR
    CMAKE_Fortran_COMPILER_ID STREQUAL "PGI") OR OMP_TGT)
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

if(${OMP_TGT})
  # OMP TGT unit tests
  define_test(unit/test_vecadd.f90 1 omp_tgt)
  define_test(unit/test_setget_field.f90 1 omp_tgt)
  define_test(unit/test_omptgt_allocator.f90 1 omp_tgt)

  # OMP TGT verification tests
  define_verification_test(verification/test_time_integrator.f90 1 omp_tgt)
endif()

# Performance tests (benchmarks, excluded from CI)
define_performance_test(performance/perf_thom.f90 1 omp)
define_performance_test(performance/perf_omp_tridiag.f90 ${CMAKE_CTEST_NPROCS} omp)

# CUDA backend tests
# Parenthesise the compiler-ID check so NOT OMP_TGT gates both PGI and NVHPC;
# without it CMake precedence (NOT > AND > OR) leaves an OpenMP-target build
# able to pick up CUDA tests.
if((${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
    ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC") AND
   NOT OMP_TGT)

  # CUDA unit tests
  define_test(unit/test_reordering.f90 1 cuda)
  define_test(unit/test_vecadd.f90 1 cuda)
  define_test(unit/test_kinetic_energy.f90 1 cuda)
  define_test(unit/test_setget_field.f90 1 cuda)
  define_test(unit/test_cuda_allocator.f90 1 cuda)
  define_test(unit/test_cuda_reorder.f90 1 cuda)

  # CUDA verification tests
  define_verification_test(verification/test_thom.f90 1 cuda)
  define_verification_test(verification/test_time_integrator.f90 1 cuda)
  define_verification_test(verification/test_fft.f90 1 cuda)
    define_verification_test(verification/test_poisson_bc.f90 1 cuda)
  define_verification_test(verification/test_poisson_bc.f90 2 cuda --config 100)
  define_verification_test(verification/test_cuda_tridiag.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_verification_test(verification/test_cuda_penta.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_verification_test(verification/test_cuda_transeq.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_verification_test(verification/test_les_field.f90 1 cuda)

  # CUDA performance tests
  define_performance_test(performance/perf_thom.f90 1 cuda)
  define_performance_test(performance/perf_cuda_reorder.f90 1 cuda)
  define_performance_test(performance/perf_cuda_tridiag.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_performance_test(performance/perf_cuda_penta.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_performance_test(performance/perf_cuda_transeq.f90 ${CMAKE_CTEST_NPROCS} cuda)

  if(WITH_ADIOS2)
    define_test(unit/test_adios2_read_write.f90 4 cuda)
    define_test(unit/test_io_session.f90 4 cuda)
    define_test(unit/test_ab_checkpoint.f90 1 cuda)
  endif()
endif()

set(EXECUTABLE_OUTPUT_PATH "bin")
