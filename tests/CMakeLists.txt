
set(CMAKE_CTEST_ARGUMENTS "--output-on-failure" CACHE STRING "comma separated list of arguments to pass to ctest")
set(CMAKE_CTEST_NPROCS "4" CACHE STRING "number of process used in tests")

# Shared test utilities library
add_library(x3d2_test_utils STATIC common/test_utils.f90)
target_link_libraries(x3d2_test_utils PRIVATE x3d2)
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

  add_executable(${test_name} ${testfile})

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
    endif()
    find_package(OpenMP REQUIRED)
    target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  target_link_libraries(${test_name} PRIVATE x3d2)
  target_link_libraries(${test_name} PRIVATE x3d2_backends)
  target_link_libraries(${test_name} PRIVATE x3d2_test_utils)

if(WITH_2DECOMPFFT)
  find_package(decomp2d REQUIRED)
  include_directories(${decomp2d_INCLUDE_DIRS})
  target_link_libraries(decomp2d)
  target_link_libraries(${test_name} PRIVATE decomp2d)
  target_compile_options(x3d2 PRIVATE "-DWITH_2DECOMPFFT")
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

# Unit tests (fast, isolated component tests)
define_test(unit/test_allocator.f90 1 omp)
define_test(unit/test_reordering.f90 1 omp)
define_test(unit/test_reordering.f90 2 omp)
define_test(unit/test_reordering.f90 4 omp)
define_test(unit/test_mesh.f90 4 omp)
define_test(unit/test_reorder_map.f90 1 omp)
define_test(unit/test_vecadd.f90 1 omp)
define_test(unit/test_scalar_product.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_sum_intox.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_setget_field.f90 1 omp)
define_test(unit/test_snapshot_species_fields.f90 1 omp)
define_test(unit/test_statistics.f90 1 omp)

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
if(WITH_2DECOMPFFT)
  define_verification_test(verification/test_fft.f90 4 omp)
endif()

# Performance tests (benchmarks, excluded from CI)
define_performance_test(performance/perf_thom.f90 1 omp)
define_performance_test(performance/perf_omp_tridiag.f90 ${CMAKE_CTEST_NPROCS} omp)

# CUDA backend tests
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")

  # CUDA unit tests
  define_test(unit/test_reordering.f90 1 cuda)
  define_test(unit/test_vecadd.f90 1 cuda)
  define_test(unit/test_cuda_allocator.f90 1 cuda)
  define_test(unit/test_cuda_reorder.f90 1 cuda)

  # CUDA verification tests
  define_verification_test(verification/test_thom.f90 1 cuda)
  define_verification_test(verification/test_time_integrator.f90 1 cuda)
  define_verification_test(verification/test_fft.f90 1 cuda)
  define_verification_test(verification/test_poisson_bc.f90 1 cuda)
  define_verification_test(verification/test_cuda_tridiag.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_verification_test(verification/test_cuda_penta.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_verification_test(verification/test_cuda_transeq.f90 ${CMAKE_CTEST_NPROCS} cuda)

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
