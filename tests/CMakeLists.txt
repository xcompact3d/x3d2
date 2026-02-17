
set(CMAKE_CTEST_ARGUMENTS "--output-on-failure" CACHE STRING "comma separated list of arguments to pass to ctest")
set(CMAKE_CTEST_NPROCS "4" CACHE STRING "number of process used in tests")

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

  # default all tests defined via this function as 'unit' tests
  set_tests_properties(${test_name} PROPERTIES LABELS "unit")
endfunction()

function(define_performance_test testfile np backend)
  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "performance")
endfunction()

function(define_regression_test testfile np backend)
  define_test(${testfile} ${np} ${backend} ${ARGN})
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})
  set_tests_properties(${test_name} PROPERTIES LABELS "regression")
endfunction()

# omp backend shared tests
if(WITH_2DECOMPFFT)
  define_test(unit/test_fft.f90 4 omp)
endif()
define_test(unit/test_time_integrator.f90 1 omp)
define_test(unit/test_reordering.f90 1 omp)
define_test(unit/test_reordering.f90 2 omp)
define_test(unit/test_reordering.f90 4 omp)
define_test(unit/test_mesh.f90 4 omp)
define_test(unit/test_reorder_map.f90 1 omp)
define_test(unit/test_vecadd.f90 1 omp)
define_test(unit/test_scalar_product.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_sum_intox.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_setget_field.f90 1 omp)
define_test(unit/test_statistics.f90 1 omp)

# omp backend-specific tests
define_test(unit/test_omp_tridiag.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_omp_thom.f90 1 omp)
define_test(unit/test_omp_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_omp_transeq_species.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(unit/test_omp_dist_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)

if(WITH_ADIOS2)
  define_test(unit/test_adios2_read_write.f90 4 omp)
  define_test(unit/test_io_session.f90 4 omp)
  define_test(unit/test_ab_checkpoint.f90 1 omp)
endif()

# cuda backend shared tests
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  # shared tests
  define_test(unit/test_reordering.f90 1 cuda)
  define_test(unit/test_vecadd.f90 1 cuda)
  define_test(unit/test_time_integrator.f90 1 cuda)
  define_test(unit/test_fft.f90 1 cuda)
  define_test(unit/test_poisson_bc.f90 1 cuda)

  # cuda backend-specific tests
  define_test(unit/test_cuda_allocator.f90 1 cuda)
  define_test(unit/test_cuda_reorder.f90 1 cuda)
  define_test(unit/test_cuda_thom.f90 1 cuda)
  define_test(unit/test_cuda_tridiag.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_test(unit/test_cuda_transeq.f90 ${CMAKE_CTEST_NPROCS} cuda)

  if(WITH_ADIOS2)
    define_test(unit/test_adios2_read_write.f90 4 cuda)
    define_test(unit/test_io_session.f90 4 cuda)
    define_test(unit/test_ab_checkpoint.f90 1 cuda)
    if(WITH_ADIOS2_GPU_AWARE)
      define_test(unit/test_cuda_gpu_aware_io.f90 4 cuda)
    endif()
  endif()
endif()

# pFUnit integration 
find_package(pFUnit QUIET)

if(pFUnit_FOUND)
  # helper to define pFUnit tests
  function(define_pfunit_test test_source np)
    get_filename_component(test_name ${test_source} NAME_WE)
    
    # pFUnit macro to handle preprocessing .pf -> .f90
    add_pfunit_ctest(${test_name}
      TEST_SOURCE ${test_source}
      LINK_LIBRARIES x3d2
      TEST_NUM_PROCS ${np}
    )
    # label as unit test
    set_tests_properties(${test_name} PROPERTIES LABELS "unit;pfunit")
  endfunction()

  # register test
  define_pfunit_test(unit/test_allocator.pf 1)
  define_pfunit_test(unit/test_pfunit_demo.pf 1)
endif()

set(EXECUTABLE_OUTPUT_PATH "bin")
