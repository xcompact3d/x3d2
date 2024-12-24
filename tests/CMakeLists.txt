
set(CMAKE_CTEST_ARGUMENTS "--output-on-failure" CACHE STRING "comma separated list of arguments to pass to ctest")
set(CMAKE_CTEST_NPROCS "4" CACHE STRING "number of process used in tests")

function(define_test testfile np backend)
  get_filename_component(test_name ${testfile} NAME_WE)
  string(CONCAT test_name ${test_name} _ ${backend} _ ${np})

  add_executable(${test_name} ${testfile})

  if(${backend} STREQUAL cuda)
    target_compile_options(${test_name} PRIVATE "-cpp")
    target_compile_options(${test_name} PRIVATE "-cuda")
    target_compile_options(${test_name} PRIVATE "-DCUDA")
    target_link_options(${test_name} INTERFACE "-cuda")
  else()
    target_compile_options(${test_name} PRIVATE "-cpp")
    find_package(OpenMP REQUIRED)
    target_compile_options(${test_name} PRIVATE "-cpp")
    target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  target_link_libraries(${test_name} PRIVATE x3d2)

  add_test(NAME ${test_name} COMMAND sh -c "mpirun --oversubscribe -np ${np} ${test_name}")

endfunction()

define_test(test_allocator.f90 1 omp)
define_test(test_mesh.f90 4 omp)
define_test(test_reordering.f90 1 omp)
define_test(test_reordering.f90 2 omp)
define_test(test_reordering.f90 4 omp)
define_test(test_reorder_map.f90 1 omp)
define_test(test_time_integrator.f90 1 omp)
define_test(omp/test_omp_tridiag.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(omp/test_thom.f90 1 omp)
define_test(omp/test_omp_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(omp/test_omp_dist_transeq.f90 ${CMAKE_CTEST_NPROCS} omp)

define_test(test_scalar_product.f90 ${CMAKE_CTEST_NPROCS} omp)
define_test(test_sum_intox.f90 ${CMAKE_CTEST_NPROCS} omp)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
    ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  define_test(test_reordering.f90 1 cuda)
  define_test(test_time_integrator.f90 1 cuda)
  define_test(cuda/test_cuda_allocator.f90 1 cuda)
  define_test(cuda/test_cuda_reorder.f90 1 cuda)
  define_test(cuda/test_thom.f90 1 cuda)
  define_test(cuda/test_cuda_tridiag.f90 ${CMAKE_CTEST_NPROCS} cuda)
  define_test(cuda/test_cuda_transeq.f90 ${CMAKE_CTEST_NPROCS} cuda)
endif()

