set(TESTSRC
  test_allocator.f90
  omp/test_omp_tridiag.f90
  omp/test_omp_transeq.f90
  omp/test_omp_dist_transeq.f90
)
set(CUDATESTSRC
  cuda/test_cuda_allocator.f90
  cuda/test_cuda_reorder.f90
  cuda/test_cuda_tridiag.f90
  cuda/test_cuda_transeq.f90
)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  foreach(testfile IN LISTS CUDATESTSRC)
    get_filename_component(test_name ${testfile} NAME_WE)

    add_executable(${test_name} ${testfile})
    target_compile_options(${test_name} PRIVATE "-cuda")
    target_compile_options(${test_name} PRIVATE "-O3")
    target_compile_options(${test_name} PRIVATE "-fast")
    target_link_options(${test_name} INTERFACE "-cuda")
    target_link_libraries(${test_name} PRIVATE x3d2)

    add_test(NAME ${test_name} COMMAND sh -c "mpirun -np 1 ${test_name}")
  endforeach()
endif()

list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

foreach(testfile IN LISTS TESTSRC)
  get_filename_component(test_name ${testfile} NAME_WE)

  add_executable(${test_name} ${testfile})

  target_compile_options(${test_name} PRIVATE "-O3")

  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(${test_name} PRIVATE "-ffast-math")
  endif()

  target_link_libraries(${test_name} PRIVATE x3d2)

  find_package(OpenMP REQUIRED)
  target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_Fortran)

  add_test(NAME ${test_name} COMMAND sh -c "mpirun -np 4 ${test_name}")
endforeach()
