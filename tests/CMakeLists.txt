set(TESTSRC
  test_allocator.f90
  test_tridiagonal.f90
)
set(CUDATESTSRC
  cuda/test_cuda_allocator.f90
)

foreach(testfile IN LISTS testsrc)
  get_filename_component(test_name ${testfile} NAME_WE)

  add_executable(${test_name} ${testfile})
  target_link_libraries(${test_name} PRIVATE x3d2)

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  foreach(testfile IN LISTS CUDATESTSRC)
  get_filename_component(test_name ${testfile} NAME_WE)

  add_executable(${test_name} ${testfile})
  target_link_libraries(${test_name} PRIVATE x3d2)
  target_compile_options(${test_name} "-cuda")

  add_test(NAME ${test_name} COMMAND ${test_name})
endif()
