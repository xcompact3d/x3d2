set(srctest_allocator ${CMAKE_SOURCE_DIR}/src/allocator.f90)
set(srctest_tridiagonal ${CMAKE_SOURCE_DIR}/src/thomas.f90)
# set(srcstencil
#   ${CMAKE_SOURCE_DIR}../src/stencil.f90
#   ${CMAKE_SOURCE_DIR}../src/stencil_definitions.f90
# )
set(testsrc
  ${CMAKE_CURRENT_SOURCE_DIR}/test_allocator.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/test_tridiagonal.f90
)
foreach(testfile IN LISTS testsrc)
  get_filename_component(test_name ${testfile} NAME_WE)
  add_executable(${test_name} ${testfile} ${src${test_name}})
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
