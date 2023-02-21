file(GLOB TESTSRC
  LIST_DIRECTORIES false
  *.f90
)
foreach(testsrc IN LISTS TESTSRC)
  get_filename_component(test_name ${testsrc} NAME_WE)
  add_executable(${test_name} ${testsrc})
  target_link_libraries(${test_name} PRIVATE x3d2)

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
