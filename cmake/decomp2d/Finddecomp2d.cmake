# - Find the 2decomp-fft library
find_package(decomp2d CONFIG
             PATHS ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt
             QUIET)
if (decomp2d_FOUND)
  message(STATUS "2decomp-fft FOUND")
else(decomp2d_FOUND)
  message(STATUS "2decomp-fft PATH not available we'll try to download and install")

  # Keep the driver project separate from the external project's BINARY_DIR.
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download")

  if (SINGLE_PREC)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/downloadBuild2decompSinglePrec.cmake.in decomp2d-download/CMakeLists.txt)
  else()
    configure_file(${CMAKE_CURRENT_LIST_DIR}/downloadBuild2decomp.cmake.in decomp2d-download/CMakeLists.txt)
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download )
  if(result)
      message(FATAL_ERROR "CMake step for 2decomp-fft failed: ${result}")
  else()
      message("CMake step for 2decomp-fft completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
         RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download )
  if(result)
      message(FATAL_ERROR "Build step for 2decomp-fft failed: ${result}")
  endif()
  set(D2D_ROOT ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt)

  find_package(decomp2d REQUIRED CONFIG
          PATHS ${D2D_ROOT}
          NO_DEFAULT_PATH)
endif(decomp2d_FOUND)

# Make the imported target global so sibling directories (e.g. tests/) can
# link it.
if(TARGET decomp2d)
  get_target_property(_decomp2d_is_global decomp2d IMPORTED_GLOBAL)
  if(NOT _decomp2d_is_global)
    set_target_properties(decomp2d PROPERTIES IMPORTED_GLOBAL TRUE)
  endif()
endif()
