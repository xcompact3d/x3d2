# FindADIOS2.cmake
#
# This module tries to locate an ADIOS2 installation.
# It first searches for a local installation in a specified directory.
# If ADIOS2 is not found, it invokes an ExternalProject-based script
# (downloadBuildAdios2.cmake.in) to download, build, and install ADIOS2 locally.
#
# After the install, it re-invokes find_package to load the locally built ADIOS2.

if(DEFINED _ADIOS2_FIND_GUARD)
  return()
endif()
set(_ADIOS2_FIND_GUARD ON)

find_package(ADIOS2 QUIET
  PATHS "${CMAKE_BINARY_DIR}/adios2-build/adios2-install/lib/cmake/adios2"
)

if(ADIOS2_FOUND)
  message(STATUS "ADIOS2 found at: ${ADIOS2_DIR}")
else()
  message(STATUS "ADIOS2 not found locally; attempting auto-download and install")

  # copy the external project script into a subdirectory of the build tree
  configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/downloadBuildAdios2.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/adios2-build/CMakeLists.txt"
    COPYONLY
  )

  # ensure that the build directory for the external project exists
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/adios2-build")

  # run cmake to configure the ADIOS2 external project
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/adios2-build"
    RESULT_VARIABLE adios2_cmake_res
  )
  if(adios2_cmake_res)
    message(FATAL_ERROR "CMake configuration step for ADIOS2 failed with code: ${adios2_cmake_res}")
  else()
    message(STATUS "CMake configuration step for ADIOS2 completed (exit code: ${adios2_cmake_res}).")
  endif()

  # build & install ADIOS2 using the external project target
  # note: our external project script (downloadBuildAdios2.cmake.in) defines a target
  # named "downloadBuildAdios2" which performs the install step
  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build . --target downloadBuildAdios2
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/adios2-build"
    RESULT_VARIABLE adios2_build_res
  )
  if(adios2_build_res)
    message(FATAL_ERROR "Build step for ADIOS2 failed with code: ${adios2_build_res}")
  endif()

  # set the local ADIOS2 root to where our external project installed ADIOS2
  set(ADIOS2_ROOT "${CMAKE_CURRENT_BINARY_DIR}/adios2-build/adios2-install")

  # re-run find_package to load ADIOS2 from the local install
  find_package(ADIOS2 REQUIRED PATHS "${ADIOS2_ROOT}" NO_DEFAULT_PATH)
  if(ADIOS2_FOUND)
    message(STATUS "ADIOS2 was auto-downloaded and installed at: ${ADIOS2_DIR}")
  else()
    message(FATAL_ERROR "ADIOS2 not found after auto-build!")
  endif()
endif()
