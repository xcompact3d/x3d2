# FindADIOS2.cmake
#
# This module locates and configures ADIOS2 for use in this project.
#
# If USE_SYSTEM_ADIOS2=ON:
#   - It searches for an existing ADIOS2 installation in the specified
#     ADIOS2_ROOT_DIR or in standard system paths.
#
# If USE_SYSTEM_ADIOS2=OFF (default):
#   - It downloads, builds, and installs ADIOS2 locally using the
#     ExternalProject module.
#
# After ADIOS2 is located or built, it makes the libraries and include
# files available to the project through the standard find_package mechanism.
#
# Options:
#   USE_SYSTEM_ADIOS2 - Set to ON to use a system-installed version of ADIOS2
#   ADIOS2_ROOT_DIR - Directory where ADIOS2 is installed (when USE_SYSTEM_ADIOS2=ON)
option(USE_SYSTEM_ADIOS2 "Use system-installed ADIOS2" OFF)
set(ADIOS2_ROOT_DIR "" CACHE PATH "Directory where ADIOS2 is installed (optional)")

if(NOT USE_SYSTEM_ADIOS2 AND ADIOS2_ROOT_DIR)
  message(WARNING "ADIOS2_ROOT_DIR is set but will be ignored because USE_SYSTEM_ADIOS2=OFF")
endif()

if(DEFINED _ADIOS2_FIND_GUARD)
  return()
endif()
set(_ADIOS2_FIND_GUARD ON)

if (USE_SYSTEM_ADIOS2)
  message(STATUS "Searching for system-installed ADIOS2...")

  if (ADIOS2_ROOT_DIR)
    message(STATUS "Looking for ADIOS2 in user-specified directory: ${ADIOS2_ROOT_DIR}")
    find_package(ADIOS2 QUIET PATHS "${ADIOS2_ROOT_DIR}" NO_DEFAULT_PATH)
  else()
    find_package(ADIOS2 QUIET
      PATHS "${CMAKE_BINARY_DIR}/adios2-build/adios2-install/lib/cmake/adios2"
    )
  endif()

  if (ADIOS2_FOUND)
    message(STATUS "ADIOS2 found at: ${ADIOS2_DIR}")
  else()
    message(FATAL_ERROR "USE_SYSTEM_ADIOS2 is ON but ADIOS2 was not found. "
    "Please install ADIOS2, specify ADIOS2_ROOT_DIR, "
    "or set USE_SYSTEM_ADIOS2=OFF to auto-download.")
  endif()
else()
  message(STATUS "USE_SYSTEM_ADIOS2 is OFF; building ADIOS2 from source...")

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
