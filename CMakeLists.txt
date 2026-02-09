set(cmake_minimum_version 3.10)
cmake_minimum_required(VERSION ${cmake_minimum_version})
project(x3d2 LANGUAGES Fortran CXX C)
enable_testing()

set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)
option(WITH_ADIOS2_GPU_AWARE "Enable ADIOS2 GPU-aware I/O (requires CUDA-enabled ADIOS2)" OFF)

if(SINGLE_PREC)
  add_definitions("-DSINGLE_PREC")
endif()

if(WITH_ADIOS2)
  set(X3D2_ADIOS2_CUDA OFF)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
    set(X3D2_ADIOS2_CUDA ON)
  endif()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/adios2")
  cmake_policy(SET CMP0074 NEW)
  find_package(ADIOS2)
  
  # check if ADIOS2 was built with CUDA support
  if(TARGET adios2::fortran_mpi AND X3D2_ADIOS2_CUDA AND WITH_ADIOS2_GPU_AWARE)
    # try to detect if ADIOS2 has CUDA support
    get_target_property(ADIOS2_INCLUDE_DIRS adios2::fortran_mpi INTERFACE_INCLUDE_DIRECTORIES)
    if(ADIOS2_INCLUDE_DIRS)
      message(STATUS "Checking ADIOS2 for CUDA support...")
      # check if ADIOS2 was compiled with CUDA
      # this is a heuristic check - we assume GPU-aware if CUDA compiler is used with ADIOS2
      message(STATUS "ADIOS2 GPU-aware I/O enabled (user requested)")
      add_definitions("-DADIOS2_GPU_AWARE")
    endif()
  endif()
else()
  message(STATUS "ADIOS2 is disabled")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
add_subdirectory(src)
add_subdirectory(tests)
