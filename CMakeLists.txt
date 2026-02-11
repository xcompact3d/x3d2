set(cmake_minimum_version 3.10)
cmake_minimum_required(VERSION ${cmake_minimum_version})
project(x3d2 LANGUAGES Fortran CXX C)
enable_testing()

set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)
option(WITH_ADIOS2_GPU_AWARE "Enable GPU-aware I/O with ADIOS2 (requires X3D2_ADIOS2_CUDA=ON)" OFF)

if(SINGLE_PREC)
  add_definitions("-DSINGLE_PREC")
endif()

if(WITH_ADIOS2)
  # X3D2_ADIOS2_CUDA controls whether ADIOS2 is built with CUDA support
  # It can be set explicitly or auto-detected from compiler
  if(NOT DEFINED X3D2_ADIOS2_CUDA)
    set(X3D2_ADIOS2_CUDA OFF)
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
       ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
      set(X3D2_ADIOS2_CUDA ON)
    endif()
  endif()
  
  message(STATUS "X3D2_ADIOS2_CUDA is set to: ${X3D2_ADIOS2_CUDA}")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/adios2")
  cmake_policy(SET CMP0074 NEW)
  find_package(ADIOS2)
  
  if(WITH_ADIOS2_GPU_AWARE)
    if(NOT TARGET adios2::core_cuda)
      message(FATAL_ERROR
        "WITH_ADIOS2_GPU_AWARE=ON requires ADIOS2 built with CUDA support.\n"
        "Use the NVHPC/PGI compiler (auto-detected) or set -DX3D2_ADIOS2_CUDA=ON")
    endif()
    message(STATUS "ADIOS2 GPU-aware I/O enabled")
    add_definitions("-DX3D2_ADIOS2_CUDA")
  endif()
else()
  message(STATUS "ADIOS2 is disabled")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
add_subdirectory(src)
add_subdirectory(tests)
