set(cmake_minimum_version 3.10)
cmake_minimum_required(VERSION ${cmake_minimum_version})
project(x3d2 LANGUAGES Fortran CXX C)
enable_testing()

set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)

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
else()
  message(STATUS "ADIOS2 is disabled")
endif()

set(OMP_TGT OFF CACHE BOOL
  "Enable OpenMP target offloading.")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
add_subdirectory(src)
add_subdirectory(tests)
