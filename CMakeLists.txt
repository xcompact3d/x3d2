set(cmake_minimum_version 3.13)
cmake_minimum_required(VERSION ${cmake_minimum_version})
project(x3d2 LANGUAGES Fortran CXX C)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/decomp2d")

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
set(OMP_TGT_ARCH "" CACHE STRING "Target architecture for OpenMP offloading, e.g. gfx942")
if (${OMP_TGT} AND NOT
    ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  if ("${OMP_TGT_ARCH}" STREQUAL "")
    message(FATAL_ERROR
      "OMP_TGT=ON requires setting OMP_TGT_ARCH=<arch> for this compiler")
  endif()
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Resolve OpenMP once here so src/ and tests/ share a single import target.
find_package(OpenMP REQUIRED)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR
   CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  # FindOpenMP probing nvfortran with -cuda pulls the CUDA runtime into the OpenMP
  # target with cudadevrt not last -> nvlink "unexpected object after cudadevrt".
  # Drop them; the per-target -cuda link option relinks them with cudadevrt last.
  get_target_property(_omp_fortran_libs OpenMP::OpenMP_Fortran INTERFACE_LINK_LIBRARIES)
  if(_omp_fortran_libs)
    list(FILTER _omp_fortran_libs EXCLUDE REGEX "cuda")
    set_target_properties(OpenMP::OpenMP_Fortran PROPERTIES
      INTERFACE_LINK_LIBRARIES "${_omp_fortran_libs}")
  endif()
endif()

# Give the imported target the link option FindOpenMP omits: every target that
# links OpenMP::OpenMP_Fortran picks it up, executables and tests alike.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND OpenMP_Fortran_FLAGS)
  set_property(TARGET OpenMP::OpenMP_Fortran APPEND PROPERTY
    INTERFACE_LINK_OPTIONS "SHELL:${OpenMP_Fortran_FLAGS}")
endif()

add_subdirectory(src)
add_subdirectory(tests)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SelectBackend)
