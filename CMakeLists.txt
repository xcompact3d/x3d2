# Set mininum CMake requirement
set(cmake_minimum_version 3.13)
cmake_minimum_required(VERSION ${cmake_minimum_version})

# Define project languages
project(x3d2 LANGUAGES Fortran CXX C)

# Require OpenMP 
find_package(OpenMP REQUIRED)

# Toggleable build type (default Release)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose type of build" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")

# Single precision option (default OFF)
option(SINGLE_PREC "Enable single precision" OFF)
if(SINGLE_PREC)
  add_definitions("-DSINGLE_PREC")
endif()

# Backend for GPU offloading options (default OFF)
set(ENABLE_BACKEND OFF CACHE STRING "")
set_property(CACHE ENABLE_BACKEND PROPERTY STRINGS OFF CUDA OMP_TGT)

# ADIOS2 option (default OFF)
option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)

# FFT option (default on)
set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/decomp2d")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

if(WITH_ADIOS2)
  if (${ENABLE_BACKEND} STREQUAL "CUDA")
    set(ADIOS2_USE_CUDA ON)
  else()
    set(ADIOS2_USE_CUDA OFF)
  endif()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/adios2")
  cmake_policy(SET CMP0074 NEW)
  find_package(ADIOS2)
else()
  message(STATUS "ADIOS2 is disabled")
endif()

# Backend specific checks
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    message(FATAL_ERROR "CUDA backend require NVHPC and PGI compiler")
  endif()
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  if(OpenMP_VERSION VERSION_LESS "5.0")
    message(FATAL_ERROR "OpenMP version ${OpenMP_VERSION} is detected, but at least 5.0 is required with OMP_TGT.")
  endif()
  set(OMP_TGT_ARCH "" CACHE STRING "Target architecture for OpenMP offloading, e.g. gfx942")
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    if ("${OMP_TGT_ARCH}" STREQUAL "")
      message(FATAL_ERROR
        "OMP_TGT=ON requires setting OMP_TGT_ARCH=<arch> for this compiler")
    endif()
  endif()
 if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
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
endif()

# Give the imported target the link option FindOpenMP omits: every target that
# links OpenMP::OpenMP_Fortran picks it up, executables and tests alike.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND OpenMP_Fortran_FLAGS)
  set_property(TARGET OpenMP::OpenMP_Fortran APPEND PROPERTY
    INTERFACE_LINK_OPTIONS "SHELL:${OpenMP_Fortran_FLAGS}")
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
