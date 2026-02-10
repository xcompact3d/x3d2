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
  
  # Automatically enable GPU-aware I/O when using CUDA backend
  if(TARGET adios2::fortran_mpi AND X3D2_ADIOS2_CUDA)
    # Check if ADIOS2 has CUDA support
    get_target_property(ADIOS2_LINK_LIBS adios2::fortran_mpi INTERFACE_LINK_LIBRARIES)
    if(ADIOS2_LINK_LIBS)
      string(FIND "${ADIOS2_LINK_LIBS}" "CUDA" ADIOS2_HAS_CUDA)
      if(NOT ADIOS2_HAS_CUDA EQUAL -1)
        message(STATUS "ADIOS2 GPU-aware I/O enabled (CUDA backend detected)")
        add_definitions("-DADIOS2_GPU_AWARE")
      else()
        message(STATUS "ADIOS2 found but without CUDA support - using host-staged I/O")
      endif()
    else()
      # Fallback: assume GPU-aware is available when using CUDA compiler
      message(STATUS "ADIOS2 GPU-aware I/O enabled (CUDA compiler detected)")
      add_definitions("-DADIOS2_GPU_AWARE")
    endif()
  endif()
else()
  message(STATUS "ADIOS2 is disabled")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
add_subdirectory(src)
add_subdirectory(tests)
