# Set mininum CMake requirement
cmake_minimum_required(VERSION 3.13)

# Define project languages
project(x3d2 LANGUAGES Fortran CXX C)

# Require OpenMP 
find_package(OpenMP REQUIRED)

option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)

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

if(WITH_ADIOS2)
  # Allow command-line override, otherwise auto-detect based on compiler
  if(NOT DEFINED X3D2_ADIOS2_CUDA)
    set(X3D2_ADIOS2_CUDA OFF)
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
       ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
      set(X3D2_ADIOS2_CUDA ON)
    endif()
  endif()
  
  # Enable CUDA language if GPU-aware I/O is requested
  if(X3D2_ADIOS2_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES "80" CACHE STRING "CUDA architectures")
    endif()
    message(STATUS "CUDA Toolkit found: ${CUDAToolkit_VERSION}")
    message(STATUS "CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
  endif()
  
  message(STATUS "X3D2_ADIOS2_CUDA is set to: ${X3D2_ADIOS2_CUDA}")
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
        message(STATUS "ADIOS2 GPU-aware I/O enabled (CUDA backend detected in ADIOS2)")
        add_definitions("-DADIOS2_GPU_AWARE")
      else()
        message(STATUS "ADIOS2 found but without CUDA support - using host-staged I/O")
        message(STATUS "To enable GPU-aware I/O, rebuild ADIOS2 with -DADIOS2_USE_CUDA=ON")
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

# Backend specific checks
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    message(FATAL_ERROR "CUDA backend require NVHPC and PGI compiler")
  endif()
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  if("${OpenMP_Fortran_VERSION}" VERSION_LESS "4.5")
    message(FATAL_ERROR "OpenMP Fortran version ${OpenMP_Fortran_VERSION} is detected, but at least 4.5 is required with OMP_TGT.")
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
    # FindOpenMP supplies -mp, which targets the CPU with NVHPC.  Compile and
    # link target regions for the GPU when this backend is selected.
    set_target_properties(OpenMP::OpenMP_Fortran PROPERTIES
      INTERFACE_COMPILE_OPTIONS "-mp=gpu"
      INTERFACE_LINK_OPTIONS "-mp=gpu")
    # FindOpenMP can pull CUDA runtime libraries into the target in an order
    # that nvlink rejects.  Let -mp=gpu add them in the correct order.
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
include(ThirdPartyADIOS2)
include(ThirdPartyDecomp2D)
add_subdirectory(src)
add_subdirectory(tests)
