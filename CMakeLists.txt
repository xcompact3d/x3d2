# Set mininum CMake requirement
cmake_minimum_required(VERSION 3.13)

# Define project languages
project(x3d2 LANGUAGES Fortran CXX C)

# Require OpenMP 
find_package(OpenMP REQUIRED)

# MPI option (default ON).  A build without it is a single serial rank: see
# src/mpi.f90, which stands in for the MPI entities the code uses.
option(WITH_MPI "Enable MPI" ON)
if(WITH_MPI)
  find_package(MPI REQUIRED)
  add_definitions("-DMPI")
endif()

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

# FFT option (default ON)
option(WITH_2DECOMPFFT "Enable Poisson based FFT solver on the OpenMP backend" ON)
if(WITH_2DECOMPFFT AND NOT WITH_MPI)
  message(STATUS
    "2decomp-fft requires MPI: disabling WITH_2DECOMPFFT for this serial build")
  set(WITH_2DECOMPFFT OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Backend specific checks
if (${ENABLE_BACKEND} STREQUAL "CUDA" OR ${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  set(BACKEND_ARCH "" CACHE STRING "Target architecture for GPU offloading")
  if ("${BACKEND_ARCH}" STREQUAL "")
    message(FATAL_ERROR
      "CUDA and OMP_TGT backend requires setting BACKEND_ARCH=<arch>")
  endif()
endif()
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    message(FATAL_ERROR "CUDA backend require NVHPC and PGI compiler")
  endif()
  IF(NOT BACKEND_ARCH MATCHES "sm_")
    message(FATAL_ERROR "For Nvidia GPU, BACKEND_ARCH must be specified as sm_**")
  ENDIF()
  # nvfortran names compute capabilities ccXX, not sm_XX: -gpu=cc80 for sm_80.
  string(REPLACE "sm_" "cc" NVHPC_GPU_ARCH "${BACKEND_ARCH}")
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  if("${OpenMP_Fortran_VERSION}" VERSION_LESS "4.5")
    message(FATAL_ERROR "OpenMP Fortran version ${OpenMP_Fortran_VERSION} is detected, but at least 4.5 is required with OMP_TGT.")
  endif()
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    IF(NOT BACKEND_ARCH MATCHES "gfx*")
      message(FATAL_ERROR "For AMD GPU, BACKEND_ARCH must be specified as gfx***")
    ENDIF()
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    IF(NOT BACKEND_ARCH MATCHES "sm_")
      message(FATAL_ERROR "For Nvidia GPU, BACKEND_ARCH must be specified as sm_**")
    ENDIF()
    # nvfortran names compute capabilities ccXX, not sm_XX: -gpu=cc80 for sm_80.
    string(REPLACE "sm_" "cc" NVHPC_GPU_ARCH "${BACKEND_ARCH}")
    # FindOpenMP supplies -mp, which targets the CPU with NVHPC.  Compile and
    # link target regions for the GPU when this backend is selected.
    set_target_properties(OpenMP::OpenMP_Fortran PROPERTIES
      INTERFACE_COMPILE_OPTIONS "-mp=gpu;-gpu=${NVHPC_GPU_ARCH}"
      INTERFACE_LINK_OPTIONS "-mp=gpu;-gpu=${NVHPC_GPU_ARCH}")
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

# Flang links libomp.so out of its own installation, which is not on the
# loader's default search path, and neither the driver nor FindOpenMP (it
# reports only -fopenmp here) leaves a library path behind for CMake to turn
# into an RPATH.  Binaries then link but fail to start.  Locate the runtime and
# put its directory on the RPATH instead, so running the build does not depend
# on LD_LIBRARY_PATH being set.
set(FLANG_OPENMP_RUNTIME_DIR "" CACHE PATH
  "Directory holding flang's libomp.so; detected automatically when empty")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang" OR
   CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
  set(_flang_libomp_dir "${FLANG_OPENMP_RUNTIME_DIR}")
  if(_flang_libomp_dir)
    if(NOT EXISTS "${_flang_libomp_dir}/libomp.so")
      message(FATAL_ERROR
        "FLANG_OPENMP_RUNTIME_DIR is set to ${_flang_libomp_dir}, but "
        "libomp.so is not there. Point it at the directory holding the flang "
        "OpenMP runtime, or leave it empty to detect it automatically.")
    endif()
  else()
    # Ask the driver first: it knows about resource and per-target directories
    # that cannot be guessed from the install prefix alone.  It answers with the
    # bare name when it cannot find the file, and some builds reject the option
    # outright, so the result has to be checked rather than trusted.
    execute_process(
      COMMAND "${CMAKE_Fortran_COMPILER}" -print-file-name=libomp.so
      OUTPUT_VARIABLE _flang_libomp
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(IS_ABSOLUTE "${_flang_libomp}" AND EXISTS "${_flang_libomp}")
      get_filename_component(_flang_libomp_dir "${_flang_libomp}" DIRECTORY)
    else()
      # Fall back to searching the toolchain the driver was invoked from, which
      # covers both a plain LLVM install (lib/, or lib/<triple>/ when OpenMP is
      # built as an LLVM_ENABLE_RUNTIMES project) and the ROCm layout.
      get_filename_component(_flang_root "${CMAKE_Fortran_COMPILER}" DIRECTORY)
      get_filename_component(_flang_root "${_flang_root}" DIRECTORY)
      file(GLOB _flang_libomp_hints "${_flang_root}/lib/*")
      find_library(FLANG_OPENMP_RUNTIME
        NAMES omp
        HINTS "${_flang_root}/lib" "${_flang_root}/lib64" ${_flang_libomp_hints}
        NO_DEFAULT_PATH
        DOC "Path to flang's OpenMP runtime library")
      if(FLANG_OPENMP_RUNTIME)
        get_filename_component(_flang_libomp_dir "${FLANG_OPENMP_RUNTIME}" DIRECTORY)
      endif()
      unset(_flang_libomp_hints)
      unset(_flang_root)
    endif()
    unset(_flang_libomp)
  endif()
  if(_flang_libomp_dir)
    set_property(TARGET OpenMP::OpenMP_Fortran APPEND PROPERTY
      INTERFACE_LINK_OPTIONS "LINKER:-rpath,${_flang_libomp_dir}")
    message(STATUS "Flang OpenMP runtime found in ${_flang_libomp_dir}")
  else()
    message(WARNING
      "Could not locate libomp.so for ${CMAKE_Fortran_COMPILER}. Point "
      "FLANG_OPENMP_RUNTIME_DIR at the directory holding it, or the built "
      "binaries will only run with that directory on LD_LIBRARY_PATH.")
  endif()
  unset(_flang_libomp_dir)
endif()

enable_testing()
include(ThirdPartyADIOS2)
include(ThirdPartyDecomp2D)
add_subdirectory(src)
add_subdirectory(tests)
