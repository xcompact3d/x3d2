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
  # The version string tells the implementations apart below.
  set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  find_package(MPI REQUIRED)
  add_definitions("-DMPI")

  # Only Open MPI needs --oversubscribe, and only it accepts the flag: it
  # refuses to start more ranks than the node has slots, where the Hydra
  # launchers of Intel MPI and MPICH oversubscribe unasked and reject the
  # option outright.  The multi-rank tests ask for up to 4 ranks, so the
  # launcher in tests/CMakeLists.txt takes the flag from here.
  set(MPIEXEC_OVERSUBSCRIBE_FLAG "")
  if(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES "Open MPI")
    set(MPIEXEC_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()

  # An implementation is free to leave routines out of the `mpi` module: Intel
  # MPI gives no explicit interface to the ones taking a choice buffer, nor to
  # MPI_Cart_coords, so those names cannot be re-exported from m_mpi by use
  # association.  Ask this MPI whether it has them and let src/mpi.f90 declare
  # them external where it does not, rather than guessing from the compiler.
  include(CheckFortranSourceCompiles)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_Fortran)
  check_fortran_source_compiles("
program main
  use mpi, only: MPI_Allreduce, MPI_Alltoall, MPI_Bcast, MPI_Cart_coords, &
                 MPI_Irecv, MPI_Isend, MPI_Reduce, MPI_Sendrecv, MPI_Waitall
  implicit none
end program main"
    MPI_MODULE_HAS_INTERFACES SRC_EXT f90)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT MPI_MODULE_HAS_INTERFACES)
    add_definitions("-DMPI_MODULE_LACKS_INTERFACES")
  endif()
endif()

# Toggleable build type (default Release)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose type of build" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")

# Single precision option (default OFF)
option(SINGLE_PREC "Enable single precision" OFF)
option(X3D2_FAST_MATH "Use fast-math Release flags (GNU, Flang, Intel and NVHPC)" OFF)
if(SINGLE_PREC)
  add_definitions("-DSINGLE_PREC")
endif()

# Backend for GPU offloading options (default OFF)
set(ENABLE_BACKEND OFF CACHE STRING "")
set_property(CACHE ENABLE_BACKEND PROPERTY STRINGS OFF CUDA OMP_TGT)

# ADIOS2 option (default OFF)
option(WITH_ADIOS2 "Enable ADIOS2 for I/O" OFF)
option(WITH_ADIOS2_GPU_AWARE "Enable GPU-aware ADIOS2 I/O with the CUDA backend" OFF)

if(WITH_ADIOS2_GPU_AWARE)
  if(NOT WITH_ADIOS2)
    message(FATAL_ERROR "WITH_ADIOS2_GPU_AWARE=ON requires WITH_ADIOS2=ON")
  endif()
  if(NOT ENABLE_BACKEND STREQUAL "CUDA")
    message(FATAL_ERROR "WITH_ADIOS2_GPU_AWARE=ON requires ENABLE_BACKEND=CUDA")
  endif()
  add_compile_definitions(X3D2_ADIOS2_CUDA)
endif()

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
  set(BACKEND_ARCH "" CACHE STRING "Target GPU architecture, e.g. cc80 for NVIDIA, gfx942 for AMD, or pvc / spir64 for Intel")
  if ("${BACKEND_ARCH}" STREQUAL "")
    message(FATAL_ERROR
       "CUDA and OMP_TGT backends require BACKEND_ARCH=<arch> (replaces CUDA_ARCH / OMP_TGT_ARCH), e.g. cc80, gfx942 or pvc")
  endif()
endif()
# NVHPC picks the CUDA toolkit matching the installed driver. On a node without
# a GPU (driver version 0) it falls back to its configured default, which may
# not be bundled with the SDK. Naming a bundled version (e.g. 12.6) is a
# fallback for that case and lets the code compile without a GPU.
set(BACKEND_CUDA_VERSION "" CACHE STRING "Fallback CUDA toolkit version for NVHPC -gpu=cudaX.Y (empty: compiler default)")
set(NVHPC_GPU_FLAG "-gpu=${BACKEND_ARCH}")
if (NOT "${BACKEND_CUDA_VERSION}" STREQUAL "")
  set(NVHPC_GPU_FLAG "${NVHPC_GPU_FLAG},cuda${BACKEND_CUDA_VERSION}")
endif()
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    message(FATAL_ERROR "CUDA backend require NVHPC and PGI compiler")
  endif()
  if(NOT BACKEND_ARCH MATCHES "^cc[0-9]+$")
    message(FATAL_ERROR "For Nvidia GPU, BACKEND_ARCH must be specified as ccXX")
  endif()
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  # Classic ifort has no device offload; ifx (IntelLLVM) does, see below.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(FATAL_ERROR
      "The OMP_TGT backend needs ifx (IntelLLVM): classic ifort cannot offload. "
      "Build the CPU backend with ENABLE_BACKEND=OFF.")
  endif()
  # FindOpenMP's version probe fails to compile with ifx (OpenMP_Fortran_VERSION
  # comes back empty), and every ifx that offloads implements OpenMP 5.x.
  if("${OpenMP_Fortran_VERSION}" STREQUAL "" AND
     CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(STATUS "OpenMP Fortran version not detected with ifx: assuming 5.x")
  elseif("${OpenMP_Fortran_VERSION}" VERSION_LESS "4.5")
    message(FATAL_ERROR "OpenMP Fortran version ${OpenMP_Fortran_VERSION} is detected, but at least 4.5 is required with OMP_TGT.")
  endif()
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
     ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    if(NOT BACKEND_ARCH MATCHES "^gfx[0-9]+[a-f]?$")
      message(FATAL_ERROR "For AMD GPU, BACKEND_ARCH must be specified as gfxXXX")
    endif()
    if(NOT "${BACKEND_CUDA_VERSION}" STREQUAL "")
      message(WARNING "BACKEND_CUDA_VERSION applies to NVHPC only and is ignored with ${CMAKE_Fortran_COMPILER_ID}")
    endif()
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    # BACKEND_ARCH=spir64 builds SPIR-V that the GPU driver compiles at first
    # launch (JIT) and runs on any Intel GPU.  Any other value names the device
    # to compile for ahead of time with ocloc, e.g. pvc (Data Center GPU Max)
    # or bmg (Arc B-series); `ocloc compile --help` lists the names.
    if(NOT BACKEND_ARCH MATCHES "^[A-Za-z0-9_.-]+$")
      message(FATAL_ERROR
        "For Intel GPU, BACKEND_ARCH must be spir64 (JIT) or an ocloc device name, e.g. pvc")
    endif()
    if(NOT "${BACKEND_CUDA_VERSION}" STREQUAL "")
      message(WARNING "BACKEND_CUDA_VERSION applies to NVHPC only and is ignored with ${CMAKE_Fortran_COMPILER_ID}")
    endif()
    if(BACKEND_ARCH STREQUAL "spir64")
      set(_intel_offload_flags "-fopenmp-targets=spir64")
    else()
      set(_intel_offload_flags
        "-fopenmp-targets=spir64_gen"
        "SHELL:-Xopenmp-target-backend \"-device ${BACKEND_ARCH}\"")
    endif()
    # The offload targets must reach every compile and the final link, so hang
    # them on the imported target as the NVHPC branch below does.
    set_property(TARGET OpenMP::OpenMP_Fortran APPEND PROPERTY
      INTERFACE_COMPILE_OPTIONS ${_intel_offload_flags})
    set_property(TARGET OpenMP::OpenMP_Fortran APPEND PROPERTY
      INTERFACE_LINK_OPTIONS ${_intel_offload_flags})
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    if(NOT BACKEND_ARCH MATCHES "^cc[0-9]+$")
      message(FATAL_ERROR "For Nvidia GPU, BACKEND_ARCH must be specified as ccXX")
    endif()
    # FindOpenMP supplies -mp, which targets the CPU with NVHPC.  Compile and
    # link target regions for the GPU when this backend is selected.
    set_target_properties(OpenMP::OpenMP_Fortran PROPERTIES
      INTERFACE_COMPILE_OPTIONS "-mp=gpu;${NVHPC_GPU_FLAG}"
      INTERFACE_LINK_OPTIONS "-mp=gpu;${NVHPC_GPU_FLAG}")
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
mark_as_advanced(FLANG_OPENMP_RUNTIME_DIR)
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
if(WITH_ADIOS2_GPU_AWARE AND ADIOS2_FOUND AND NOT TARGET adios2::core_cuda)
  message(FATAL_ERROR "WITH_ADIOS2_GPU_AWARE=ON requires CUDA-enabled ADIOS2")
endif()
add_subdirectory(src)
add_subdirectory(tests)
