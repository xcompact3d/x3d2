# 1. Set source files
set(SRC
  allocator.f90
  common.f90
  config.f90
  mpi.f90
  field.f90
  mesh.f90
  ordering.f90
  poisson_fft.f90
  solver.f90
  les/les.f90
  postprocess/postprocess.f90
  postprocess/monitoring.f90
  postprocess/scalar_series.f90
  tdsops.f90
  time_integrator.f90
  mesh.f90
  mesh_content.f90
  field.f90
  vector_calculus.f90
  backend/backend.f90
  backend/omp/common.f90
  backend/omp/kernels/distributed.f90
  backend/omp/kernels/spectral_processing.f90
  backend/omp/kernels/thomas.f90
  backend/omp/sendrecv.f90
  backend/omp/exec_dist.f90
  backend/omp/exec_thom.f90
  case/base_case.f90
  case/channel.f90
  case/cylinder.f90
  case/generic.f90
  case/tgv.f90
  module/ibm.f90
)
set(BACKENDSRC
  backend/omp/backend.f90
)
# CUDA specific
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  list(APPEND SRC
    backend/cuda/allocator.f90
    backend/cuda/common.f90
    backend/cuda/exec_dist.f90
    backend/cuda/exec_thom.f90
    backend/cuda/kernels/distributed.f90
    backend/cuda/kernels/fieldops.f90
    backend/cuda/kernels/reorder.f90
    backend/cuda/kernels/spectral_processing.f90
    backend/cuda/kernels/thomas.f90
    backend/cuda/poisson_fft.f90
    backend/cuda/sendrecv.f90
    backend/cuda/tdsops.f90
  )
  list(APPEND BACKENDSRC
    backend/cuda/backend.f90)
# OMP TGT specific
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  list(APPEND SRC
    backend/omp/target/allocator.f90
  )
  list(APPEND BACKENDSRC
    backend/omp/target/backend.f90)
endif()

# 2Decomp
if(WITH_2DECOMPFFT)
  list(APPEND SRC
    backend/omp/poisson_fft.f90
    decomp/decomp_2decompfft.f90
  )
else()
  list(APPEND SRC
    decomp/decomp_dummy.f90
  )
endif()

list(APPEND SRC
  io/nvtx.f90
  io/io_base.f90
  io/io_session.f90
  io/io_field_utils.f90
  io/checkpoint_manager.f90
  io/checkpoint_state.f90
  io/snapshot_manager.f90
  io/stats.f90
  io/io_manager.f90
)

# ADIOS2
if(WITH_ADIOS2)
  list(APPEND SRC
    io/adios2/io.f90
  )
  if(WITH_ADIOS2_GPU_AWARE)
    list(APPEND SRC
      io/adios2/adios2_device_put.c
    )
  endif()
else()
  list(APPEND SRC
  io/dummy/io.f90
  )
endif()

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")

# 2. Create the main library
add_library(x3d2 STATIC ${SRC})
set_target_properties(x3d2 PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2")
target_include_directories(x3d2 INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2)

add_library(x3d2_backends STATIC ${BACKENDSRC})
set_target_properties(x3d2_backends PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_backends")
target_include_directories(x3d2_backends INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}/modules/x3d2_backends)
target_link_libraries(x3d2_backends PRIVATE x3d2)

add_executable(xcompact xcompact.f90)
target_link_libraries(xcompact PRIVATE x3d2_backends)
target_link_libraries(xcompact PRIVATE x3d2)

# 3. Set compiler specific flags and options
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
         ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
         ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018")
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
      ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wno-ambiguous-structure-constructor -Werror")
  else()
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Wno-unused-dummy-argument")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -floop-optimize -march=native")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" OR
       ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "-cpp")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS "-eF -M878") # -M878 suppresses WARNING multiple module includes (not useful)
  set(CMAKE_Fortran_FLAGS_DEBUG "-G2 -O0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  target_link_options(x3d2 INTERFACE -h omp)
  target_link_options(x3d2_backends INTERFACE -h omp)
else()
  message(FATAL_ERROR "Unsupported compiler")
endif()

# 4. Set backend specific flags and options
if (${ENABLE_BACKEND} STREQUAL "CUDA")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds")
  # Opt-in: -Ktrap=fp traps FP exceptions process-wide, turning benign FP ops in
  # linked libs (cuFFTMp PTX-JIT, HPC-X/UCX at MPI_Init) into fatal SIGFPEs.
  option(CUDA_DEBUG_FP_TRAP "Trap FP exceptions in the CUDA debug build (-Ktrap=fp)" OFF)
  if(CUDA_DEBUG_FP_TRAP)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Ktrap=fp")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  target_link_options(x3d2 INTERFACE "-cuda")
  # cuFFTMp distributes a transform over MPI ranks and its nvfortran wrapper
  # calls MPI_Comm_f2c, so it cannot link into a build without MPI.  A serial
  # build is one rank and links plain cuFFT instead: see create_fft_plan in
  # backend/cuda/poisson_fft.f90 for the matching runtime path.
  if(WITH_MPI)
    target_link_options(x3d2 INTERFACE "-cudalib=cufftmp")
  else()
    target_link_options(x3d2 INTERFACE "-cudalib=cufft")
  endif()
  target_compile_options(xcompact PRIVATE "-DCUDA")
  # NVTX ranges label the I/O paths in Nsight Systems (see io/nvtx.f90).
  target_compile_definitions(x3d2 PRIVATE X3D2_NVTX)
  find_package(CUDAToolkit QUIET)
  if(CUDAToolkit_FOUND AND TARGET CUDA::nvToolsExt)
    target_link_libraries(x3d2 PUBLIC CUDA::nvToolsExt)
  else()
    target_link_libraries(x3d2 PUBLIC nvToolsExt)
  endif()
elseif (${ENABLE_BACKEND} STREQUAL "OMP_TGT")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
    add_definitions("-DOMP_TGT_NVIDIA")
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
         ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang" OR
         ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp --offload-arch=${OMP_TGT_ARCH}")
    add_definitions("-DOMP_TGT_AMD") # Assume AMD GPU
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
    add_definitions("-DOMP_TGT_AMD") # Assume AMD GPU
  endif()
endif()

# OpenMP is resolved once at the top level (see the root CMakeLists.txt).
target_link_libraries(x3d2 PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(x3d2_backends PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(xcompact PRIVATE OpenMP::OpenMP_Fortran)

if(WITH_MPI)
  target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
  target_link_libraries(x3d2_backends PRIVATE MPI::MPI_Fortran)
  target_link_libraries(xcompact PRIVATE MPI::MPI_Fortran)
endif()

# 5. Set thirdparty libary dependencies
if(WITH_2DECOMPFFT)
  message(STATUS "Using the FFT poisson solver with 2decomp&fft support on the OpenMP backend")
  target_link_libraries(x3d2 PRIVATE decomp2d)
  target_link_libraries(x3d2_backends PRIVATE decomp2d)
  target_compile_options(x3d2 PRIVATE "-DWITH_2DECOMPFFT")
  target_compile_options(x3d2_backends PRIVATE "-DWITH_2DECOMPFFT")
endif()

if(WITH_ADIOS2)
  target_link_libraries(x3d2 PRIVATE ${ADIOS2_FORTRAN_TARGET})
  target_link_libraries(xcompact PRIVATE ${ADIOS2_FORTRAN_TARGET})
  target_compile_options(x3d2 PRIVATE "-DWITH_ADIOS2")

  # GPU-aware I/O requires C API (called via iso_c_binding)
  if(WITH_ADIOS2_GPU_AWARE)
    target_link_libraries(x3d2 PRIVATE ${ADIOS2_C_TARGET})
    target_link_libraries(xcompact PRIVATE ${ADIOS2_C_TARGET})
  endif()
endif()
