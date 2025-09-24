set(SRC
  allocator.f90
  common.f90
  config.f90
  field.f90
  mesh.f90
  ordering.f90
  poisson_fft.f90
  solver.f90
  tdsops.f90
  time_integrator.f90
  ordering.f90
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
  case/generic.f90
  case/tgv.f90
  module/ibm.f90
)
set(CUDASRC
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
set(BACKENDSRC
  backend/omp/backend.f90
  )
set(2DECOMPFFTSRC
  backend/omp/poisson_fft.f90
  decomp/decomp_2decompfft.f90
)
set(DUMMYDECOMPSRC
  decomp/decomp_dummy.f90
)

set(IOSRC
  io/io_base.f90
  io/io_dummy.f90
  io/io_service.f90
  io/checkpoint_io.f90
)

if(WITH_ADIOS2)
  list(APPEND IOSRC io/io_adios2.f90)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  list(APPEND SRC ${CUDASRC})
  list(APPEND BACKENDSRC backend/cuda/backend.f90)
endif()

# Decide whether 2decomp&fft is supported by the build
find_package(decomp2d)

if(WITH_2DECOMPFFT)
  list(APPEND SRC ${2DECOMPFFTSRC})
else()
  list(APPEND SRC ${DUMMYDECOMPSRC})
endif()

list(APPEND SRC ${IOSRC})

add_library(x3d2 STATIC ${SRC})
target_include_directories(x3d2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

add_library(x3d2_backends STATIC ${BACKENDSRC})
target_include_directories(x3d2_backends INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(x3d2_backends PRIVATE x3d2)

add_executable(xcompact xcompact.f90)
target_link_libraries(xcompact PRIVATE x3d2)
target_link_libraries(xcompact PRIVATE x3d2_backends)

# if CUDA compiler
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  set(CMAKE_Fortran_FLAGS "-cpp -cuda")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds -Mchkptr -Ktrap=fp")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  target_link_options(x3d2 INTERFACE "-cuda")
  target_link_options(x3d2 INTERFACE "-cudalib=cufftmp")

  target_compile_options(xcompact PRIVATE "-DCUDA")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
       ${CMAKE_Fortran_COMPILER_ID} STREQUAL "LLVMFlang")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Wno-unused-dummy-argument")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS "-eF -M878") # -M878 suppresses WARNING multiple module includes (not useful)
  set(CMAKE_Fortran_FLAGS_DEBUG "-G2 -O0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  target_link_options(x3d2 INTERFACE -h omp)
endif()

if(WITH_2DECOMPFFT)
  message(STATUS "Using the FFT poisson solver with 2decomp&fft support on the OpenMP backend")
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
  find_package(decomp2d REQUIRED)
  include_directories(${decomp2d_INCLUDE_DIRS})
  target_link_libraries(decomp2d)
  target_link_libraries(x3d2 PRIVATE decomp2d)
  target_compile_options(x3d2 PRIVATE "-DWITH_2DECOMPFFT")
  target_compile_options(x3d2_backends PRIVATE "-DWITH_2DECOMPFFT")
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(x3d2 PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(x3d2_backends PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(xcompact PRIVATE OpenMP::OpenMP_Fortran)

find_package(MPI REQUIRED)
target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
target_link_libraries(xcompact PRIVATE MPI::MPI_Fortran)

if(WITH_ADIOS2)
  target_link_libraries(x3d2 PRIVATE adios2::fortran_mpi)
  target_link_libraries(xcompact PRIVATE adios2::fortran_mpi)
  target_compile_options(x3d2 PRIVATE "-DWITH_ADIOS2")
  target_compile_options(x3d2_backends PRIVATE "-DWITH_ADIOS2")
endif()
