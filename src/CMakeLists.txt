set(SRC
  allocator.f90
  backend.f90
  common.f90
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
  decomp/decomp.f90
  field.f90
  vector_calculus.f90
  omp/backend.f90
  omp/common.f90
  omp/kernels/distributed.f90
  omp/kernels/spectral_processing.f90
  omp/sendrecv.f90
  omp/exec_dist.f90
)
set(CUDASRC
  cuda/allocator.f90
  cuda/backend.f90
  cuda/common.f90
  cuda/exec_dist.f90
  cuda/exec_thom.f90
  cuda/kernels/distributed.f90
  cuda/kernels/reorder.f90
  cuda/kernels/spectral_processing.f90
  cuda/kernels/thomas.f90
  cuda/poisson_fft.f90
  cuda/sendrecv.f90
  cuda/tdsops.f90
)
set(2DECOMPFFTSRC
  omp/poisson_fft.f90
  decomp/decomp_2decompfft.f90
)
set(GENERICDECOMPSRC
  decomp/decomp_dummy.f90
)

if(${BACKEND} STREQUAL "CUDA")
  list(APPEND SRC ${CUDASRC})
endif()

# Decide whether 2decomp&fft is supported by the build
find_package(decomp2d)
set(2DECOMPFFT_SUPPORT (${POISSON_SOLVER} STREQUAL "FFT" AND decomp2d_FOUND))

if(2DECOMPFFT_SUPPORT)
  list(APPEND SRC ${2DECOMPFFTSRC})
else()
  list(APPEND SRC ${GENERICDECOMPSRC})
endif()

add_library(x3d2 STATIC ${SRC})
target_include_directories(x3d2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

add_executable(xcompact xcompact.f90)
target_link_libraries(xcompact PRIVATE x3d2)

if(${COMPILER} STREQUAL "CUDA")
  set(CMAKE_Fortran_FLAGS "-cpp -cuda")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds -Mchkptr -Ktrap=fp")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  target_link_options(x3d2 INTERFACE "-cuda")
  target_link_options(x3d2 INTERFACE "-cudalib=cufftmp")

  target_compile_options(xcompact PRIVATE "-DCUDA")
elseif(${COMPILER} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Wno-unused-dummy-argument")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math")
endif()

if(2DECOMPFFT_SUPPORT)
  message(STATUS "Using the FFT poisson solver with 2decomp&fft support")
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
  find_package(decomp2d REQUIRED)
  include_directories(${decomp2d_INCLUDE_DIRS})
  target_link_libraries(decomp2d)
  target_link_libraries(x3d2 PRIVATE decomp2d)
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(x3d2 PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(xcompact PRIVATE OpenMP::OpenMP_Fortran)

find_package(MPI REQUIRED)
target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
target_link_libraries(xcompact PRIVATE MPI::MPI_Fortran)
