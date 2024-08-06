set(SRC
  allocator.f90
  backend.f90
  common.f90
  poisson_fft.f90
  solver.f90
  tdsops.f90
  time_integrator.f90
  ordering.f90
  mesh.f90
  field.f90
  omp/backend.f90
  omp/common.f90
  omp/kernels/distributed.f90
  omp/poisson_fft.f90
  omp/sendrecv.f90
  omp/exec_dist.f90
  omp/exec_thom.f90
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

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  list(APPEND SRC ${CUDASRC})
endif()

add_library(x3d2 STATIC ${SRC})
target_include_directories(x3d2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

add_executable(xcompact xcompact.f90)
target_link_libraries(xcompact PRIVATE x3d2)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
 
  set(CMAKE_Fortran_FLAGS "-cpp -cuda")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds -Mchkptr -Ktrap=fp")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  target_link_options(x3d2 INTERFACE "-cuda")
  target_link_options(x3d2 INTERFACE "-cudalib=cufftmp")

  target_compile_options(xcompact PRIVATE "-DCUDA")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Wno-unused-dummy-argument")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math")
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(x3d2 PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(xcompact PRIVATE OpenMP::OpenMP_Fortran)

find_package(MPI REQUIRED)
target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
target_link_libraries(xcompact PRIVATE MPI::MPI_Fortran)

