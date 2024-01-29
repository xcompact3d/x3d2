set(SRC
  allocator.f90
  backend.f90
  common.f90
  solver.f90
  tdsops.f90
  time_integrator.f90
  omp/backend.f90
  omp/common.f90
  omp/kernels_dist.f90
  omp/sendrecv.f90
  omp/exec_dist.f90
)
set(CUDASRC
  cuda/backend.f90
  cuda/common.f90
  cuda/allocator.f90
  cuda/exec_dist.f90
  cuda/kernels_dist.f90
  cuda/kernels_trans.f90
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

target_compile_options(x3d2 PRIVATE "-O3")
target_compile_options(xcompact PRIVATE "-O3")
target_compile_options(xcompact PRIVATE "-cpp")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  target_compile_options(x3d2 PRIVATE "-cuda")
  target_compile_options(x3d2 PRIVATE "-fast")
  target_link_options(x3d2 INTERFACE "-cuda")
  target_compile_options(xcompact PRIVATE "-cuda")
  target_compile_options(xcompact PRIVATE "-fast")

  target_compile_options(xcompact PRIVATE "-DCUDA")
  #  target_link_options(xcompact INTERFACE "-cuda")
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  target_compile_options(x3d2 PRIVATE "-ffast-math")
  target_compile_options(xcompact PRIVATE "-ffast-math")
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(x3d2 PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(xcompact PRIVATE OpenMP::OpenMP_Fortran)

find_package(MPI REQUIRED)
target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
target_link_libraries(xcompact PRIVATE MPI::MPI_Fortran)

