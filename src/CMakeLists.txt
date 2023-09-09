set(SRC
  allocator.f90
  common.f90
  derparams.f90
  thomas.f90
  stencil.f90
  stencil_definitions.f90
  diffengine.f90
  vector3d.f90
  vector3d_simd.f90
)
set(CUDASRC
  cuda/cuda_allocator.f90
  cuda/kernels_dist.f90
)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  list(APPEND SRC ${CUDASRC})
endif()

add_library(x3d2 STATIC ${SRC})
target_include_directories(x3d2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  target_compile_options(x3d2 PRIVATE "-cuda")
  target_compile_options(x3d2 PRIVATE "-O3")
  target_compile_options(x3d2 PRIVATE "-fast")
  target_link_options(x3d2 INTERFACE "-cuda")
endif()

find_package(MPI REQUIRED)
target_link_libraries(x3d2 PRIVATE MPI::MPI_Fortran)
