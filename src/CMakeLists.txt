set(SRC
  allocator.f90
  thomas.f90
)
set(CUDASRC
  cuda/cuda_allocator.f90
)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  list(APPEND SRC ${CUDASRC})
  target_compile_options(x3d2 "-cuda")
endif()

add_library(x3d2 STATIC ${SRC})
