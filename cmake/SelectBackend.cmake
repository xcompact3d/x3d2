# Defines the backend selection for x3d2.
#
# Currently available backends
# - GPU/CUDAFortran (nvfort only)
# - CPU/OMP (all)

set(X3D2_BACKEND_OPTIONS "OMP/CPU")
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
    ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  set(X3D2_BACKEND_DEFAULT "CUDA/GPU")
  list(APPEND X3D2_BACKEND_OPTIONS ${X3D2_BACKEND_DEFAULT})
else()
  set(X3D2_BACKEND_DEFAULT "OMP/CPU")
endif()

set(X3D2_BACKEND
  ${X3D2_BACKEND_DEFAULT}
  CACHE
  STRING
  "Select the backend. Current options are ${X3D2_BACKEND_OPTIONS}")
set_property(CACHE
  X3D2_BACKEND
  PROPERTY STRINGS
  ${X3D2_BACKEND_OPTIONS})
message(STATUS "Using ${X3D2_BACKEND} backend")
