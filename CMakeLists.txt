cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran)
enable_testing()

#
# Set the Poisson solver choice
#
set(POISSON_SOLVER FFT CACHE STRING
  "Select the Poisson solver: FFT or ITER")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
   set(BACKEND CUDA)
else()
   set(BACKEND OMP)
endif()

add_subdirectory(src)
add_subdirectory(tests)
