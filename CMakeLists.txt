cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran)
enable_testing()

#
# Set the Poisson solver choice
#
set(POISSON_SOLVER FFT CACHE STRING
  "Select the Poisson solver: FFT or ITER")

set(BACKEND None CACHE STRING
  "Select the Backend: CUDA or OMP")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
   set(COMPILER CUDA)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   set(COMPILER GNU)
else()
   set(COMPILER OTHER)
endif()

add_subdirectory(src)
add_subdirectory(tests)
