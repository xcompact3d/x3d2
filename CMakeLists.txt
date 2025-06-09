cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran)
enable_testing()

#
# Set the Poisson solver choice
#
set(POISSON_SOLVER FFT CACHE STRING
  "Select the Poisson solver: FFT or ITER")
set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

add_subdirectory(src)
add_subdirectory(tests)
