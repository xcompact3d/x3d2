cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran)
enable_testing()

set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

add_subdirectory(src)
add_subdirectory(tests)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SelectBackend)
