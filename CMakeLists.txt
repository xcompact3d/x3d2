cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran)
set(CMAKE_Fortran_STANDARD 2003)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
enable_testing()

add_subdirectory(src)
if (CUDA)
  add_subdirectory(src/cuda)
endif()

add_subdirectory(tests)
