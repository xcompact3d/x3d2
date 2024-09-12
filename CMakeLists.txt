cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran CXX C)
enable_testing()

# TODO switch adios between system install and this installed version

include(FetchContent)

FetchContent_Declare(
  Adios2
  GIT_REPOSITORY https://github.com/ornladios/ADIOS2.git
  GIT_TAG        v2.10.0-rc1
)
FetchContent_MakeAvailable(Adios2)

add_subdirectory(src)
add_subdirectory(tests)
