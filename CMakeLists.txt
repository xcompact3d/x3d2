cmake_minimum_required(VERSION 3.10)
project(x3d2 LANGUAGES Fortran CXX C)
enable_testing()

# option to toggle between system-installed ADIOS2 and fetching it
option(USE_SYSTEM_ADIO2 "Use system-installed ADIOS2" OFF)

set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

# detect MPI for parallel support
#find_package(MPI COMPONENTS Fortran REQUIRED)

# handling ADIOS2: system-installed vs. FetchContent
#if(USE_SYSTEM_ADIOS2)
#    find_package(ADIOS2 REQUIRED COMPONENTS Fortran MPI)
#else()
#    include(FetchContent)
#    FetchContent_Declare(
#        ADIOS2
#        GIT_REPOSITORY https://github.com/ornladios/ADIOS2.git
#        GIT_TAG        v2.10.2
#    )
#    FetchContent_MakeAvailable(ADIOS2)
#endif()

add_subdirectory(src)
add_subdirectory(tests)
