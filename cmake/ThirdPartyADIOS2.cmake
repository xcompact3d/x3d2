# Provide the ADIOS2 Fortran/MPI bindings as the target `adios2::fortran_mpi`.
#
# Set USE_SYSTEM_ADIOS2=ON, optionally with -DADIOS2_ROOT_DIR=<prefix>, to link
# a pre-built install, otherwise ADIOS2 is downloaded and built into the build
# tree as part of the build rather than of the configuration.  Either way,
# linking `adios2::fortran_mpi` is all a consumer has to do.

if(WITH_ADIOS2)
  set(adios2_version "v2.12.1")

  option(USE_SYSTEM_ADIOS2 "Use system-installed ADIOS2" OFF)
  set(ADIOS2_ROOT_DIR "" CACHE PATH
    "Directory where ADIOS2 is installed, only used when USE_SYSTEM_ADIOS2=ON")

  mark_as_advanced(ADIOS2_ROOT_DIR)

  # find_package() caches ADIOS2_DIR, so a previous backend or a stale system
  # install can otherwise satisfy this lookup even when the search path changed.
  unset(ADIOS2_DIR CACHE)
  unset(ADIOS2_DIR)
  mark_as_advanced(ADIOS2_DIR)

  if(USE_SYSTEM_ADIOS2)
    if(ADIOS2_ROOT_DIR)
      message(STATUS "Looking for ADIOS2 in ${ADIOS2_ROOT_DIR}")
      find_package(ADIOS2 CONFIG
                   COMPONENTS Fortran MPI
                   PATHS "${ADIOS2_ROOT_DIR}"
                   NO_DEFAULT_PATH
                   QUIET)
    else()
      find_package(ADIOS2 CONFIG COMPONENTS Fortran MPI QUIET)
    endif()

    if(ADIOS2_FOUND)
      message(STATUS "ADIOS2 FOUND in ${ADIOS2_DIR}")
    else()
      message(FATAL_ERROR "USE_SYSTEM_ADIOS2 is ON but ADIOS2 was not found. "
        "Please install ADIOS2, specify ADIOS2_ROOT_DIR, "
        "or set USE_SYSTEM_ADIOS2=OFF to build it from source.")
    endif()
  else()
    if(ADIOS2_ROOT_DIR)
      message(WARNING "ADIOS2_ROOT_DIR is set but will be ignored because USE_SYSTEM_ADIOS2=OFF")
    endif()

    # A CUDA-enabled ADIOS2 cannot be reused by a CPU build, and vice versa, so
    # the two get separate install trees.
    if(${ENABLE_BACKEND} STREQUAL "CUDA")
      set(adios2_config_suffix "cuda")
    else()
      set(adios2_config_suffix "cpu")
    endif()

    set(adios2_install_dir
      "${CMAKE_CURRENT_BINARY_DIR}/adios2-${adios2_config_suffix}-${adios2_version}")

    find_package(ADIOS2 CONFIG
                 COMPONENTS Fortran MPI
                 PATHS "${adios2_install_dir}"
                 NO_DEFAULT_PATH
                 QUIET)
    if(ADIOS2_FOUND)
      message(STATUS "ADIOS2 FOUND in ${adios2_install_dir}")
    else(ADIOS2_FOUND)
      message(STATUS "Building ADIOS2 from source")

      if(${ENABLE_BACKEND} STREQUAL "CUDA")
        # ADIOS2 defaults CMAKE_CUDA_ARCHITECTURES to 52 (Maxwell) when unset,
        # which CUDA >= 13 no longer supports, so an architecture is always
        # forwarded.  Native detection compiles for the GPU of the build machine.
        set(CUDA_ARCH "native" CACHE STRING
          "CUDA architecture(s) for the ADIOS2 build, e.g. 80")

        set(adios2_cuda_args
          "-DADIOS2_USE_CUDA=ON"
          "-DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCH}")
      else()
        set(adios2_cuda_args "-DADIOS2_USE_CUDA=OFF")
      endif()

      # The install tree is only populated during the build, so the libraries
      # have to be named up front.  Pin the layout the external project installs
      # into: GNUInstallDirs picks lib64 over lib on some distributions.
      set(adios2_fortran_library
        "${adios2_install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}adios2_fortran${CMAKE_SHARED_LIBRARY_SUFFIX}")
      set(adios2_fortran_mpi_library
        "${adios2_install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}adios2_fortran_mpi${CMAKE_SHARED_LIBRARY_SUFFIX}")
      if(WITH_ADIOS2_GPU_AWARE)
        set(adios2_c_mpi_library
          "${adios2_install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}adios2_c_mpi${CMAKE_SHARED_LIBRARY_SUFFIX}")
      endif()

      include(ExternalProject)

      ExternalProject_Add(adios2-${adios2_version}
        GIT_REPOSITORY    "https://github.com/ornladios/ADIOS2.git"
        GIT_TAG           "${adios2_version}"
        GIT_SHALLOW       TRUE
        SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/adios2-src"
        BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/adios2-build"
        INSTALL_DIR       "${adios2_install_dir}"
        CMAKE_ARGS        "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
                          "-DCMAKE_INSTALL_LIBDIR=lib"
                          "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                          "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
                          "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                          "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
                          "-DBUILD_SHARED_LIBS=ON"
                          "-DBUILD_TESTING=OFF"
                          "-DADIOS2_BUILD_EXAMPLES=OFF"
                          "-DADIOS2_USE_Fortran=ON"
                          "-DADIOS2_USE_MPI=ON"
                          "-DADIOS2_USE_HDF5=OFF"
                          "-DADIOS2_USE_CURL=OFF"
                          ${adios2_cuda_args}
        BUILD_BYPRODUCTS  "${adios2_fortran_library}" "${adios2_fortran_mpi_library}"
                          ${adios2_c_mpi_library}
        TEST_COMMAND      ""
        # The git update step carries no stamp, so it is always out of date and
        # drags the configure, build and install steps with it on every `make`.
        # Disconnecting it pins the checkout to the tag already cloned; re-fetch
        # deliberately with the `adios2-<version>-update` target.
        UPDATE_DISCONNECTED TRUE
      )

      # The module files land here during the build, but the directory has to
      # exist before it can be used as an include directory.
      file(MAKE_DIRECTORY "${adios2_install_dir}/include/adios2/fortran")

      # Stands in for the package config that is not installed yet, exporting
      # what it would export.  Its dependency on the external project is
      # followed by everything linking it, so ADIOS2 is built first.  The alias
      # gives consumers the same name the package config provides, and is global
      # already, so sibling directories (e.g. tests/) can link it.
      add_library(adios2_fortran_mpi INTERFACE)
      target_include_directories(adios2_fortran_mpi INTERFACE
        "${adios2_install_dir}/include/adios2/fortran")
      target_link_libraries(adios2_fortran_mpi INTERFACE
        "${adios2_fortran_mpi_library}" "${adios2_fortran_library}" MPI::MPI_Fortran)
      target_compile_definitions(adios2_fortran_mpi INTERFACE ADIOS2_USE_MPI)
      add_dependencies(adios2_fortran_mpi adios2-${adios2_version})
      add_library(adios2::fortran_mpi ALIAS adios2_fortran_mpi)
      if(WITH_ADIOS2_GPU_AWARE)
        add_library(adios2_c_mpi INTERFACE)
        target_link_libraries(adios2_c_mpi INTERFACE
          "${adios2_c_mpi_library}" MPI::MPI_C)
        add_dependencies(adios2_c_mpi adios2-${adios2_version})
        add_library(adios2::c_mpi ALIAS adios2_c_mpi)
      endif()
    endif(ADIOS2_FOUND)
  endif()

  # find_package() imports the targets into this directory only, so promote them
  # to let sibling directories (e.g. tests/) link them.  The interface library
  # built above needs no promotion: its alias is global already.
  if(ADIOS2_FOUND)
    foreach(adios2_target adios2::c adios2::fortran adios2::fortran_mpi)
      if(TARGET ${adios2_target})
        get_target_property(adios2_target_is_global ${adios2_target} IMPORTED_GLOBAL)
        if(NOT adios2_target_is_global)
          set_target_properties(${adios2_target} PROPERTIES IMPORTED_GLOBAL TRUE)
        endif()
      endif()
    endforeach()
  endif()
else()
  message(STATUS "ADIOS2 is disabled")
endif()
