# Provide the 2decomp-fft library as the target `decomp2d`.
#
# Point -Ddecomp2d_install_dir=<prefix> at a pre-built install to use it,
# otherwise the library is downloaded and built into the build tree.  Either
# way, linking `decomp2d` is all a consumer has to do.

# Decomp2D version
if(WITH_2DECOMPFFT)
  set(decomp2d_version "v2.0.3.2")

  # Set install directory, if pre-build directory not specified by user.
  if(NOT decomp2d_install_dir)
    if(SINGLE_PREC)
        set(decomp2d_install_dir "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt-sp-${decomp2d_version}")
    else()
        set(decomp2d_install_dir "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt-dp-${decomp2d_version}")
    endif()
  else()
    message("Found decomp2d install directory: ${decomp2d_install_dir}")
  endif()

  # find_package() caches decomp2d_DIR, so a previous precision can otherwise
  # satisfy this lookup even when decomp2d_install_dir has changed.
  unset(decomp2d_DIR CACHE)
  unset(decomp2d_DIR)
  mark_as_advanced(decomp2d_DIR)

  find_package(decomp2d CONFIG
               PATHS ${decomp2d_install_dir}
               NO_DEFAULT_PATH
               QUIET)
  if (decomp2d_FOUND)
    message(STATUS "2decomp-fft FOUND in ${decomp2d_install_dir}")
  else(decomp2d_FOUND)
    message(STATUS "Building 2decomp-fft from source")

    if(SINGLE_PREC)
      set(DOUBLE_PRECISION OFF)
    else()
      set(DOUBLE_PRECISION ON)
    endif()

    # The install tree is only populated during the build, so the library has to
    # be named up front.  Pin the layout the external project installs into:
    # GNUInstallDirs picks lib64 over lib on some distributions.
    set(decomp2d_library
      "${decomp2d_install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}decomp2d${CMAKE_STATIC_LIBRARY_SUFFIX}")

    include(ExternalProject)

    ExternalProject_Add(2decomp-${decomp2d_version}
      GIT_REPOSITORY    "https://github.com/xcompact3d/2decomp-fft"
      GIT_TAG           "${decomp2d_version}"
      GIT_SHALLOW       TRUE
      SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-src"
      BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build"
      INSTALL_DIR       "${decomp2d_install_dir}"
      CMAKE_ARGS        "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
                        "-DCMAKE_INSTALL_LIBDIR=lib"
                        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                        "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
                        "-DBUILD_SHARED_LIBS=OFF"
                        "-DDOUBLE_PRECISION=${DOUBLE_PRECISION}"
      BUILD_BYPRODUCTS  "${decomp2d_library}"
      TEST_COMMAND      ""
      # The git update step carries no stamp, so it is always out of date and
      # drags the configure, build and install steps with it on every `make`.
      # Disconnecting it pins the checkout to the tag already cloned; re-fetch
      # deliberately with the `2decomp-<version>-update` target.
      UPDATE_DISCONNECTED TRUE
    )

    # The module files land here during the build, but the directory has to exist
    # before it can be used as an include directory.
    file(MAKE_DIRECTORY "${decomp2d_install_dir}/include")

    # Stands in for the package config that is not installed yet, exporting what
    # it would export.  Its dependency on the external project is followed by
    # everything linking it, so the library is built first.
    add_library(decomp2d INTERFACE)
    target_include_directories(decomp2d INTERFACE "${decomp2d_install_dir}/include")
    target_link_libraries(decomp2d INTERFACE "${decomp2d_library}" MPI::MPI_Fortran)
    add_dependencies(decomp2d 2decomp-${decomp2d_version})
  endif(decomp2d_FOUND)

  # find_package() imports the target into this directory only, so promote it to
  # let sibling directories (e.g. tests/) link it.  The interface library built
  # above needs no promotion: it is global already.
  if(decomp2d_FOUND)
    get_target_property(_decomp2d_is_global decomp2d IMPORTED_GLOBAL)
    if(NOT _decomp2d_is_global)
      set_target_properties(decomp2d PROPERTIES IMPORTED_GLOBAL TRUE)
    endif()
  endif()
else()
  message(STATUS "2decomp-fft is disabled")
endif()
