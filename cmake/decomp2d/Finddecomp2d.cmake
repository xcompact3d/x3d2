# - Find the 2decomp-fft library
set(decomp2d_git_tag "v2.0.3")
string(REPLACE "/" "-" decomp2d_git_tag_dir "${decomp2d_git_tag}")

if(SINGLE_PREC)
  set(decomp2d_install_dir "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt-sp-${decomp2d_git_tag_dir}")
else()
  set(decomp2d_install_dir "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-opt-dp-${decomp2d_git_tag_dir}")
endif()

# find_package() caches decomp2d_DIR, so a previous precision can otherwise
# satisfy this lookup even when decomp2d_install_dir has changed.
unset(decomp2d_DIR CACHE)
unset(decomp2d_DIR)

find_package(decomp2d CONFIG
             PATHS ${decomp2d_install_dir}
             NO_DEFAULT_PATH
             QUIET)
if (decomp2d_FOUND)
  message(STATUS "2decomp-fft FOUND in ${decomp2d_install_dir}")
else(decomp2d_FOUND)
  message(STATUS "2decomp-fft PATH not available we'll try to download and install")

  # Keep the driver project separate from the external project's BINARY_DIR.
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download")

  if (SINGLE_PREC)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/downloadBuild2decompSinglePrec.cmake.in decomp2d-download/CMakeLists.txt)
  else()
    configure_file(${CMAKE_CURRENT_LIST_DIR}/downloadBuild2decomp.cmake.in decomp2d-download/CMakeLists.txt)
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download )
  if(result)
      message(FATAL_ERROR "CMake step for 2decomp-fft failed: ${result}")
  else()
      message("CMake step for 2decomp-fft completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
         RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-download )
  if(result)
      message(FATAL_ERROR "Build step for 2decomp-fft failed: ${result}")
  endif()
  set(D2D_ROOT ${decomp2d_install_dir})

  find_package(decomp2d REQUIRED CONFIG
          PATHS ${D2D_ROOT}
          NO_DEFAULT_PATH)
endif(decomp2d_FOUND)

# Make the imported target global so sibling directories (e.g. tests/) can
# link it.
if(TARGET decomp2d)
  get_target_property(_decomp2d_is_global decomp2d IMPORTED_GLOBAL)
  if(NOT _decomp2d_is_global)
    set_target_properties(decomp2d PROPERTIES IMPORTED_GLOBAL TRUE)
  endif()
endif()
