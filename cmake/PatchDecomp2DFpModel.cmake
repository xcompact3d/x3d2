# Run with cmake -P from the decomp2d ExternalProject patch step.
# Adds the Intel fp-model to the base Fortran flags of 2decomp-fft, which
# otherwise replaces CMAKE_Fortran_FLAGS and keeps the Intel fast default.
# Inputs: SRC_DIR (decomp2d source tree) and FP_MODEL (the flag to add).
set(_file "${SRC_DIR}/cmake/compilers/D2D_flags_intel.cmake")
set(_anchor "set(D2D_FFLAGS_DEBUG ")
set(_begin "# x3d2: fp-model begin")
set(_end "# x3d2: fp-model end")

file(READ "${_file}" _content)
# Drop the block from an earlier run, so a changed FP_MODEL replaces it.
string(REGEX REPLACE "${_begin}\n.*${_end}\n" "" _content "${_content}")

string(FIND "${_content}" "${_anchor}" _pos)
if(_pos EQUAL -1)
  message(FATAL_ERROR "decomp2d fp-model patch: anchor not found in ${_file}")
endif()
# D2D_Compilers.cmake copies D2D_FFLAGS into the cached CMAKE_Fortran_FLAGS only
# while FLAGS_SET is false.  On a resumed build the cache already holds the old
# model, so replace it there too.
set(_block "${_begin}
set(D2D_FFLAGS \"\${D2D_FFLAGS} ${FP_MODEL}\")
if(FLAGS_SET)
  string(REGEX REPLACE \" ?-fp-model=[a-z0-9]+\" \"\" _x3d2_flags \"\${CMAKE_Fortran_FLAGS}\")
  set(CMAKE_Fortran_FLAGS \"\${_x3d2_flags} ${FP_MODEL}\" CACHE STRING
    \"Base FFLAGS for build\" FORCE)
endif()
${_end}
")
string(REPLACE "${_anchor}" "${_block}${_anchor}" _content "${_content}")
file(WRITE "${_file}" "${_content}")
