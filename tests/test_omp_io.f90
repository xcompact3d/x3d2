program test_io
  use iso_fortran_env, only: stderr => error_unit

  implicit none

  logical :: allpass

  allpass = .true.

  print *, "Hello World!"
end program
