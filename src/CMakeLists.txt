file(GLOB SRC
  LIST_DIRECTORIES false
  *.f90
)
add_library(x3d2 STATIC ${SRC})
