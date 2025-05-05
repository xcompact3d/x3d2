# Defines the backend selection for x3d2.
#
# Currently available backends
# - GPU/CUDAFortran (nvfort only)
# - CPU/OMP (all)

set(X3D2_BACKEND
  "OMP/CPU"
  CACHE
  STRING
  "Select the backend. Current options are 'OMP/CPU' and 'CUDA/GPU'")
set_property(CACHE
  X3D2_BACKEND
  PROPERTY STRINGS
  "OMP/CPU" "CUDA/GPU")
