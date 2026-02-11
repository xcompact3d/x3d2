/*
 * C wrapper for GPU-aware ADIOS2 put operation.
 *
 * Required because the Fortran generic adios2_put interface does not
 * accept CUDA device arrays. This wrapper receives the raw C handles
 * (engine%f2c, variable%f2c) from the Fortran side and calls the C API
 * adios2_put directly with the device pointer.
 *
 * Memory space must be set to GPU on the Fortran side via the native
 * adios2_set_memory_space binding BEFORE calling this wrapper.
 */

#include <adios2_c.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void adios2_put_gpu(
    const int64_t *engine_f2c,
    const int64_t *variable_f2c,
    const void *data,
    int *ierr)
{
    adios2_engine *engine = (adios2_engine *)(*engine_f2c);
    adios2_variable *variable = (adios2_variable *)(*variable_f2c);
    *ierr = (int)adios2_put(engine, variable, data, adios2_mode_sync);
}

#ifdef __cplusplus
}
#endif
