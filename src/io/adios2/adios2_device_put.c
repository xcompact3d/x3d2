/*
 * GPU-aware ADIOS2 put for the Fortran I/O backend.
 *
 * The Fortran generic adios2_put does not accept CUDA device arrays, so
 * m_io_backend calls the C API instead. ADIOS2's Fortran bindings store each
 * C handle as an integer(8) (the `f2c` component) holding the handle's
 * address. Converting it back to a pointer here, rather than reinterpreting
 * it as a c_ptr in Fortran, keeps the Fortran interface ISO C interoperable:
 * the handles cross as int64_t values and the data as a void pointer. This is
 * the same conversion ADIOS2's own Fortran-to-C layer performs.
 *
 * The variable's memory space must already be set to GPU with the Fortran
 * adios2_set_memory_space before calling this.
 */

#include <adios2_c.h>
#include <stdint.h>

_Static_assert(sizeof(void *) <= sizeof(int64_t),
               "ADIOS2 Fortran f2c handles must be able to hold a pointer");

int x3d2_adios2_put_device(int64_t engine_f2c, int64_t variable_f2c,
                           const void *data, int mode)
{
  adios2_engine *engine = (adios2_engine *)(intptr_t)engine_f2c;
  adios2_variable *variable = (adios2_variable *)(intptr_t)variable_f2c;

  return (int)adios2_put(engine, variable, data, (adios2_mode)mode);
}
