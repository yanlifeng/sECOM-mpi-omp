/* Initialize the *same* I/O library before gfortran main enables FP traps.
   No change to solver flags, data, or arithmetic. Both modes use this shim. */
#include <hdf5.h>
#include <stdlib.h>
#include <fenv.h>
__attribute__((constructor)) static void secom_h5_preinit(void) {
    fenv_t env;
    if (feholdexcept(&env) != 0 || H5open() < 0 || fesetenv(&env) != 0)
        abort();
}
