#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(checkboundary)(double, int, int, double);
extern void F77_NAME(dmmf)(double, int, int, double, int, int, double, double, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
extern void F77_NAME(sinkfill)(double, int, int, double, double, double, double);

static const R_FortranMethodDef FortranEntries[] = {
    {"checkboundary", (DL_FUNC) &F77_NAME(checkboundary),  4},
    {"dmmf",          (DL_FUNC) &F77_NAME(dmmf),          60},
    {"sinkfill",      (DL_FUNC) &F77_NAME(sinkfill),       7},
    {NULL, NULL, 0}
};

void R_init_DMMF(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
