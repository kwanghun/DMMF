#include "DMMF.h"
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Fortran calls */
static const R_FortranMethodDef FortranEntries[] = {
    {"dmmf",          (DL_FUNC) &F77_SUB(dmmf),          60},
    {"ke",            (DL_FUNC) &F77_SUB(ke),             5},
    {"mdinf",         (DL_FUNC) &F77_SUB(mdinf),          1},
    {"slope",         (DL_FUNC) &F77_SUB(slope),          3},
    {"sinkfill",      (DL_FUNC) &F77_SUB(sinkfill),       7},
    {"checkboundary", (DL_FUNC) &F77_SUB(checkboundary),  4},
    {NULL, NULL, 0}
};

void R_init_DMMF(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
