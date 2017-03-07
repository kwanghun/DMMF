#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(checkboundary)(void *, void *, void *, void *);
extern void F77_NAME(dmmf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sinkfill)(void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"checkboundary", (DL_FUNC) &F77_NAME(checkboundary),  4},
    {"dmmf",          (DL_FUNC) &F77_NAME(dmmf),          60},
    {"sinkfill",      (DL_FUNC) &F77_NAME(sinkfill),       7},
    {NULL, NULL, 0}
};

void R_init_DMMF(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
