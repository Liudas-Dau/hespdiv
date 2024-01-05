#include <R.h>
#include <Rinternals.h>
#include "hespdiv.h"

#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[]  = {
/*    {"pipbb", (DL_FUNC) &pipbb, 3},
    {"between", (DL_FUNC) &between, 3}, */
    {"setup_poly_minmax", (DL_FUNC) &setup_poly_minmax, 1},
    {"InPoly", (DL_FUNC) &InPoly, 2},
/* RSB 091203 */
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
/*    {"insiders", (DL_FUNC) &insiders, 2}, */
    {"R_point_in_polygon_sp", (DL_FUNC) &R_point_in_polygon_sp, 4},
/* RSB 091203 */
    {NULL, NULL, 0}
};


void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_hespdiv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

}
