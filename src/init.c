#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ADMMSLOPE_ADMM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ADMMSLOPE_FastProxSL1(SEXP, SEXP);
extern SEXP _ADMMSLOPE_rcpparma_bothproducts(SEXP);
extern SEXP _ADMMSLOPE_rcpparma_hello_world();
extern SEXP _ADMMSLOPE_rcpparma_innerproduct(SEXP);
extern SEXP _ADMMSLOPE_rcpparma_outerproduct(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ADMMSLOPE_ADMM",                  (DL_FUNC) &_ADMMSLOPE_ADMM,                  5},
    {"_ADMMSLOPE_FastProxSL1",           (DL_FUNC) &_ADMMSLOPE_FastProxSL1,           2},
    {"_ADMMSLOPE_rcpparma_bothproducts", (DL_FUNC) &_ADMMSLOPE_rcpparma_bothproducts, 1},
    {"_ADMMSLOPE_rcpparma_hello_world",  (DL_FUNC) &_ADMMSLOPE_rcpparma_hello_world,  0},
    {"_ADMMSLOPE_rcpparma_innerproduct", (DL_FUNC) &_ADMMSLOPE_rcpparma_innerproduct, 1},
    {"_ADMMSLOPE_rcpparma_outerproduct", (DL_FUNC) &_ADMMSLOPE_rcpparma_outerproduct, 1},
    {NULL, NULL, 0}
};

void R_init_ADMMSLOPE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}