/* mex.h -- SHIM, NOT MATLAB'S HEADER.
 *
 * Declares just enough of the MEX API for a compiler to check
 * mex/ngfKernelGateway.c on a machine that has no MATLAB. Nothing here is
 * implemented and nothing links: the check is `cc -c`, and it catches the
 * typos, the wrong argument counts and the undeclared identifiers that would
 * otherwise wait for the owner's host. A Tier B session cannot run MEX, but it
 * has no excuse for shipping C that does not compile.
 *
 * This file is on the include path ONLY for that check
 * (mex/test/shim). NonGravForces.buildMex never sees it.
 *
 * Created 2026-09-08 by Claude Opus 5 (Anthropic).
 */
#ifndef NGF_MEX_SHIM_H
#define NGF_MEX_SHIM_H

#include <stddef.h>

typedef size_t mwSize;
typedef size_t mwIndex;
typedef struct mxArray_tag mxArray;

typedef enum { mxREAL = 0, mxCOMPLEX = 1 } mxComplexity;
typedef enum { mxDOUBLE_CLASS = 6 } mxClassID;

#if defined(__GNUC__)
#define NGF_SHIM_NORETURN __attribute__((noreturn))
#else
#define NGF_SHIM_NORETURN
#endif

NGF_SHIM_NORETURN void mexErrMsgIdAndTxt(const char *id, const char *fmt, ...);

const mxArray *mxGetField(const mxArray *a, mwIndex i, const char *name);
void           mxSetField(mxArray *a, mwIndex i, const char *name, mxArray *v);
int            mxIsDouble(const mxArray *a);
int            mxIsComplex(const mxArray *a);
int            mxIsEmpty(const mxArray *a);
int            mxIsStruct(const mxArray *a);
mwSize         mxGetM(const mxArray *a);
mwSize         mxGetN(const mxArray *a);
mwSize         mxGetNumberOfElements(const mxArray *a);
double        *mxGetPr(const mxArray *a);
double         mxGetScalar(const mxArray *a);
mxArray       *mxCreateDoubleMatrix(mwSize m, mwSize n, mxComplexity c);
mxArray       *mxCreateNumericArray(mwSize nd, const mwSize *dims, mxClassID id, mxComplexity c);
mxArray       *mxCreateStructMatrix(mwSize m, mwSize n, int nfields, const char **names);
void           mxDestroyArray(mxArray *a);
void          *mxCalloc(mwSize n, mwSize size);

#endif /* NGF_MEX_SHIM_H */
