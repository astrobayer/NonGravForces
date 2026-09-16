/* NGFMSIS21GATEWAY  MATLAB gateway to the unmodified NRLMSIS 2.1 Fortran (NRL).
 *
 * Built by NonGravForces.buildMex into +NonGravForces/+env/private/msis21Gateway.<mexext>,
 * visible to NonGravForces.env.density and to nothing else. The model is called
 * through mex/ngfMsis21Shim.F90, two bind(C) wrappers of this toolbox's own
 * around MSISINIT and MSISCALC.
 *
 * THIS FILE CONTAINS NO MODEL. It marshals K epochs and returns the model's own
 * SI outputs in its own species order; NonGravForces.env.density reorders them
 * into const.speciesName and replaces the model's "missing" marker.
 *
 * CALL
 *   [d, t] = msis21Gateway(parmDir, anomalous, doy, sec, alt, glat, glong, f107a, f107, ap)
 *
 *   parmDir    char         directory of msis21.parm WITH a trailing separator; with the
 *                           file name at most 256 characters (see the shim)
 *   anomalous  1x1 double   nonzero: the mass density includes anomalous oxygen
 *   doy .. f107  1xK double day of year, UT seconds, geodetic altitude [km], geodetic
 *                           latitude and longitude [deg], 81-day centred F10.7, previous-day F10.7
 *   ap         1xK double   daily Ap                 -> legacy switch 9 = 1
 *              7xK double   NRLMSISE ap history      -> legacy switch 9 = -1
 *   d          10xK double  rho [kg/m^3], then N2 O2 O He H Ar N anomalous-O NO [m^-3]
 *   t          2xK double   [exospheric temperature; temperature at altitude] K
 *
 * MSISINIT is run on the first call and again whenever the parameter directory, the
 * ap mode or the mass-density flag changes. EVERY INPUT IS SCREENED FOR NaN/Inf
 * HERE, before the model sees it: a NaN altitude reaches a STOP in msis_dfn.F90
 * (pwmp), and a Fortran STOP ends the MATLAB process, not the call. env.density
 * screens first as well; this is the second line because the first one is not
 * in this file.
 *
 * Created 2026-09-15 by Claude Opus 5 (Anthropic).
 */
#include "mex.h"

#include <math.h>
#include <string.h>

extern void ngf_msis21_init(const char *dir, int n, int aphistory, int anomalous, int *status);
extern void ngf_msis21_calc(double day, double utsec, double z, double lat, double lon,
                            double f107a, double f107, const double *ap,
                            double *tn, double *dn, double *tex);

#define GW_ID "NonGravForces:msis21Gateway:"
#define NARG 10

static char lastDir[512] = "";
static int lastMode = -99;
static int lastAnomalous = -99;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char dir[512];
    const double *doy, *sec, *alt, *glat, *glong, *f107a, *f107, *ap;
    double *dOut, *tOut, apk[7], tn, tex, dn[10];
    mwSize K, k, apRows, nAll;
    int i, j, mode, anomalous, status;

    if (nrhs != NARG) {
        mexErrMsgIdAndTxt(GW_ID "inputArity",
            "msis21Gateway takes exactly %d arguments, got %d -> call it through NonGravForces.env.density", NARG, nrhs);
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt(GW_ID "inputArity", "msis21Gateway returns at most 2 outputs -> [d, t]");
    }
    if (!mxIsChar(prhs[0]) || mxGetString(prhs[0], dir, sizeof dir) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputPath", "parmDir must be a char row shorter than %d -> pass the folder of msis21.parm", (int)sizeof dir);
    }
    for (i = 1; i < NARG; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(GW_ID "inputClass", "argument %d must be a real double array -> cast before the call", i + 1);
        }
    }
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputShape", "anomalous must be a scalar -> pass 1 or 0");
    }
    K = mxGetNumberOfElements(prhs[2]);
    for (i = 3; i < NARG - 1; i++) {
        if (mxGetNumberOfElements(prhs[i]) != K) {
            mexErrMsgIdAndTxt(GW_ID "inputShape",
                "argument %d has %d element(s), argument 3 has %d -> broadcast in MATLAB before the call",
                i + 1, (int)mxGetNumberOfElements(prhs[i]), (int)K);
        }
    }
    apRows = mxGetM(prhs[NARG - 1]);
    if (!(apRows == 1 || apRows == 7) || mxGetN(prhs[NARG - 1]) != K) {
        mexErrMsgIdAndTxt(GW_ID "inputShape", "ap must be 1x%d (daily) or 7x%d (history) -> got %dx%d",
            (int)K, (int)K, (int)apRows, (int)mxGetN(prhs[NARG - 1]));
    }
    for (i = 1; i < NARG; i++) {
        const double *x = mxGetPr(prhs[i]);
        nAll = mxGetNumberOfElements(prhs[i]);
        for (k = 0; k < nAll; k++) {
            if (!isfinite(x[k])) {
                mexErrMsgIdAndTxt(GW_ID "numericNaN",
                    "argument %d holds NaN/Inf at element %d -> NRLMSIS 2.1 would STOP the MATLAB process on it; screen the inputs",
                    i + 1, (int)k + 1);
            }
        }
    }
    mode = (apRows == 7) ? 1 : 0;
    anomalous = (mxGetScalar(prhs[1]) != 0.0) ? 1 : 0;
    if (strcmp(dir, lastDir) != 0 || mode != lastMode || anomalous != lastAnomalous) {
        ngf_msis21_init(dir, (int)strlen(dir), mode, anomalous, &status);
        if (status == 3) {
            mexErrMsgIdAndTxt(GW_ID "dataMissing", "msis21.parm is not in %s -> run NonGravForces.buildMex, which copies it beside the gateway", dir);
        }
        if (status != 0) {
            mexErrMsgIdAndTxt(GW_ID "inputPath",
                "the parameter file name is longer than MSISINIT can hold (%d characters in the folder) -> install the toolbox under a shorter path",
                (int)strlen(dir));
        }
        strncpy(lastDir, dir, sizeof lastDir - 1);
        lastMode = mode;
        lastAnomalous = anomalous;
    }

    doy = mxGetPr(prhs[2]); sec = mxGetPr(prhs[3]); alt = mxGetPr(prhs[4]);
    glat = mxGetPr(prhs[5]); glong = mxGetPr(prhs[6]);
    f107a = mxGetPr(prhs[7]); f107 = mxGetPr(prhs[8]); ap = mxGetPr(prhs[9]);

    plhs[0] = mxCreateDoubleMatrix(10, K, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2, K, mxREAL);
    dOut = mxGetPr(plhs[0]);
    tOut = mxGetPr(plhs[1]);

    for (k = 0; k < K; k++) {
        for (j = 0; j < 7; j++) {
            apk[j] = (apRows == 7) ? ap[7 * k + j] : ap[k];
        }
        ngf_msis21_calc(doy[k], sec[k], alt[k], glat[k], glong[k], f107a[k], f107[k], apk, &tn, dn, &tex);
        for (j = 0; j < 10; j++) {
            dOut[10 * k + j] = dn[j];
        }
        tOut[2 * k] = tex;
        tOut[2 * k + 1] = tn;
    }
}
