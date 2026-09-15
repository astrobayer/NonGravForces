/* NGFHWM14GATEWAY  MATLAB gateway to the unmodified HWM14 Fortran (NRL).
 *
 * Built by NonGravForces.buildMex into +NonGravForces/+env/private/hwm14Gateway.<mexext>,
 * visible to NonGravForces.env.wind and to nothing else. The model is called through
 * mex/ngfHwm14Shim.F90, a bind(C) wrapper of this toolbox's own around HWM14.
 *
 * THIS FILE CONTAINS NO MODEL. It marshals K points and returns the model's two wind
 * components; env.wind turns them into a GCRS vector.
 *
 * CALL
 *   w = hwm14Gateway(dataDir, iyd, sec, alt, glat, glon, ap2)
 *
 *   dataDir  char        folder holding hwm123114.bin, dwm07b104i.dat and gd2qd.dat
 *   iyd      1xK double  yyddd (integer valued)
 *   sec      1xK double  UT seconds
 *   alt      1xK double  geodetic altitude [km]
 *   glat     1xK double  geodetic latitude [deg]
 *   glon     1xK double  geodetic longitude [deg]
 *   ap2      1xK double  current 3-hour ap; negative for the quiet-time wind alone
 *   w        2xK double  [meridional (+north); zonal (+east)] [m/s]
 *
 * TWO THINGS IN THE NRL CODE THIS FILE HAS TO WORK AROUND, without editing it:
 *   * The data files are found by FINDANDOPEN, which looks in the current directory,
 *     then in HWMPATH (a character(128)), then in ../Meta, and executes STOP when it
 *     finds nothing -- which ends the MATLAB process. So this gateway checks that the
 *     three files exist, makes dataDir the PROCESS current directory for the one call
 *     that loads them (no length limit, no environment variable), and restores the
 *     previous directory before it returns, whatever happened (mex/ngfHwm14Load.c,
 *     shared with the native right-hand side).
 *   * The model loads its tables once, on the first call, and keeps them for the
 *     life of the process. A different dataDir afterwards is therefore refused, not
 *     silently ignored.
 *   * HWM14 NEVER RETURNS for a negative altitude (measured 2026-09-15: -1e-6 km and -1 km
 *     both ran into a 10 s timeout in a standalone build; 742 km and 5000 km answered).
 *     A negative altitude is therefore refused here, before the model sees it.
 *   Every input is screened for NaN/Inf here, before the model sees it.
 *
 * Created 2026-09-15 by Claude Opus 5 (Anthropic).
 */
#include "mex.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

extern void ngf_hwm14_calc(int iyd, double sec, double alt, double glat, double glon, double ap2, double *w);
extern int ngfHwm14Load(const char *dir, char *missing, int missingSize);   /* mex/ngfHwm14Load.c */

#define GW_ID "NonGravForces:hwm14Gateway:"
#define NARG 7

static char loadedDir[1024] = "";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char dir[1024];
    const double *iyd, *sec, *alt, *glat, *glon, *ap2;
    double *out, w[2];
    mwSize K, k, nAll;
    int i;

    if (nrhs != NARG) {
        mexErrMsgIdAndTxt(GW_ID "inputArity", "hwm14Gateway takes exactly %d arguments, got %d -> call it through NonGravForces.env.wind", NARG, nrhs);
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt(GW_ID "inputArity", "hwm14Gateway returns one output -> w");
    }
    if (!mxIsChar(prhs[0]) || mxGetString(prhs[0], dir, sizeof dir) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputPath", "dataDir must be a char row shorter than %d -> pass the folder of the HWM14 data files", (int)sizeof dir);
    }
    for (i = 1; i < NARG; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(GW_ID "inputClass", "argument %d must be a real double array -> cast before the call", i + 1);
        }
    }
    K = mxGetNumberOfElements(prhs[1]);
    for (i = 2; i < NARG; i++) {
        if (mxGetNumberOfElements(prhs[i]) != K) {
            mexErrMsgIdAndTxt(GW_ID "inputShape", "argument %d has %d element(s), argument 2 has %d -> broadcast in MATLAB before the call",
                i + 1, (int)mxGetNumberOfElements(prhs[i]), (int)K);
        }
    }
    for (i = 1; i < NARG; i++) {
        const double *x = mxGetPr(prhs[i]);
        nAll = mxGetNumberOfElements(prhs[i]);
        for (k = 0; k < nAll; k++) {
            if (!isfinite(x[k])) {
                mexErrMsgIdAndTxt(GW_ID "numericNaN", "argument %d holds NaN/Inf at element %d -> screen the inputs", i + 1, (int)k + 1);
            }
        }
    }
    {
        const double *h = mxGetPr(prhs[3]);
        for (k = 0; k < K; k++) {
            if (!(h[k] >= 0.0)) {
                mexErrMsgIdAndTxt(GW_ID "inputRange", "altitude %g km at element %d is below the ellipsoid -> HWM14 does not return for a negative altitude; screen the positions", h[k], (int)k + 1);
            }
        }
    }
    if (loadedDir[0] == '\0') {
        char missing[64] = "";
        int loaded = ngfHwm14Load(dir, missing, (int)sizeof missing);
        if (loaded == 1) {
            mexErrMsgIdAndTxt(GW_ID "dataMissing", "%s is not in %s -> run NonGravForces.buildMex, which copies the HWM14 data beside the gateway", missing, dir);
        }
        if (loaded == 2) {
            mexErrMsgIdAndTxt(GW_ID "inputPath", "cannot enter %s, or read the current directory -> check the folder of the HWM14 data files", dir);
        }
        strncpy(loadedDir, dir, sizeof loadedDir - 1);
        if (loaded == 3) {
            mexErrMsgIdAndTxt(GW_ID "inputPath", "HWM14 was loaded but the previous directory could not be restored -> cd back by hand");
        }
    } else if (strcmp(loadedDir, dir) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputPath", "HWM14 was loaded from %s and keeps its tables for the life of the process -> clear the gateway (clear mex) before loading it from %s", loadedDir, dir);
    }

    iyd = mxGetPr(prhs[1]); sec = mxGetPr(prhs[2]); alt = mxGetPr(prhs[3]);
    glat = mxGetPr(prhs[4]); glon = mxGetPr(prhs[5]); ap2 = mxGetPr(prhs[6]);
    plhs[0] = mxCreateDoubleMatrix(2, K, mxREAL);
    out = mxGetPr(plhs[0]);
    for (k = 0; k < K; k++) {
        ngf_hwm14_calc((int)iyd[k], sec[k], alt[k], glat[k], glon[k], ap2[k], w);
        out[2 * k] = w[0];
        out[2 * k + 1] = w[1];
    }
}
