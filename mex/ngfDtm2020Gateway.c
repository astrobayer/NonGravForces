/* NGFDTM2020GATEWAY  MATLAB gateway to CNES's unmodified DTM2020 (operational: F10.7, Kp).
 *
 * Built by NonGravForces.buildMex into +NonGravForces/+env/private/dtm2020Gateway.<mexext>,
 * visible to NonGravForces.env.density and to nothing else. The model is called
 * through mex/ngfDtm2020Shim.F90, two bind(C) wrappers of this toolbox's own around
 * LECDTM and DTM3.
 *
 * THIS FILE CONTAINS NO MODEL. It marshals K epochs and returns the model's own
 * outputs in its own units (g/cm^3) and order; NonGravForces.env.density converts
 * them to SI and to const.speciesName.
 *
 * CALL
 *   [d, t] = dtm2020Gateway(coeffFile, doy, alt, lst, glat, glong, f107, f107bar, kp3, kp24)
 *
 *   coeffFile  char         full name of DTM_2020_F107_Kp.dat
 *   doy .. kp24  1xK double day of year, altitude [km], local time [rad], latitude and
 *                           longitude [rad], F10.7 at t - 24 h, mean F10.7 of the last
 *                           81 days, Kp delayed by 3 h, mean Kp of the last 24 h
 *   d          8xK double   ro, then H He O N2 O2 N [g/cm^3], then the mean molecular mass [g]
 *   t          2xK double   [exospheric temperature; temperature at altitude] K
 *
 * LECDTM is run on the first call and again whenever the coefficient file changes. EVERY
 * INPUT IS SCREENED FOR NaN/Inf HERE, and the altitude for the model's 120 km floor:
 * env.density screens first as well; this is the second line because the first one is
 * not in this file.
 *
 * Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
 */
#include "mex.h"

#include <math.h>
#include <string.h>

extern void ngf_dtm2020_init(const char *path, int n, int *status);
extern void ngf_dtm2020_calc(double day, double f107, double f107bar, double kp3, double kp24,
                             double alt, double hl, double lat, double lon,
                             double *tz, double *tinf, double *ro, double *d, double *wmm);

#define GW_ID "NonGravForces:dtm2020Gateway:"
#define NARG 10

static char lastFile[1024] = "";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char file[1024];
    const double *x[NARG];
    double *dOut, *tOut, tz, tinf, ro, d[6], wmm;
    mwSize K, k;
    int i, j, status;

    if (nrhs != NARG) {
        mexErrMsgIdAndTxt(GW_ID "inputArity",
            "dtm2020Gateway takes exactly %d arguments, got %d -> call it through NonGravForces.env.density", NARG, nrhs);
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt(GW_ID "inputArity", "dtm2020Gateway returns at most 2 outputs -> [d, t]");
    }
    if (!mxIsChar(prhs[0]) || mxGetString(prhs[0], file, sizeof file) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputPath", "coeffFile must be a char row shorter than %d -> pass the full name of DTM_2020_F107_Kp.dat", (int)sizeof file);
    }
    K = mxGetNumberOfElements(prhs[1]);
    for (i = 1; i < NARG; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt(GW_ID "inputClass", "argument %d must be a real double array -> cast before the call", i + 1);
        }
        if (mxGetNumberOfElements(prhs[i]) != K) {
            mexErrMsgIdAndTxt(GW_ID "inputShape",
                "argument %d has %d element(s), argument 2 has %d -> broadcast in MATLAB before the call",
                i + 1, (int)mxGetNumberOfElements(prhs[i]), (int)K);
        }
        x[i] = mxGetPr(prhs[i]);
        for (k = 0; k < K; k++) {
            if (!isfinite(x[i][k])) {
                mexErrMsgIdAndTxt(GW_ID "numericNaN", "argument %d holds NaN/Inf at element %d -> screen the inputs", i + 1, (int)k + 1);
            }
        }
    }
    for (k = 0; k < K; k++) {
        if (x[2][k] < 120.0) {
            mexErrMsgIdAndTxt(GW_ID "inputRange", "altitude %g km at element %d is below DTM2020's 120 km floor -> use another density model there",
                x[2][k], (int)k + 1);
        }
    }
    if (strcmp(file, lastFile) != 0) {
        ngf_dtm2020_init(file, (int)strlen(file), &status);
        if (status != 0) {
            mexErrMsgIdAndTxt(GW_ID "dataMissing", "cannot open %s -> run setupNonGravForces and NonGravForces.buildMex, which copies it beside the gateway", file);
        }
        strncpy(lastFile, file, sizeof lastFile - 1);
    }

    plhs[0] = mxCreateDoubleMatrix(8, K, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2, K, mxREAL);
    dOut = mxGetPr(plhs[0]);
    tOut = mxGetPr(plhs[1]);
    for (k = 0; k < K; k++) {
        ngf_dtm2020_calc(x[1][k], x[6][k], x[7][k], x[8][k], x[9][k], x[2][k], x[3][k], x[4][k], x[5][k],
                         &tz, &tinf, &ro, d, &wmm);
        dOut[8 * k] = ro;
        for (j = 0; j < 6; j++) {
            dOut[8 * k + 1 + j] = d[j];
        }
        dOut[8 * k + 7] = wmm;
        tOut[2 * k] = tinf;
        tOut[2 * k + 1] = tz;
    }
}
