/* NGFMSISGATEWAY  MATLAB gateway to the vendored NRLMSISE-00 C (Brodowski).
 *
 * Built by NonGravForces.buildMex into +NonGravForces/+env/private/msisGateway.<mexext>,
 * visible to NonGravForces.env.density and to nothing else.
 *
 * THIS FILE CONTAINS NO MODEL. It marshals K epochs into the unmodified
 * public-domain implementation in mex/nrlmsise00/nrlmsise-00.c and marshals the
 * results back in the reference's own units (cgs), so that the published
 * DOCUMENTATION table can be compared without any arithmetic in between.
 * NonGravForces.env.density converts to SI.
 *
 * CALL
 *   [d, t] = msisGateway(doy, sec, alt, glat, glong, lst, f107a, f107, ap, anomalous)
 *
 *   doy .. f107  1xK double: day of year, UT seconds of day, altitude [km],
 *                geodetic latitude and longitude [deg], local solar time [h],
 *                81-day centred F10.7, previous-day F10.7.
 *   ap           1xK double  daily Ap             -> switches[9] = 1
 *                7xK double  NRLMSISE ap history  -> switches[9] = -1
 *   anomalous    1x1 double  non-zero: gtd7d (anomalous oxygen in the mass
 *                density), zero: gtd7.
 *   d            9xK double  cm^-3 for He O N2 O2 Ar H N O_anomalous (rows
 *                1-5, 7-9) and g/cm^3 total mass density (row 6): the
 *                reference's own d[0..8] order.
 *   t            2xK double  K: exospheric temperature, temperature at alt.
 *
 * Every argument is validated for class and size; MATLAB does all
 * broadcasting before the call, where it is visible. Identifiers follow the
 * NonGravForces taxonomy (input*).
 *
 * Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
 */
#include "mex.h"
#include "nrlmsise-00.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double *doy, *sec, *alt, *glat, *glong, *lst, *f107a, *f107, *ap;
    double *dOut, *tOut;
    double anomalous;
    mwSize K, k, apRows;
    int i, j;
    struct nrlmsise_input input;
    struct nrlmsise_output output;
    struct nrlmsise_flags flags;
    struct ap_array aph;

    if (nrhs != 10) {
        mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputArity",
            "msisGateway takes exactly 10 arguments, got %d -> call it through NonGravForces.env.density", nrhs);
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputArity",
            "msisGateway returns at most 2 outputs -> [d, t]");
    }
    for (i = 0; i < 10; i++) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputClass",
                "argument %d must be a real double array -> cast before the call", i + 1);
        }
    }
    K = mxGetNumberOfElements(prhs[0]);
    for (i = 1; i < 8; i++) {
        if (mxGetNumberOfElements(prhs[i]) != K) {
            mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputShape",
                "argument %d has %d element(s), argument 1 has %d -> broadcast in MATLAB before the call",
                i + 1, (int)mxGetNumberOfElements(prhs[i]), (int)K);
        }
    }
    apRows = mxGetM(prhs[8]);
    if (!(apRows == 1 || apRows == 7) || mxGetN(prhs[8]) != K) {
        mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputShape",
            "ap must be 1x%d (daily) or 7x%d (history), got %dx%d -> see NonGravForces.io.spaceWeather query",
            (int)K, (int)K, (int)apRows, (int)mxGetN(prhs[8]));
    }
    if (mxGetNumberOfElements(prhs[9]) != 1) {
        mexErrMsgIdAndTxt("NonGravForces:msisGateway:inputShape",
            "anomalous must be a scalar -> pass 0 (gtd7) or 1 (gtd7d)");
    }

    doy   = mxGetPr(prhs[0]);
    sec   = mxGetPr(prhs[1]);
    alt   = mxGetPr(prhs[2]);
    glat  = mxGetPr(prhs[3]);
    glong = mxGetPr(prhs[4]);
    lst   = mxGetPr(prhs[5]);
    f107a = mxGetPr(prhs[6]);
    f107  = mxGetPr(prhs[7]);
    ap    = mxGetPr(prhs[8]);
    anomalous = mxGetScalar(prhs[9]);

    plhs[0] = mxCreateDoubleMatrix(9, K, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2, K, mxREAL);
    dOut = mxGetPr(plhs[0]);
    tOut = mxGetPr(plhs[1]);

    flags.switches[0] = 0;                 /* cgs output, the reference's default */
    for (i = 1; i < 24; i++) {
        flags.switches[i] = 1;
    }
    flags.switches[9] = (apRows == 7) ? -1 : 1;

    for (k = 0; k < K; k++) {
        input.year   = 0;                  /* ignored by the model */
        input.doy    = (int)doy[k];
        input.sec    = sec[k];
        input.alt    = alt[k];
        input.g_lat  = glat[k];
        input.g_long = glong[k];
        input.lst    = lst[k];
        input.f107A  = f107a[k];
        input.f107   = f107[k];
        if (apRows == 7) {
            for (j = 0; j < 7; j++) {
                aph.a[j] = ap[7 * k + j];
            }
            input.ap   = aph.a[0];
            input.ap_a = &aph;
        } else {
            input.ap   = ap[k];
            input.ap_a = NULL;
        }
        if (anomalous != 0.0) {
            gtd7d(&input, &flags, &output);
        } else {
            gtd7(&input, &flags, &output);
        }
        for (j = 0; j < 9; j++) {
            dOut[9 * k + j] = output.d[j];
        }
        tOut[2 * k]     = output.t[0];
        tOut[2 * k + 1] = output.t[1];
    }
}
