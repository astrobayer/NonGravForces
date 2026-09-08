/* NGFKERNELGATEWAY  MATLAB gateway to the single-epoch NonGravForces kernel.
 *
 * Built by NonGravForces.buildMex into
 * +NonGravForces/+kernel/private/kernelGateway.<mexext>, visible to
 * NonGravForces.kernel.singleEpoch and to nothing else.
 *
 * THIS FILE CONTAINS NO PHYSICS. It marshals one epoch into mex/ngfKernel.c and
 * marshals the results back. Every number it needs -- the speed of light, the
 * Stefan-Boltzmann and Boltzmann constants, the molecular masses -- arrives from
 * MATLAB, where NonGravForces.const and referenceSystem.iersConstants own them
 * (rule R1): a constant written twice is a constant that will disagree once.
 *
 * CALL
 *   [F, tau, Fk, aux] = kernelGateway(macro, gas, rad, tPanel, scale, con)
 *
 *   macro   struct: A 1xN, N 3xN, C 3xN, optVis 3xN, optIR 3xN, rcm 3x1
 *   gas     struct: rho 1x1, species 8x1 (or []), T_K 1x1, vrel 3x1,
 *                   model 1x1 (0 sentman, 1 maxwell, 2 constantCd),
 *                   alpha 1x1, Tw_K 1x1, Cd 1x1, SpecularFraction 1x1,
 *                   shadow Nx1 (or [] for all 1)
 *   rad     struct: dir 3xM, P Mx1, band 1xM (0 vis, 1 ir), shadow NxM (or [])
 *                   M = 0 is legal and gives no radiation force
 *   tPanel  Nx1 double, or [] for no thermal recoil
 *   scale   2x1 double [cdScale; crScale]
 *   con     struct: cLight 1x1, sigmaSB 1x1, kB 1x1, speciesMass 8x1 [kg per
 *                   MOLECULE, i.e. const().speciesMass / const().NA]
 *
 *   F, tau  3x3 double, columns drag | radiation | thermal, body frame
 *   Fk      3xNx3 double, per-panel forces, pages in the same order
 *   aux     struct: gamma, theta_rad, Cp, Ctau, absorbed_W_m2, emitted_W_m2
 *                   (each Nx1) and speedRatio (8x1)
 *
 * Shapes are checked here and NOT re-checked in MATLAB, and the ranges are
 * checked in MATLAB and NOT re-checked here beyond what the C returns: nothing
 * is validated twice (HANDOVER E.2 point 3).
 *
 * Created 2026-09-08 by Claude Opus 5 (Anthropic).
 */
#include "mex.h"
#include "ngfKernel.h"

#include <string.h>

#define GW_ID "NonGravForces:kernelGateway:"

static const double *requireField(const mxArray *s, const char *group,
                                  const char *name, mwSize rows, mwSize cols,
                                  int allowEmpty)
{
    const mxArray *f = mxGetField(s, 0, name);
    if (f == NULL) {
        mexErrMsgIdAndTxt(GW_ID "inputField",
            "%s lacks the field '%s' -> build the argument with NonGravForces.kernel.singleEpoch",
            group, name);
    }
    if (!mxIsDouble(f) || mxIsComplex(f)) {
        mexErrMsgIdAndTxt(GW_ID "inputClass",
            "%s.%s must be a real double array -> cast before the call", group, name);
    }
    if (allowEmpty && mxIsEmpty(f)) {
        return NULL;
    }
    if (mxGetM(f) != rows || mxGetN(f) != cols) {
        mexErrMsgIdAndTxt(GW_ID "inputShape",
            "%s.%s must be %dx%d, got %dx%d -> check the plate or source count",
            group, name, (int)rows, (int)cols, (int)mxGetM(f), (int)mxGetN(f));
    }
    return mxGetPr(f);
}

static double requireScalar(const mxArray *s, const char *group, const char *name)
{
    const mxArray *f = mxGetField(s, 0, name);
    if (f == NULL || !mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputField",
            "%s.%s must be a real double scalar -> check the argument", group, name);
    }
    return mxGetScalar(f);
}

static void requireStruct(const mxArray *a, const char *name)
{
    if (!mxIsStruct(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputClass",
            "%s must be a 1x1 struct -> see the help of NonGravForces.kernel.singleEpoch", name);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *macroArg, *gasArg, *radArg, *conArg, *field;
    ngfMacro macro;
    ngfConstants con;
    ngfSources src;
    ngfEpochInput in;
    ngfEpochOutput out;
    const double *speciesMass, *scale, *band;
    double *Fk;
    int *bandInt = NULL;
    mwSize n, m, j;
    mwSize dims[3];
    int status;
    const char *auxNames[7] = {"gamma", "theta_rad", "speedRatio", "Cp", "Ctau",
                               "absorbed_W_m2", "emitted_W_m2"};
    mxArray *aux, *tmp;

    if (nrhs != 6) {
        mexErrMsgIdAndTxt(GW_ID "inputArity",
            "kernelGateway takes exactly 6 arguments, got %d -> call it through NonGravForces.kernel.singleEpoch", nrhs);
    }
    if (nlhs > 4) {
        mexErrMsgIdAndTxt(GW_ID "inputArity",
            "kernelGateway returns at most 4 outputs -> [F, tau, Fk, aux]");
    }
    macroArg = prhs[0]; gasArg = prhs[1]; radArg = prhs[2]; conArg = prhs[5];
    requireStruct(macroArg, "macro");
    requireStruct(gasArg, "gas");
    requireStruct(radArg, "rad");
    requireStruct(conArg, "con");

    /* ---- plate count from macro.A, then every other field against it ---- */
    field = mxGetField(macroArg, 0, "A");
    if (field == NULL || !mxIsDouble(field) || mxIsComplex(field) || mxGetM(field) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputField", "macro.A must be a real 1xN double -> the plate areas");
    }
    n = mxGetN(field);
    if (n < 1) {
        mexErrMsgIdAndTxt(GW_ID "inputShape", "macro.A is empty -> a macro model needs at least one plate");
    }
    macro.n = (int)n;
    macro.area = mxGetPr(field);
    macro.normal   = requireField(macroArg, "macro", "N", 3, n, 0);
    macro.centroid = requireField(macroArg, "macro", "C", 3, n, 0);
    macro.optVis   = requireField(macroArg, "macro", "optVis", 3, n, 0);
    macro.optIR    = requireField(macroArg, "macro", "optIR", 3, n, 0);
    in.rcm         = requireField(macroArg, "macro", "rcm", 3, 1, 0);

    /* ---- constants, all of them from MATLAB ---- */
    con.cLight  = requireScalar(conArg, "con", "cLight");
    con.sigmaSB = requireScalar(conArg, "con", "sigmaSB");
    con.kB      = requireScalar(conArg, "con", "kB");
    speciesMass = requireField(conArg, "con", "speciesMass", NGF_SPECIES, 1, 0);
    for (j = 0; j < NGF_SPECIES; j++) {
        con.speciesMass[j] = speciesMass[j];
    }

    /* ---- gas ---- */
    in.rho      = requireScalar(gasArg, "gas", "rho");
    in.tGas     = requireScalar(gasArg, "gas", "T_K");
    in.alpha    = requireScalar(gasArg, "gas", "alpha");
    in.tw       = requireScalar(gasArg, "gas", "Tw_K");
    in.cd       = requireScalar(gasArg, "gas", "Cd");
    in.specular = requireScalar(gasArg, "gas", "SpecularFraction");
    in.gsiModel = (ngfGsiModel)((int)requireScalar(gasArg, "gas", "model"));
    in.vrel      = requireField(gasArg, "gas", "vrel", 3, 1, 0);
    in.species   = requireField(gasArg, "gas", "species", NGF_SPECIES, 1, 1);
    in.shadowGas = requireField(gasArg, "gas", "shadow", n, 1, 1);

    /* ---- radiation: M = 0 is legal ---- */
    field = mxGetField(radArg, 0, "dir");
    if (field == NULL || !mxIsDouble(field) || mxIsComplex(field) ||
        (!mxIsEmpty(field) && mxGetM(field) != 3)) {
        mexErrMsgIdAndTxt(GW_ID "inputField", "rad.dir must be a real 3xM double -> unit vectors to the sources");
    }
    m = mxIsEmpty(field) ? 0 : mxGetN(field);
    src.m = (int)m;
    src.dir  = (m == 0) ? NULL : mxGetPr(field);
    src.pres = (m == 0) ? NULL : requireField(radArg, "rad", "P", m, 1, 0);
    band     = (m == 0) ? NULL : requireField(radArg, "rad", "band", 1, m, 0);
    if (m > 0) {
        bandInt = (int *)mxCalloc(m, sizeof(int));
        for (j = 0; j < m; j++) {
            bandInt[j] = (int)band[j];
        }
    }
    src.band = bandInt;
    in.src = &src;
    in.shadowRad = (m == 0) ? NULL : requireField(radArg, "rad", "shadow", n, m, 1);

    /* ---- thermal, scales ---- */
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt(GW_ID "inputClass", "tPanel must be a real double array -> Nx1 or []");
    }
    if (mxIsEmpty(prhs[3])) {
        in.tPanel = NULL;
    } else {
        if (mxGetM(prhs[3]) != n || mxGetN(prhs[3]) != 1) {
            mexErrMsgIdAndTxt(GW_ID "inputShape", "tPanel must be %dx1 or [], got %dx%d",
                (int)n, (int)mxGetM(prhs[3]), (int)mxGetN(prhs[3]));
        }
        in.tPanel = mxGetPr(prhs[3]);
    }
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 2) {
        mexErrMsgIdAndTxt(GW_ID "inputShape", "scale must be a real 2x1 double [cdScale; crScale]");
    }
    scale = mxGetPr(prhs[4]);
    in.cdScale = scale[0];
    in.crScale = scale[1];
    in.macro = &macro;
    in.con = &con;

    /* ---- outputs ---- */
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
    dims[0] = 3; dims[1] = n; dims[2] = 3;
    Fk = (double *)mxCalloc(3 * n * 3, sizeof(double));
    memset(&out, 0, sizeof out);
    out.F = mxGetPr(plhs[0]);
    out.tau = mxGetPr(plhs[1]);
    out.Fk = Fk;

    aux = mxCreateStructMatrix(1, 1, 7, auxNames);
    for (j = 0; j < 7; j++) {
        mwSize len = (j == 2) ? (mwSize)NGF_SPECIES : n;
        tmp = mxCreateDoubleMatrix(len, 1, mxREAL);
        mxSetField(aux, 0, auxNames[j], tmp);
    }
    out.gamma    = mxGetPr(mxGetField(aux, 0, "gamma"));
    out.theta    = mxGetPr(mxGetField(aux, 0, "theta_rad"));
    out.sRatio   = mxGetPr(mxGetField(aux, 0, "speedRatio"));
    out.cpMix    = mxGetPr(mxGetField(aux, 0, "Cp"));
    out.ctauMix  = mxGetPr(mxGetField(aux, 0, "Ctau"));
    out.absorbed = mxGetPr(mxGetField(aux, 0, "absorbed_W_m2"));
    out.emitted  = mxGetPr(mxGetField(aux, 0, "emitted_W_m2"));

    status = ngfSingleEpoch(&in, &out);
    if (status != NGF_OK) {
        mxDestroyArray(aux);
        mexErrMsgIdAndTxt(GW_ID "inputRange",
            "the C kernel refused this epoch: %s (code %d) -> check the epoch's ranges; MATLAB validates them before the call, so reaching this means one slipped through",
            ngfStatusText(status), status);
    }

    plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    memcpy(mxGetPr(plhs[2]), Fk, 3 * n * 3 * sizeof(double));
    plhs[3] = aux;
}
