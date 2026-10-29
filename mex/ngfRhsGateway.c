/* NGFRHSGATEWAY  MATLAB gateway to the native non-gravitational right-hand side (LP-C, WP-C9).
 *
 * Built by NonGravForces.buildMex into +NonGravForces/+app/private/rhsGateway.<mexext>,
 * visible to NonGravForces.app.DynamicsSource and to nothing else.
 *
 * CALLS
 *   [handle, entry] = rhsGateway('create', setup)   setup: the struct DynamicsSource.prepare builds
 *   A    = rhsGateway('accel', handle, t, R, V)        R, V, A: 3 x m
 *   diag = rhsGateway('diagnose', handle, t, r, v)     37 x 1 (see ngfRhs.h)
 *   rhsGateway('base', handle, baseEntry)              uint64 address of a gj_native_entry, or 0
 *   rhsGateway('destroy', handle)
 *   abi  = rhsGateway('abi')                           GJ_NATIVE_ABI_VERSION this was built against
 *   n    = rhsGateway('live')                          contexts alive in this process
 *
 * HANDLES ARE ADDRESSES AND ARE CHECKED AGAINST A REGISTRY before use, so a stale or
 * invented handle is refused by name instead of dereferenced. While any context is alive
 * the MEX is locked (mexLock), because gjToolbox's stepper holds a pointer into it and
 * "clear mex" must not unload the code under that pointer.
 *
 * Created 2026-09-16 by Claude Opus 5 (Anthropic).
 */
#include "mex.h"
#include "ngfRhs.h"
#include <stdint.h>
#include "gj_native.h"

#include <string.h>

#define GW_ID "NonGravForces:rhsGateway:"
#define NGF_LIVE_MAX 64

static ngfRhsContext *live[NGF_LIVE_MAX];

static int registryFind(ngfRhsContext *c)
{
    int i;
    for (i = 0; i < NGF_LIVE_MAX; i++) {
        if (c != NULL && live[i] == c) {
            return i;
        }
    }
    return -1;
}

static ngfRhsContext *handleArg(const mxArray *a)
{
    ngfRhsContext *c;
    if (!mxIsUint64(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputHandle", "the handle must be a uint64 scalar -> use the one 'create' returned");
    }
    c = (ngfRhsContext *)(uintptr_t)(*(uint64_T *)mxGetData(a));
    if (registryFind(c) < 0) {
        mexErrMsgIdAndTxt(GW_ID "inputHandle", "the handle names no live context (destroyed, or from another process) -> prepare the source again");
    }
    return c;
}

static const mxArray *field(const mxArray *s, const char *name, int required)
{
    const mxArray *f = mxGetField(s, 0, name);
    if (f == NULL && required) {
        mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup has no field %s -> build it with NonGravForces.app.DynamicsSource.prepare", name);
    }
    return f;
}

static double scalar(const mxArray *s, const char *name)
{
    const mxArray *f = field(s, name, 1);
    if (!mxIsDouble(f) || mxGetNumberOfElements(f) != 1) {
        mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup.%s must be a double scalar", name);
    }
    return mxGetScalar(f);
}

/* a double array of exactly count elements; count < 0 accepts any size; an empty field gives NULL */
static const double *array(const mxArray *s, const char *name, long count)
{
    const mxArray *f = field(s, name, 1);
    if (mxIsEmpty(f)) {
        return NULL;
    }
    if (!mxIsDouble(f) || mxIsComplex(f)) {
        mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup.%s must be a real double array", name);
    }
    if (count >= 0 && (long)mxGetNumberOfElements(f) != count) {
        mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup.%s has %d elements, %ld expected", name, (int)mxGetNumberOfElements(f), count);
    }
    return mxGetPr(f);
}

static void text(const mxArray *s, const char *name, char *buf, int size)
{
    const mxArray *f = field(s, name, 1);
    if (!mxIsChar(f) || mxGetString(f, buf, (mwSize)size) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup.%s must be a char row shorter than %d", name, size);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char cmd[32];
    if (nrhs < 1 || !mxIsChar(prhs[0]) || mxGetString(prhs[0], cmd, sizeof cmd) != 0) {
        mexErrMsgIdAndTxt(GW_ID "inputCommand", "the first argument must be a command -> create, accel, diagnose, base, destroy, abi, live");
    }
    if (strcmp(cmd, "abi") == 0) {
        plhs[0] = mxCreateDoubleScalar((double)GJ_NATIVE_ABI_VERSION);
        return;
    }
    if (strcmp(cmd, "live") == 0) {
        int i, n = 0;
        for (i = 0; i < NGF_LIVE_MAX; i++) {
            n += (live[i] != NULL);
        }
        plhs[0] = mxCreateDoubleScalar((double)n);
        return;
    }
    if (strcmp(cmd, "create") == 0) {
        ngfRhsSetup S;
        ngfRhsContext *c = NULL;
        const mxArray *st;
        const double *sm;
        char msis21[1024], hwm14[1024];
        int status, slot, n, nT, nb, cells, i;
        if (nrhs != 2 || !mxIsStruct(prhs[1])) {
            mexErrMsgIdAndTxt(GW_ID "inputSetup", "create takes one setup struct");
        }
        st = prhs[1];
        memset(&S, 0, sizeof S);
        S.cLight = scalar(st, "cLight"); S.sigmaSB = scalar(st, "sigmaSB"); S.kB = scalar(st, "kB");
        S.au = scalar(st, "au"); S.Rsun = scalar(st, "Rsun"); S.Rmoon = scalar(st, "Rmoon");
        S.aEarth = scalar(st, "aEarth"); S.bEarth = scalar(st, "bEarth"); S.fEarth = scalar(st, "fEarth");
        S.e2Earth = scalar(st, "e2Earth"); S.tsi = scalar(st, "tsi"); S.amu = scalar(st, "amu");
        /* THE GATEWAY COPIES BY NAME, so a constant added to the C struct and to the MATLAB packing
           is still absent here until this line names it. eV and torr arrived as uninitialised memory
           on 2026-09-17 and sesam's Langmuir constant became garbage; the seventh parity
           configuration caught it at 6.5e-02 against a 1e-13 bound, which is what that configuration
           was added for one commit earlier. */
        S.eV = scalar(st, "eV"); S.torr = scalar(st, "torr");
        sm = array(st, "speciesMass", NGF_SPECIES);
        for (i = 0; i < NGF_SPECIES; i++) { S.speciesMass[i] = sm[i]; }
        n = (int)scalar(st, "n"); S.n = n;
        S.area = array(st, "area", n); S.normal = array(st, "normal", 3L * n); S.centroid = array(st, "centroid", 3L * n);
        S.optVis = array(st, "optVis", 3L * n); S.optIR = array(st, "optIR", 3L * n);
        S.internalShare = array(st, "internalShare", -1);
        sm = array(st, "rcm", 3); S.rcm[0] = sm[0]; S.rcm[1] = sm[1]; S.rcm[2] = sm[2];
        S.gsi = (int)scalar(st, "gsi"); S.cd = scalar(st, "cd"); S.specular = scalar(st, "specular"); S.twK = scalar(st, "twK");
        S.cdScale = scalar(st, "cdScale"); S.crScale = scalar(st, "crScale");
        S.accommodation = (int)scalar(st, "accommodation"); S.alphaFixed = scalar(st, "alphaFixed"); S.surfaceMassAmu = scalar(st, "surfaceMassAmu");
        S.density = (int)scalar(st, "density"); S.expRho0 = scalar(st, "expRho0"); S.expH0 = scalar(st, "expH0"); S.expH = scalar(st, "expH");
        S.wind = (int)scalar(st, "wind"); S.eclipse = (int)scalar(st, "eclipse"); S.lunar = (int)scalar(st, "lunar");
        S.erp = (int)scalar(st, "erp"); S.nRings = (int)scalar(st, "nRings"); S.nSectors = (int)scalar(st, "nSectors");
        S.thermal = (int)scalar(st, "thermal"); S.internalPower = scalar(st, "internalPower");
        S.shadow = (int)scalar(st, "shadow"); S.lutResolution = (int)scalar(st, "lutResolution");
        cells = 6 * (S.lutResolution + 1) * (S.lutResolution + 1);
        S.lut = array(st, "lut", (S.shadow == NGF_SHADOW_LUT) ? (long)n * cells : -1);
        S.nChi = (int)scalar(st, "nChi"); S.nArc = (int)scalar(st, "nArc");
        S.xChi = array(st, "xChi", -1); S.wChi = array(st, "wChi", -1); S.xArc = array(st, "xArc", -1); S.wArc = array(st, "wArc", -1);
        S.flowTable = (int)scalar(st, "flowTable");
        S.ftN = (int)scalar(st, "ftN"); S.ftNs = (int)scalar(st, "ftNs"); S.ftNShaded = (int)scalar(st, "ftNShaded");
        S.ftHalfWidthTan = scalar(st, "ftHalfWidthTan"); S.ftT0 = scalar(st, "ftT0"); S.ftH = scalar(st, "ftH");
        S.ftS0 = scalar(st, "ftS0"); S.ftHs = scalar(st, "ftHs");
        S.ftShaded = array(st, "ftShaded", S.flowTable ? (long)S.ftNShaded : -1);
        S.ftPlus = array(st, "ftPlus", S.flowTable ? 24L * S.ftNShaded * S.ftN * S.ftN * S.ftNs : -1);
        S.ftMinus = array(st, "ftMinus", S.flowTable ? 24L * S.ftNShaded * S.ftN * S.ftN * S.ftNs : -1);
        S.t0 = scalar(st, "t0"); S.dt = scalar(st, "dt"); nT = (int)scalar(st, "nT"); S.nT = nT;
        S.R = array(st, "R", 9L * nT); S.Rdot = array(st, "Rdot", 9L * nT); S.sun = array(st, "sun", 3L * nT);
        S.moon = array(st, "moon", -1); S.eclipticNormal = array(st, "eclipticNormal", -1); S.mass = array(st, "mass", nT);
        S.attitude = (int)scalar(st, "attitude"); S.alongVelocity = (int)scalar(st, "alongVelocity"); S.senseTranspose = (int)scalar(st, "senseTranspose");
        S.nq = (int)scalar(st, "nq"); S.tq = array(st, "tq", -1); S.q = array(st, "q", -1);
        nb = (int)scalar(st, "nb"); S.nb = nb; S.tb0 = scalar(st, "tb0");
        S.alignedEdges = array(st, "alignedEdges", nb + 1L);
        S.f107 = array(st, "f107", nb); S.f107a = array(st, "f107a", nb); S.ap7 = array(st, "ap7", 7L * nb); S.ap3 = array(st, "ap3", nb);
        S.doy = array(st, "doy", nb); S.yy = array(st, "yy", nb); S.utcSecAtBin = array(st, "utcSecAtBin", nb); S.utcDay = array(st, "utcDay", nb);
        S.nLat = (int)scalar(st, "nLat"); S.nLon = (int)scalar(st, "nLon"); S.nMonth = (int)scalar(st, "nMonth");
        S.latC = array(st, "latC", -1); S.lonC = array(st, "lonC", -1); S.tMonth = array(st, "tMonth", -1);
        S.rsut = array(st, "rsut", -1); S.rlut = array(st, "rlut", -1); S.rsdt = array(st, "rsdt", -1);
        S.ceresHeight = scalar(st, "ceresHeight");
        S.nScaleDay = (int)scalar(st, "nScaleDay"); S.scaleDay0 = scalar(st, "scaleDay0");
        S.scale = array(st, "scale", S.nScaleDay > 0 ? 6L * S.nScaleDay : -1);
        if (S.nScaleDay > 0 && S.scale == NULL) {
            mexErrMsgIdAndTxt(GW_ID "inputSetup", "setup.scale is empty while nScaleDay is %d", S.nScaleDay);
        }
        text(st, "msis21Dir", msis21, sizeof msis21); S.msis21Dir = msis21;
        text(st, "hwm14Dir", hwm14, sizeof hwm14); S.hwm14Dir = hwm14;
        S.maxColumns = (int)scalar(st, "maxColumns");
        for (slot = 0; slot < NGF_LIVE_MAX && live[slot] != NULL; slot++) {
        }
        if (slot == NGF_LIVE_MAX) {
            mexErrMsgIdAndTxt(GW_ID "inputSetup", "%d contexts are alive -> release sources before preparing more", NGF_LIVE_MAX);
        }
        status = ngfRhsCreate(&S, &c);
        if (status != NGF_OK) {
            mexErrMsgIdAndTxt(GW_ID "inputSetup", "the native context could not be created (status %d: %s) -> check the setup", status, ngfStatusText(status));
        }
        live[slot] = c;
        mexLock();
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
        *(uint64_T *)mxGetData(plhs[0]) = (uint64_T)(uintptr_t)c;
        if (nlhs > 1) {
            plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
            *(uint64_T *)mxGetData(plhs[1]) = (uint64_T)(uintptr_t)ngfRhsEntry(c);
        }
        return;
    }
    if (strcmp(cmd, "destroy") == 0) {
        ngfRhsContext *c;
        int slot;
        if (nrhs != 2) {
            mexErrMsgIdAndTxt(GW_ID "inputCommand", "destroy takes the handle");
        }
        c = handleArg(prhs[1]);
        slot = registryFind(c);
        live[slot] = NULL;
        ngfRhsDestroy(c);
        mexUnlock();
        return;
    }
    if (strcmp(cmd, "base") == 0) {
        ngfRhsContext *c;
        uint64_T addr;
        if (nrhs != 3 || !mxIsUint64(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
            mexErrMsgIdAndTxt(GW_ID "inputCommand", "base takes the handle and a uint64 entry address (0 for none)");
        }
        c = handleArg(prhs[1]);
        addr = *(uint64_T *)mxGetData(prhs[2]);
        if (addr != 0 && gj_native_check((const gj_native_entry *)(uintptr_t)addr, 3, 0) != GJ_NATIVE_OK) {
            mexErrMsgIdAndTxt(GW_ID "inputEntry", "the base entry fails gj_native_check for a 3-state model -> publish a conforming entry");
        }
        ngfRhsSetBase(c, (void *)(uintptr_t)addr);
        return;
    }
    if (strcmp(cmd, "accel") == 0 || strcmp(cmd, "diagnose") == 0) {
        ngfRhsContext *c;
        double t;
        mwSize m;
        int status, diag = (strcmp(cmd, "diagnose") == 0);
        if (nrhs != 5) {
            mexErrMsgIdAndTxt(GW_ID "inputCommand", "%s takes the handle, t, R and V", cmd);
        }
        c = handleArg(prhs[1]);
        if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) ||
            mxGetM(prhs[3]) != 3 || mxGetM(prhs[4]) != 3 || mxGetN(prhs[3]) != mxGetN(prhs[4]) || mxGetN(prhs[3]) < 1) {
            mexErrMsgIdAndTxt(GW_ID "inputShape", "t must be a double scalar and R, V 3 x m double arrays of the same m");
        }
        t = mxGetScalar(prhs[2]);
        m = mxGetN(prhs[3]);
        if (diag) {
            if (m != 1) {
                mexErrMsgIdAndTxt(GW_ID "inputShape", "diagnose takes one state (3 x 1)");
            }
            plhs[0] = mxCreateDoubleMatrix(NGF_RHS_DIAG, 1, mxREAL);
            status = ngfRhsDiagnose(c, t, mxGetPr(prhs[3]), mxGetPr(prhs[4]), mxGetPr(plhs[0]));
        } else {
            plhs[0] = mxCreateDoubleMatrix(3, m, mxREAL);
            status = ngfRhsAccel(c, t, mxGetPr(prhs[3]), mxGetPr(prhs[4]), mxGetPr(plhs[0]), (int)m);
        }
        if (status != NGF_OK) {
            mexErrMsgIdAndTxt((status == NGF_DOMAIN) ? GW_ID "dataRange" : GW_ID "numericModel",
                "the native evaluation at t = %.6f s refused (status %d): %s", t, status, ngfRhsLastError(c));
        }
        return;
    }
    mexErrMsgIdAndTxt(GW_ID "inputCommand", "unknown command %s -> create, accel, diagnose, base, destroy, abi, live", cmd);
}
