/*
 * ngfRhs.c
 * --------
 * Implementation of the native non-gravitational right-hand side (LP-C, WP-C9).
 * See ngfRhs.h for what it is. Each block names the MATLAB unit it follows and,
 * where the last digits depend on it, keeps that unit's association order.
 *
 * Created 2026-09-16 by Claude Opus 5 (Anthropic).
 */
#include "ngfRhs.h"
#include "nrlmsise-00.h"
#include "gj_native.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NGF_PI 3.14159265358979323846
#define NGF_REALMIN 2.2250738585072014e-308
#define NGF_BIN_S 10800.0
#define NGF_OWN_MAX 64

/* the Fortran shims and the HWM14 loader */
extern void ngf_msis21_init(const char *dir, int n, int aphistory, int anomalous, int *status);
extern void ngf_msis21_calc(double day, double utsec, double z, double lat, double lon,
                            double f107a, double f107, const double *ap, double *tn, double *dn, double *tex);
extern void ngf_hwm14_calc(int iyd, double sec, double alt, double glat, double glon, double ap2, double *w);
extern int ngfHwm14Load(const char *dir, char *missing, int missingSize);

struct ngfRhsContext {
    ngfRhsSetup s;                   /* array members re-pointed at the context's own copies */
    ngfMacro macro;
    ngfConstants con;
    double *own[NGF_OWN_MAX];
    int nOwn;
    int *lutAllOnes;
    int *band;
    int mTiles, mMax, ringMax;
    double *srcDir, *srcP, *fRad, *tileN, *tileA, *tileAlb, *tileLw;
    double *fN, *fP, *fT, *ring, *beam, *tPanel, *Fk, *baseA;
    const gj_native_entry *base;
    gj_native_entry entry;
    char msis21Dir[1024], hwm14Dir[1024];
    int msis21Ready, hwm14Ready;
    char lastError[256];
};

/* ========================================================================= utilities */
static int ngfOwnCopy(ngfRhsContext *c, const double **field, size_t count, int required)
{
    double *y;
    if (*field == NULL || count == 0) {
        *field = NULL;
        return required ? NGF_NULL_ARGUMENT : NGF_OK;
    }
    if (c->nOwn >= NGF_OWN_MAX) {
        return NGF_ALLOC;
    }
    y = (double *)malloc(count * sizeof(double));
    if (y == NULL) {
        return NGF_ALLOC;
    }
    memcpy(y, *field, count * sizeof(double));
    c->own[c->nOwn++] = y;
    *field = y;
    return NGF_OK;
}

static double *ngfOwnZeros(ngfRhsContext *c, size_t count)
{
    double *y;
    if (c->nOwn >= NGF_OWN_MAX) {
        return NULL;
    }
    y = (double *)calloc(count > 0 ? count : 1, sizeof(double));
    if (y != NULL) {
        c->own[c->nOwn++] = y;
    }
    return y;
}

static void cross3(const double *a, const double *b, double *o)
{
    o[0] = a[1] * b[2] - a[2] * b[1];
    o[1] = a[2] * b[0] - a[0] * b[2];
    o[2] = a[0] * b[1] - a[1] * b[0];
}

static double dot3(const double *a, const double *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double norm3(const double *a)
{
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

static void mul33(const double *M, const double *x, double *y)      /* y = M x, M column-major */
{
    y[0] = M[0] * x[0] + M[3] * x[1] + M[6] * x[2];
    y[1] = M[1] * x[0] + M[4] * x[1] + M[7] * x[2];
    y[2] = M[2] * x[0] + M[5] * x[1] + M[8] * x[2];
}

static void mul33t(const double *M, const double *x, double *y)     /* y = M' x */
{
    y[0] = M[0] * x[0] + M[1] * x[1] + M[2] * x[2];
    y[1] = M[3] * x[0] + M[4] * x[1] + M[5] * x[2];
    y[2] = M[6] * x[0] + M[7] * x[1] + M[8] * x[2];
}

static double deg(double rad)
{
    return rad * (180.0 / NGF_PI);                   /* MATLAB rad2deg */
}

/* erfcx(y) = exp(y^2) erfc(y), y >= 0, without overflow (MATLAB erfcx) */
static double ngfErfcx(double y)
{
    double iy2, term, sum;
    int k;
    if (y < 26.0) {
        return exp(y * y) * erfc(y);
    }
    iy2 = 1.0 / (2.0 * y * y);
    term = 1.0;
    sum = 1.0;
    for (k = 1; k <= 8; k++) {
        term = -term * (2.0 * k - 1.0) * iy2;
        sum = sum + term;
    }
    return sum / (y * sqrt(NGF_PI));
}

static int ngfFail(ngfRhsContext *c, int status, const char *what)
{
    snprintf(c->lastError, sizeof c->lastError, "%s", what);
    return status;
}

/* ========================================================================= tables */
static int ngfTableIndex(const ngfRhsContext *c, double t, int *i, double *s)
{
    double x = (t - c->s.t0) / c->s.dt;
    int k;
    if (!(x >= 0.0) || !(x <= (double)(c->s.nT - 1))) {
        return NGF_DOMAIN;
    }
    k = (int)floor(x);
    if (k > c->s.nT - 2) {
        k = c->s.nT - 2;
    }
    *i = k;
    *s = x - (double)k;
    return NGF_OK;
}

/* cubic Lagrange through four table nodes about node i (clamped at the ends); exact at the nodes */
static void ngfLagrange(const ngfRhsContext *c, const double *tab, int width, int i, double s, double *out)
{
    int i0 = i - 1, j, k;
    double node[4], w[4];
    if (i0 < 0) {
        i0 = 0;
    }
    if (i0 > c->s.nT - 4) {
        i0 = c->s.nT - 4;
    }
    for (j = 0; j < 4; j++) {
        node[j] = (double)(i0 + j - i);
    }
    for (j = 0; j < 4; j++) {
        w[j] = 1.0;
        for (k = 0; k < 4; k++) {
            if (k != j) {
                w[j] = w[j] * ((s - node[k]) / (node[j] - node[k]));
            }
        }
    }
    for (k = 0; k < width; k++) {
        out[k] = w[0] * tab[width * i0 + k] + w[1] * tab[width * (i0 + 1) + k] +
                 w[2] * tab[width * (i0 + 2) + k] + w[3] * tab[width * (i0 + 3) + k];
    }
}

/* cubic Hermite of R with its rate Rdot: exact at the nodes */
static void ngfHermiteR(const ngfRhsContext *c, int i, double s, double *R)
{
    const double *R0 = c->s.R + 9 * i, *R1 = c->s.R + 9 * (i + 1);
    const double *D0 = c->s.Rdot + 9 * i, *D1 = c->s.Rdot + 9 * (i + 1);
    double s2 = s * s, s3 = s2 * s, dt = c->s.dt;
    double h00 = 2.0 * s3 - 3.0 * s2 + 1.0, h10 = s3 - 2.0 * s2 + s, h01 = -2.0 * s3 + 3.0 * s2, h11 = s3 - s2;
    int k;
    for (k = 0; k < 9; k++) {
        R[k] = h00 * R0[k] + (h10 * dt) * D0[k] + h01 * R1[k] + (h11 * dt) * D1[k];
    }
}

/* the 3-hour UTC bin: true edges for the calendar, the node-aligned edges for the drivers */
static int ngfBin(const ngfRhsContext *c, double t, int aligned)
{
    int lo, hi, mid;
    if (!aligned) {
        double x = floor((t - c->s.tb0) / NGF_BIN_S);
        if (!(x >= 0.0) || !(x < (double)c->s.nb)) {
            return -1;
        }
        return (int)x;
    }
    if (!(t >= c->s.alignedEdges[0]) || !(t <= c->s.alignedEdges[c->s.nb])) {
        return -1;
    }
    lo = 0;
    hi = c->s.nb;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (t >= c->s.alignedEdges[mid]) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

/* ========================================================================= attitude */
/* referenceSystem.interpolateQuaternions ("slerp", nodes normalised and sign-continued in prepare), then
 * referenceSystem.rotationFromQuaternion: BODY -> GCRS, column-major */
static int ngfAttitudeProvider(const ngfRhsContext *c, double t, double *R)
{
    const double *tq = c->s.tq, *q = c->s.q;
    int nq = c->s.nq, lo = 0, hi = nq - 1, mid, k;
    double frac, cosH, half, wFrom, wTo, qi[4], nrm, q0, q1, q2, q3, tmp;
    if (!(t >= tq[0]) || !(t <= tq[nq - 1])) {
        return NGF_DOMAIN;
    }
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (t >= tq[mid]) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    frac = (t - tq[lo]) / (tq[lo + 1] - tq[lo]);
    cosH = q[4 * lo] * q[4 * lo + 4] + q[4 * lo + 1] * q[4 * lo + 5] + q[4 * lo + 2] * q[4 * lo + 6] + q[4 * lo + 3] * q[4 * lo + 7];
    cosH = (cosH < -1.0) ? -1.0 : ((cosH > 1.0) ? 1.0 : cosH);
    half = acos(cosH);
    if (half < 1e-6) {
        wFrom = 1.0 - frac;
        wTo = frac;
    } else {
        wFrom = sin((1.0 - frac) * half) / sin(half);
        wTo = sin(frac * half) / sin(half);
    }
    for (k = 0; k < 4; k++) {
        qi[k] = wFrom * q[4 * lo + k] + wTo * q[4 * lo + 4 + k];
    }
    nrm = sqrt(qi[0] * qi[0] + qi[1] * qi[1] + qi[2] * qi[2] + qi[3] * qi[3]);
    q0 = qi[0] / nrm; q1 = qi[1] / nrm; q2 = qi[2] / nrm; q3 = qi[3] / nrm;
    R[0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    R[3] = 2.0 * (q1 * q2 + q0 * q3);
    R[6] = 2.0 * (q1 * q3 - q0 * q2);
    R[1] = 2.0 * (q1 * q2 - q0 * q3);
    R[4] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    R[7] = 2.0 * (q2 * q3 + q0 * q1);
    R[2] = 2.0 * (q1 * q3 + q0 * q2);
    R[5] = 2.0 * (q2 * q3 - q0 * q1);
    R[8] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    if (c->s.senseTranspose) {
        tmp = R[1]; R[1] = R[3]; R[3] = tmp;
        tmp = R[2]; R[2] = R[6]; R[6] = tmp;
        tmp = R[5]; R[5] = R[7]; R[7] = tmp;
    }
    return NGF_OK;
}

/* geom.AttitudeProvider.nominal: columns are the body axes in GCRS (+Z nadir, +X along or against v) */
static int ngfAttitudeNominal(int along, const double *r, const double *v, double *R)
{
    double h[3], rh[3], hh[3], a[3], xb[3], yb[3], zb[3], nr = norm3(r), nh;
    int k;
    cross3(r, v, h);
    nh = norm3(h);
    if (!(nh >= 1e-6 * nr * norm3(v))) {
        return NGF_DOMAIN;
    }
    for (k = 0; k < 3; k++) {
        zb[k] = -r[k] / nr;
        rh[k] = r[k] / nr;
        hh[k] = h[k] / nh;
    }
    cross3(hh, rh, a);
    for (k = 0; k < 3; k++) {
        xb[k] = along ? a[k] : -a[k];
    }
    cross3(zb, xb, yb);
    for (k = 0; k < 3; k++) {
        R[k] = xb[k];
        R[3 + k] = yb[k];
        R[6 + k] = zb[k];
    }
    return NGF_OK;
}

/* ========================================================================= geodetic */
/* referenceSystem.cartesianToGeodetic (Vermeille), WGS84: out = [lat rad; lon rad; h m] */
static int ngfGeodetic(const ngfRhsContext *c, const double *x, double *out)
{
    double a = c->s.aEarth, e2 = c->s.e2Earth, eq, p, q, r, s, t, u, v, w, k, D, hyp;
    eq = x[0] * x[0] + x[1] * x[1];
    p = eq / (a * a);
    q = (1.0 - e2) / (a * a) * (x[2] * x[2]);
    r = (p + q - e2 * e2) / 6.0;
    s = e2 * e2 * p * q / (4.0 * (r * r * r));
    if (!(r > 0.0) || !(s * (2.0 + s) >= 0.0)) {
        return NGF_DOMAIN;
    }
    t = cbrt(1.0 + s + sqrt(s * (2.0 + s)));
    u = r * (1.0 + t + 1.0 / t);
    v = sqrt(u * u + e2 * e2 * q);
    w = e2 * (u + v - q) / (2.0 * v);
    k = sqrt(u + v + w * w) - w;
    D = k * sqrt(eq) / (k + e2);
    hyp = sqrt(D * D + x[2] * x[2]);
    out[0] = 2.0 * atan2(x[2], D + hyp);
    out[1] = atan2(x[1], x[0]);
    out[2] = (k + e2 - 1.0) / k * hyp;
    return NGF_OK;
}

/* ========================================================================= eclipse */
/* env.eclipseFactor localBody */
static double ngfBodyShadow(int model, const double *r, const double *rs, const double *rb, double Rb, double Rsun)
{
    double d[3], sv[3], nd, ns, cx[3], nu = 1.0;
    int k;
    for (k = 0; k < 3; k++) {
        d[k] = rs[k] - r[k];
        sv[k] = rb[k] - r[k];
    }
    nd = norm3(d);
    ns = norm3(sv);
    if (!(ns < nd)) {
        return 1.0;
    }
    if (model == NGF_ECL_CYLINDRICAL) {
        double u[3], along, pv[3];
        for (k = 0; k < 3; k++) {
            u[k] = d[k] / nd;
        }
        along = dot3(sv, u);
        for (k = 0; k < 3; k++) {
            pv[k] = sv[k] - along * u[k];
        }
        return (along > 0.0 && norm3(pv) < Rb) ? 0.0 : 1.0;
    }
    {
        double aa = asin(fmin(1.0, Rsun / nd)), bb = asin(fmin(1.0, Rb / ns)), cc;
        cross3(d, sv, cx);
        cc = atan2(norm3(cx), dot3(d, sv));
        if (cc >= aa + bb) {
            nu = 1.0;
        } else if (cc <= bb - aa) {
            nu = 0.0;
        } else if (cc <= aa - bb) {
            nu = 1.0 - (bb / aa) * (bb / aa);
        } else {
            double xx = (cc * cc + aa * aa - bb * bb) / (2.0 * cc);
            double yy = sqrt(fmax(0.0, aa * aa - xx * xx));
            double A = aa * aa * atan2(yy, xx) + bb * bb * atan2(yy, cc - xx) - cc * yy;
            nu = 1.0 - A / (NGF_PI * (aa * aa));
        }
        return fmin(1.0, fmax(0.0, nu));
    }
}

/* env.eclipseFactor localSolaars (Robertson 2015, eqs 3.80-3.86, Table 3.2) */
static double ngfSolaars(const double *r, const double *rs, const double *e, double aE, double bE)
{
    static const double T[8][5] = {
        {1,  0.1715,    -0.1423,     0.01061,  -0.01443},
        {2,  0.008162,   0.3401,     0,         0},
        {1,  260.9,     -0.4661,     27.81,    -0.009437},
        {3, -0.006119,   1.176,      6.385,     0},
        {1,  87.56,     -0.09188,    19.30,    -0.01089},
        {2,  0.002047,   6.409,      0,         0},
        {1,  61.98,     -0.1629,     27.87,    -0.02217},
        {1,  6.413,     -0.0002593, -0.01479,  -0.1318}};
    double nrs = norm3(rs), Rh[3], rR, rEv[3], rEperp, rEpar2, sO, rE, rRc, a[8], num, nu;
    int k;
    for (k = 0; k < 3; k++) {
        Rh[k] = rs[k] / nrs;
    }
    rR = -dot3(r, Rh);
    for (k = 0; k < 3; k++) {
        rEv[k] = r[k] + rR * Rh[k];
    }
    rEperp = dot3(e, rEv);
    rEpar2 = fmax(0.0, (rEv[0] * rEv[0] + rEv[1] * rEv[1] + rEv[2] * rEv[2]) - rEperp * rEperp);
    sO = aE / bE;
    rE = sqrt(rEperp * rEperp + sO * sO * rEpar2) / 1e6;
    rRc = fmin(45.0, fmax(2.0, rR / 1e6));
    for (k = 0; k < 8; k++) {
        const double *b = T[k] + 1;
        if (T[k][0] == 1.0) {
            a[k] = b[0] * exp(b[1] * rRc) + b[2] * exp(b[3] * rRc);
        } else if (T[k][0] == 2.0) {
            a[k] = b[0] * rRc + b[1];
        } else {
            a[k] = b[0] * pow(rRc, b[1]) + b[2];
        }
    }
    num = 1.0 + a[0] + a[1] + a[0] * tanh(a[2] * (rE - a[3])) + a[1] * tanh(a[4] * (rE - a[5])) + tanh(a[6] * (rE - a[7]));
    nu = num / (2.0 + 2.0 * a[0] + 2.0 * a[1]);
    if (rR <= 0.0) {
        nu = 1.0;
    }
    return fmin(1.0, fmax(0.0, nu));
}

/* ========================================================================= self-shadow table */
/* geom.selfShadow localLookup: plate k (0-based), one unit direction */
static double ngfLut(const ngfRhsContext *c, int k, const double *d)
{
    static const int UA[3] = {1, 2, 0}, VA[3] = {2, 0, 1};
    const double *L = c->s.lut;
    int n = c->s.lutResolution, N = c->s.n, m = 0, face, i0, j0, c00, c10, c01, c11;
    double best = fabs(d[0]), dm, aa, bb, x, y, wx, wy, T00, T10, T01, T11, f;
    if (fabs(d[1]) > best) { m = 1; best = fabs(d[1]); }
    if (fabs(d[2]) > best) { m = 2; }
    dm = d[m];
    face = 2 * (m + 1) - (dm > 0.0);
    aa = atan(d[UA[m]] / fabs(dm));
    bb = atan(d[VA[m]] / fabs(dm));
    x = (aa + NGF_PI / 4.0) / (NGF_PI / 2.0) * n;
    y = (bb + NGF_PI / 4.0) / (NGF_PI / 2.0) * n;
    i0 = (int)fmin(fmax(floor(x), 0.0), (double)(n - 1));
    j0 = (int)fmin(fmax(floor(y), 0.0), (double)(n - 1));
    wx = fmin(fmax(x - i0, 0.0), 1.0);
    wy = fmin(fmax(y - j0, 0.0), 1.0);
    c00 = (face - 1) * (n + 1) * (n + 1) + j0 * (n + 1) + i0;
    c10 = c00 + 1;
    c01 = c00 + n + 1;
    c11 = c01 + 1;
    T00 = L[N * c00 + k]; T10 = L[N * c10 + k]; T01 = L[N * c01 + k]; T11 = L[N * c11 + k];
    f = T00 + wx * (T10 - T00) + wy * (T01 - T00) + (wx * wy) * ((T11 - T10) - (T01 - T00));
    return fmin(1.0, fmax(0.0, f));
}

/* ========================================================================= flow shadow */
/* geom.flowShadow (Method "table") for one epoch, the same arithmetic in the same order as its
   localTrilinear: returns 1 and fills fN, fP, fT when the epoch lies inside the table, 0 otherwise
   (the caller then computes it exactly, as geom.flowShadow does) */
static int ngfFlowTable(ngfRhsContext *c, const double *vB, double U, double T)
{
    const ngfRhsSetup *S = &c->s;
    int N = S->n, n = S->ftN, nSh = S->ftNShaded, P = 24 * S->ftNShaded, nn, j, p;
    long c000, c100, c010, c110, c001, c101, c011, c111;
    double v[3], ax, ty, tz, sO, sTop, fx, fy, fs, ixd, iyd, isd, wx, wy, ws;
    const double *V;
    for (j = 0; j < 3; j++) {
        v[j] = vB[j] / U;
    }
    ax = fabs(v[0]);
    if (!(ax > 0.0)) {
        return 0;
    }
    ty = v[1] / ax; tz = v[2] / ax;
    sO = U / sqrt(2.0 * S->kB * T / S->speciesMass[1]);
    sTop = S->ftS0 + (S->ftNs - 1) * S->ftHs;
    if (!(fabs(ty) <= S->ftHalfWidthTan && fabs(tz) <= S->ftHalfWidthTan && sO >= S->ftS0 && sO <= sTop)) {
        return 0;
    }
    V = (v[0] > 0.0) ? S->ftPlus : S->ftMinus;
    fx = (ty - S->ftT0) / S->ftH; ixd = fmin(fmax(floor(fx), 0.0), n - 2.0); wx = fmin(fmax(fx - ixd, 0.0), 1.0);
    fy = (tz - S->ftT0) / S->ftH; iyd = fmin(fmax(floor(fy), 0.0), n - 2.0); wy = fmin(fmax(fy - iyd, 0.0), 1.0);
    fs = (sO - S->ftS0) / S->ftHs; isd = fmin(fmax(floor(fs), 0.0), S->ftNs - 2.0); ws = fmin(fmax(fs - isd, 0.0), 1.0);
    nn = n * n;
    c000 = (long)ixd + (long)n * (long)iyd + (long)nn * (long)isd;
    c100 = c000 + 1; c010 = c000 + n; c110 = c010 + 1;
    c001 = c000 + nn; c101 = c100 + nn; c011 = c010 + nn; c111 = c110 + nn;
    for (j = 0; j < N * NGF_SPECIES; j++) {
        c->fN[j] = 1.0; c->fP[j] = 1.0; c->fT[j] = 1.0;
    }
    for (p = 0; p < P; p++) {
        double V000 = V[P * c000 + p], V001 = V[P * c001 + p], V010 = V[P * c010 + p], V011 = V[P * c011 + p];
        double a00 = V000 + wx * (V[P * c100 + p] - V000);
        double a10 = V010 + wx * (V[P * c110 + p] - V010);
        double a01 = V001 + wx * (V[P * c101 + p] - V001);
        double a11 = V011 + wx * (V[P * c111 + p] - V011);
        double b0 = a00 + wy * (a10 - a00);
        double b1 = a01 + wy * (a11 - a01);
        double f = b0 + ws * (b1 - b0);
        int sh = p % nSh, sp = (p / nSh) % NGF_SPECIES, r = p / (nSh * NGF_SPECIES);
        int k = (int)S->ftShaded[sh] - 1;
        if (r == 0) {
            c->fN[N * sp + k] = fmin(1.0, fmax(0.0, f));
        } else if (r == 1) {
            c->fP[N * sp + k] = fmin(1.0, fmax(0.0, f));
        } else {
            c->fT[N * sp + k] = f;
        }
    }
    return 1;
}

/* geom.flowShadow (Method "lut") for one epoch: fN, fP, fT laid out [N * j + k] for plate k, species j */
static void ngfFlowShadow(ngfRhsContext *c, const double *vB, double U, double T)
{
    const ngfRhsSetup *S = &c->s;
    int N = S->n, nG = S->nChi, nA = S->nArc, R = 6 * nG, k, j, g, seg, ia, ir;
    double v[3], e1[3], e2[3], a[3] = {0.0, 0.0, 0.0}, sj[NGF_SPECIES], av, ne;
    double *SN0 = c->ring, *SNO = c->ring + R, *SP0 = c->ring + 2 * R, *SPO = c->ring + 3 * R;
    double *ST0 = c->ring + 4 * R, *STO = c->ring + 5 * R, *chi = c->ring + 6 * R, *wchi = c->ring + 7 * R;
    for (j = 0; j < 3; j++) {
        v[j] = vB[j] / U;
    }
    for (j = 0; j < N * NGF_SPECIES; j++) {
        c->fN[j] = 1.0; c->fP[j] = 1.0; c->fT[j] = 1.0;
    }
    for (j = 0; j < NGF_SPECIES; j++) {
        sj[j] = U / sqrt(2.0 * S->kB * T / S->speciesMass[j]);
    }
    if (fabs(v[0]) > 0.9) { a[1] = 1.0; } else { a[0] = 1.0; }
    av = dot3(a, v);
    for (j = 0; j < 3; j++) {
        e1[j] = a[j] - av * v[j];
    }
    ne = norm3(e1);
    for (j = 0; j < 3; j++) {
        e1[j] = e1[j] / ne;
    }
    cross3(v, e1, e2);
    for (k = 0; k < N; k++) {
        const double *nk = S->normal + 3 * k;
        double nz, nx, ny, B, psi0, thetaN, chi1, chi2, vals[7], lo[6], hi[6];
        int perm[7], isH[7];
        if (c->lutAllOnes[k]) {
            continue;
        }
        nz = dot3(nk, v); nx = dot3(nk, e1); ny = dot3(nk, e2);
        B = hypot(nx, ny); psi0 = atan2(ny, nx); thetaN = atan2(B, nz);
        chi1 = fabs(NGF_PI / 2.0 - thetaN); chi2 = NGF_PI - chi1;
        vals[0] = 0.0; vals[1] = 3.0 / sj[4]; vals[2] = 3.0 / sj[1]; vals[3] = 3.0 / sj[5];
        vals[4] = chi1; vals[5] = chi2; vals[6] = NGF_PI;
        for (j = 0; j < 7; j++) {
            vals[j] = fmin(NGF_PI, fmax(0.0, vals[j]));
            perm[j] = j;
        }
        for (j = 1; j < 7; j++) {                              /* stable, as MATLAB's sort */
            double key = vals[j];
            int pk = perm[j], i2 = j - 1;
            while (i2 >= 0 && vals[i2] > key) {
                vals[i2 + 1] = vals[i2]; perm[i2 + 1] = perm[i2]; i2--;
            }
            vals[i2 + 1] = key; perm[i2 + 1] = pk;
        }
        for (j = 0; j < 7; j++) {
            isH[j] = (perm[j] == 4 || perm[j] == 5);
        }
        for (seg = 0; seg < 6; seg++) {
            lo[seg] = vals[seg];
            hi[seg] = vals[seg + 1];
        }
        for (g = 0; g < nG; g++) {                            /* ring ir = g * 6 + seg, MATLAB's reshape order */
            double tau = (S->xChi[g] + 1.0) / 2.0, wT = S->wChi[g] / 2.0;
            for (seg = 0; seg < 6; seg++) {
                int hl = isH[seg], hh = isH[seg + 1];
                double map = tau, dmap = 1.0;
                if (hl && !hh) {
                    map = tau * tau; dmap = 2.0 * tau;
                } else if (!hl && hh) {
                    map = 1.0 - (1.0 - tau) * (1.0 - tau); dmap = 2.0 * (1.0 - tau);
                } else if (hl && hh) {
                    map = 3.0 * (tau * tau) - 2.0 * (tau * tau * tau); dmap = 6.0 * tau * (1.0 - tau);
                }
                ir = g * 6 + seg;
                chi[ir] = lo[seg] + (hi[seg] - lo[seg]) * map;
                wchi[ir] = (((hi[seg] - lo[seg]) * dmap) * wT) * sin(chi[ir]);
            }
        }
        for (ir = 0; ir < R; ir++) {
            double cc = cos(chi[ir]), sc = sin(chi[ir]), A = nz * cc, bb = B * sc, delta;
            double ll = (B > 0.0) ? B : NGF_REALMIN;
            SN0[ir] = 0.0; SNO[ir] = 0.0; SP0[ir] = 0.0; SPO[ir] = 0.0; ST0[ir] = 0.0; STO[ir] = 0.0;
            delta = atan2(sqrt(fmax(bb * bb - A * A, 0.0)), -A);
            if (A >= bb) { delta = NGF_PI; }
            if (A <= -bb) { delta = 0.0; }
            for (ia = 0; ia < nA; ia++) {
                double u = delta * S->xArc[ia], wpsi = delta * S->wArc[ia], cu = cos(u);
                double cn = fmax(A + bb * cu, 0.0), psi = psi0 + u, cp = sc * cos(psi), sp = sc * sin(psi);
                double W[3], o, tang, wc, wco;
                for (j = 0; j < 3; j++) {
                    W[j] = e1[j] * cp + e2[j] * sp + v[j] * cc;
                }
                o = ngfLut(c, k, W);
                tang = cc * ll - ((nz / ll) * bb) * cu;
                wc = wpsi * cn;
                wco = (wpsi * o) * cn;
                SN0[ir] += wc;        SNO[ir] += wco;
                SP0[ir] += wc * cn;   SPO[ir] += wco * cn;
                ST0[ir] += wc * tang; STO[ir] += wco * tang;
            }
        }
        for (j = 0; j < NGF_SPECIES; j++) {
            double s = sj[j], es = exp(-(s * s)), IN0 = 0.0, INO = 0.0, IP0 = 0.0, IPO = 0.0, IT0 = 0.0, ITO = 0.0;
            for (ir = 0; ir < R; ir++) {
                double x = s * cos(chi[ir]), EG0, EJ3, EJ4, g3, g4;
                if (x >= 0.0) {
                    EG0 = (sqrt(NGF_PI) / 2.0) * exp(x * x - s * s) * erfc(-x);
                } else {
                    EG0 = (sqrt(NGF_PI) / 2.0) * es * ngfErfcx(-x);
                }
                EJ3 = es * (x * x + 1.0) / 2.0 + (x * x * x + 1.5 * x) * EG0;
                EJ4 = es * (0.5 * (x * x * x) + 1.25 * x) + (x * x * x * x + 3.0 * (x * x) + 0.75) * EG0;
                g3 = wchi[ir] * EJ3;
                g4 = wchi[ir] * EJ4;
                IN0 += g3 * SN0[ir]; INO += g3 * SNO[ir];
                IP0 += g4 * SP0[ir]; IPO += g4 * SPO[ir];
                IT0 += g4 * ST0[ir]; ITO += g4 * STO[ir];
            }
            c->fN[N * j + k] = (IN0 > 0.0) ? INO / IN0 : 1.0;
            c->fP[N * j + k] = (IP0 > 0.0) ? IPO / IP0 : 1.0;
            c->fT[N * j + k] = (!(B > 0.0) || !(fabs(IT0) > 0.0)) ? c->fN[N * j + k] : ITO / IT0;
        }
    }
}

/* ========================================================================= accommodation */
/* kernel.accommodation: alphaA (the adsorbate or the only state), coverage (-1: a single state), alphaS */
static void ngfAccommodation(const ngfRhsContext *c, const double *sp, double T, double vr, double *alphaA, double *coverage, double *alphaS)
{
    const ngfRhsSetup *S = &c->s;
    *coverage = -1.0;
    *alphaS = 0.0;
    if (S->accommodation == NGF_ACC_FIXED) {
        *alphaA = S->alphaFixed;
    } else if (S->accommodation == NGF_ACC_PILINSKI2010 ||
               S->accommodation == NGF_ACC_PILINSKI2011) {
        /* One model, two published fits: 7.50e-17 is Pilinski, Argrow & Palo (2010) eq. 8,
           4.98e-17 is Pilinski (2011) sect. 7.3. Same association order as the MATLAB. */
        double kL = (S->accommodation == NGF_ACC_PILINSKI2010) ? 7.50e-17 : 4.98e-17;
        double x = kL * sp[1] * T;
        *alphaA = x / (1.0 + x);
    } else {
        /* walker2014 and sesam share the two-state mixture and the mole-weighted mean mass, and
           differ in the pressure their Langmuir constant was fitted against and in Goodman's
           coefficient. Same association order as the MATLAB, which is what makes the two agree to
           rounding rather than to "about 1e-12". */
        static const int six[6] = {5, 0, 6, 2, 1, 3};          /* H He N N2 O O2, 0-based */
        double PO, KL, gc, msAmu, num = 0.0, den = 0.0, mean, mu;
        int i;
        if (S->accommodation == NGF_ACC_SESAM) {
            /* [P13] eq. (15): a RAM pressure, 140x the n_O kB T that walker2014 uses, with the
               free-molecular sphere coefficient whose re-emission term Appendix B removes. */
            double mO = S->speciesMass[1];
            double s2 = vr / sqrt(2.0 * S->kB * T / mO), ss = s2 * s2, s4 = ss * ss;
            double cdIn = (2.0 * ss + 1.0) / (sqrt(NGF_PI) * ss * s2) * exp(-ss)
                        + (4.0 * s4 + 4.0 * ss - 1.0) / (2.0 * s4) * erf(s2);
            double Eb = 5.7 * S->eV, kT = S->kB * T, Er = 0.5 * mO * vr * vr;
            /* eq. (26b), with the exponents combined analytically and erf(a)+erf(b) written as
               erfc(-a)-erfc(b): as printed it needs exp(Eb/kT) ~ 1e36 against exp(-(Eb+Er)/kT)
               ~ 1e-66, and the erf pair cancels to nothing once Er > Eb. */
            double ra = sqrt(Eb) - sqrt(Er), aa = ra / sqrt(kT), bb = sqrt(Er / kT);
            double g = sqrt(NGF_PI * kT * Er), s0;
            s0 = (g * (erfc(-aa) - erfc(bb)) + kT * (exp(-Er / kT) - exp(-ra * ra / kT)))
               / (g * (1.0 + erf(bb)) + kT * exp(-Er / kT));
            if (!(s0 > 0.0)) { s0 = 0.0; } else if (s0 > 1.0) { s0 = 1.0; }
            KL = s0 * (5e6 / S->torr) + 3e4 / S->torr;         /* eq. (22), Pa^-1 */
            PO = 0.5 * (sp[1] * mO) * vr * vr * cdIn;
            gc = 2.4;                                          /* eq. (19), SESAM's own */
            msAmu = 65.0;                                      /* Table C1, fixed and not an option */
        } else {
            PO = sp[1] * S->kB * T;
            KL = 1.44e6;
            gc = 3.0;                                          /* Paul et al. 2022 eq. 13 */
            msAmu = S->surfaceMassAmu;
        }
        *coverage = KL * PO / (1.0 + KL * PO);
        *alphaA = 1.0;
        for (i = 0; i < 6; i++) {
            num = num + sp[six[i]] * S->speciesMass[six[i]];
            den = den + sp[six[i]];
        }
        mean = num / ((den > NGF_REALMIN) ? den : NGF_REALMIN);
        mu = mean / (msAmu * S->amu);
        *alphaS = gc * mu / ((1.0 + mu) * (1.0 + mu));
    }
}

/* ========================================================================= gas force */
/* kernel.dragPanelForce with a flowShadow struct (flow != 0), a beam fraction (constantCd), or none */
static int ngfDrag(ngfRhsContext *c, double rho, const double *sp, double T, const double *vB, int flow, double *F)
{
    const ngfRhsSetup *S = &c->s;
    int N = S->n, k, j, i, status;
    double U = norm3(vB), v[3], alphaA = 0.0, coverage = -1.0, alphaS = 0.0;
    if (!(U > 0.0)) {
        return NGF_BAD_RANGE;
    }
    for (i = 0; i < 3; i++) {
        v[i] = vB[i] / U;
    }
    if (S->gsi != NGF_GSI_CONSTANT_CD) {
        ngfAccommodation(c, sp, T, U, &alphaA, &coverage, &alphaS);
    }
    for (k = 0; k < N; k++) {
        const double *nk = S->normal + 3 * k;
        double gk = dot3(nk, v), cr[3], lk, th, t[3], Cpm = 0.0, Ctm = 0.0, scale, fn, ft;
        cross3(nk, v, cr);
        lk = norm3(cr);
        for (i = 0; i < 3; i++) {
            t[i] = (-v[i] + gk * nk[i]);
        }
        if (lk > 0.0) {
            for (i = 0; i < 3; i++) { t[i] = t[i] / lk; }
        } else {
            t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;
        }
        th = atan2(lk, gk);
        if (S->gsi == NGF_GSI_CONSTANT_CD) {
            double inc, re, ctau;
            status = ngfGsiParts((ngfGsiModel)S->gsi, 0.0, th, 0.0, 0.0, 0.0, S->cd, S->specular, &inc, &re, &ctau);
            if (status != NGF_OK) { return status; }
            Cpm = rho * (inc + re);
            Ctm = rho * ctau;
        } else {
            for (j = 0; j < NGF_SPECIES; j++) {
                double rhoj = sp[j] * S->speciesMass[j], s, inc, re, ctj, CpRe;
                if (!(rhoj > 0.0)) {
                    continue;
                }
                s = U / sqrt(2.0 * S->kB * T / S->speciesMass[j]);
                status = ngfGsiParts((ngfGsiModel)S->gsi, s, th, alphaA, S->twK, T, S->cd, S->specular, &inc, &re, &ctj);
                if (status != NGF_OK) { return status; }
                CpRe = re;
                if (coverage >= 0.0) {
                    double incS, reS, ctS;
                    status = ngfGsiParts((ngfGsiModel)S->gsi, s, th, alphaS, S->twK, T, S->cd, S->specular, &incS, &reS, &ctS);
                    if (status != NGF_OK) { return status; }
                    CpRe = coverage * re + (1.0 - coverage) * reS;
                }
                if (flow) {
                    Cpm = Cpm + rhoj * (c->fP[N * j + k] * inc + c->fN[N * j + k] * CpRe);
                    Ctm = Ctm + rhoj * (c->fT[N * j + k] * ctj);
                } else {
                    Cpm = Cpm + rhoj * (inc + CpRe);
                    Ctm = Ctm + rhoj * ctj;
                }
            }
        }
        scale = 0.5 * (U * U) * S->area[k];
        if (S->gsi == NGF_GSI_CONSTANT_CD && S->shadow == NGF_SHADOW_LUT) {
            scale = scale * c->beam[k];
        }
        fn = -scale * Cpm;
        ft = scale * Ctm;
        for (i = 0; i < 3; i++) {
            F[3 * k + i] = nk[i] * fn + t[i] * ft;
        }
    }
    return NGF_OK;
}

/* ========================================================================= Earth tiles */
/* env.earthRadiance (knocke, ceres): tile normals (GCRS) and areas on the cap around r, and the
 * albedo and longwave exitance per tile; returns the tile count or a negative status */
static int ngfTiles(ngfRhsContext *c, const double *r, const double *Ritrs, double Searth, double utcDay, double utcSec)
{
    const ngfRhsSetup *S = &c->s;
    int nR = S->nRings, nS = S->nSectors, Mt = nR * nS, ri, si, ti, i;
    double Rt = S->aEarth + ((S->erp == NGF_ERP_CERES) ? S->ceresHeight : 0.0);
    double nr = norm3(r), rh[3], e1[3], e2[3], zx[3] = {0.0, 0.0, 1.0}, ne1, alphaMax;
    if (!(nr > Rt)) {
        return -NGF_DOMAIN;
    }
    for (i = 0; i < 3; i++) {
        rh[i] = r[i] / nr;
    }
    cross3(zx, rh, e1);
    ne1 = norm3(e1);
    if (ne1 < 1e-9) {
        e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0;
    } else {
        for (i = 0; i < 3; i++) { e1[i] = e1[i] / ne1; }
    }
    cross3(rh, e1, e2);
    alphaMax = acos(Rt / nr);
    for (ri = 0; ri < nR; ri++) {
        double aLo = ((double)ri / nR) * alphaMax, aHi = ((double)(ri + 1) / nR) * alphaMax;
        double aC = 0.5 * (aLo + aHi), area = Rt * Rt * (2.0 * NGF_PI / nS) * (cos(aLo) - cos(aHi));
        for (si = 0; si < nS; si++) {
            double ph = (2.0 * NGF_PI / nS) * ((double)(si + 1) - 0.5);
            ti = ri * nS + si;
            for (i = 0; i < 3; i++) {
                c->tileN[3 * ti + i] = cos(aC) * rh[i] + sin(aC) * (cos(ph) * e1[i] + sin(ph) * e2[i]);
            }
            c->tileA[ti] = area;
        }
    }
    if (S->erp == NGF_ERP_KNOCKE) {
        /* Knocke, Ries & Tapley (1988) Table 1, t0 = 1981-12-22 UTC (MJD 44960) */
        double days = (utcDay - 44960.0) + utcSec / 86400.0, seasonal = cos(2.0 * NGF_PI / 365.25 * days);
        for (ti = 0; ti < Mt; ti++) {
            double x = c->tileN[3 * ti + 2], P2 = 0.5 * (3.0 * (x * x) - 1.0);
            c->tileAlb[ti] = 0.34 + 0.10 * (x * seasonal) + 0.29 * P2;
            c->tileLw[ti] = (0.68 + (-0.07) * (x * seasonal) + (-0.18) * P2) * (Searth / 4.0);
        }
    } else {
        double tm = utcDay + utcSec / 86400.0, wt, f = S->fEarth;
        int k0, nLat = S->nLat, nLon = S->nLon, lo, hi, mid;
        if (!(tm >= S->tMonth[0]) || !(tm <= S->tMonth[S->nMonth - 1])) {
            return -NGF_DOMAIN;
        }
        lo = 0; hi = S->nMonth - 1;
        while (hi - lo > 1) {
            mid = (lo + hi) / 2;
            if (tm >= S->tMonth[mid]) { lo = mid; } else { hi = mid; }
        }
        k0 = lo;
        wt = (tm - S->tMonth[k0]) / (S->tMonth[k0 + 1] - S->tMonth[k0]);
        for (ti = 0; ti < Mt; ti++) {
            double nI[3], latD, lon, u, vv, wx, wy, val[3];
            int i0, j0, j1, side, fld;
            mul33t(Ritrs, c->tileN + 3 * ti, nI);
            latD = deg(atan2(nI[2], (1.0 - f) * (1.0 - f) * hypot(nI[0], nI[1])));
            lon = fmod(deg(atan2(nI[1], nI[0])), 360.0);
            if (lon < 0.0) { lon = lon + 360.0; }
            u = (fmin(fmax(latD, S->latC[0]), S->latC[nLat - 1]) - S->latC[0]) / (S->latC[1] - S->latC[0]);
            i0 = (int)fmin(floor(u), (double)(nLat - 2));
            wy = u - i0;
            vv = fmod(lon - S->lonC[0], 360.0);
            if (vv < 0.0) { vv = vv + 360.0; }
            vv = vv / (S->lonC[1] - S->lonC[0]);
            j0 = (int)floor(vv);
            wx = vv - j0;
            j0 = j0 % nLon;
            j1 = (j0 + 1) % nLon;
            for (fld = 0; fld < 3; fld++) {
                const double *G = (fld == 0) ? S->rsut : ((fld == 1) ? S->rlut : S->rsdt);
                double acc = 0.0;
                for (side = 0; side < 2; side++) {
                    int page = (k0 + side) * nLat * nLon;
                    double g00 = G[i0 + j0 * nLat + page], g10 = G[i0 + 1 + j0 * nLat + page];
                    double g01 = G[i0 + j1 * nLat + page], g11 = G[i0 + 1 + j1 * nLat + page];
                    double b = g00 + wy * (g10 - g00) + wx * (g01 - g00) + (wx * wy) * ((g11 - g10) - (g01 - g00));
                    acc = (side == 0) ? (1.0 - wt) * b : acc + wt * b;
                }
                val[fld] = acc;
            }
            c->tileAlb[ti] = (val[2] > 0.0) ? val[0] / val[2] : 0.0;
            c->tileLw[ti] = val[1];
        }
    }
    return Mt;
}

/* ========================================================================= one state */
typedef struct {
    double aD[3], aR[3], aT[3], rho, species[NGF_SPECIES], Tgas, vB[3], lla[3], nu, mass, Rb[9];
} ngfParts;

static int ngfEvalOne(ngfRhsContext *c, double t, const double *r, const double *v, double *aG, ngfParts *pt)
{
    const ngfRhsSetup *S = &c->s;
    int status, i, j, k, m, it = 0, bT, bA, M, N = S->n;
    double sT = 0.0, Ri[9], Rd[9], Rb[9], rI[3], vrel[3], vB[3], lla[3], sun[3], moon[3], ecl[3], mass;
    double rho = 0.0, sp[NGF_SPECIES], Tg = 0.0, utcSec, nu, FdSum[3], FrSum[3], FtSum[3], tau[3], aB[3];

    status = ngfTableIndex(c, t, &it, &sT);
    bT = ngfBin(c, t, 0);
    bA = ngfBin(c, t, 1);
    if (status != NGF_OK || bT < 0 || bA < 0) {
        return ngfFail(c, NGF_DOMAIN, "t lies outside the tables built at prepare (widen 'Margin_s')");
    }
    ngfHermiteR(c, it, sT, Ri);
    ngfLagrange(c, S->Rdot, 9, it, sT, Rd);
    ngfLagrange(c, S->sun, 3, it, sT, sun);
    if (S->moon != NULL) { ngfLagrange(c, S->moon, 3, it, sT, moon); }
    if (S->eclipticNormal != NULL) {
        double ne;
        ngfLagrange(c, S->eclipticNormal, 3, it, sT, ecl);
        ne = norm3(ecl);
        for (i = 0; i < 3; i++) { ecl[i] = ecl[i] / ne; }
    }
    ngfLagrange(c, S->mass, 1, it, sT, &mass);
    utcSec = S->utcSecAtBin[bT] + (t - (S->tb0 + NGF_BIN_S * bT));

    /* 1. attitude (geom.AttitudeProvider) */
    status = (S->attitude == NGF_ATT_PROVIDER) ? ngfAttitudeProvider(c, t, Rb) : ngfAttitudeNominal(S->alongVelocity, r, v, Rb);
    if (status != NGF_OK) {
        return ngfFail(c, status, "the attitude is not available at t (outside the nodes, or r x v = 0 for nominal)");
    }
    /* 2. geodetic position and relative velocity (env.wind) */
    mul33t(Ri, r, rI);
    if (ngfGeodetic(c, rI, lla) != NGF_OK) {
        return ngfFail(c, NGF_DOMAIN, "the position is near the geocentre");
    }
    for (i = 0; i < 3; i++) { vrel[i] = v[i]; }
    if (S->wind != NGF_WIND_NONE) {
        double corot[3];
        mul33(Rd, rI, corot);
        for (i = 0; i < 3; i++) { vrel[i] = v[i] - corot[i]; }
    }
    if (S->wind == NGF_WIND_HWM14) {
        double w[2], sl, cl, so, co, north[3], east[3], wI[3], wG[3];
        if (lla[2] < 0.0) {
            return ngfFail(c, NGF_DOMAIN, "hwm14 below the ellipsoid (the model does not return there)");
        }
        if (!c->hwm14Ready) {
            char miss[64] = "";
            if (ngfHwm14Load(c->hwm14Dir, miss, (int)sizeof miss) == 1) {
                return ngfFail(c, NGF_MODEL_FAIL, "an HWM14 data file is missing beside the gateway");
            }
            c->hwm14Ready = 1;
        }
        ngf_hwm14_calc((int)(S->yy[bT] * 1000.0 + S->doy[bT]), utcSec, lla[2] / 1000.0, deg(lla[0]), deg(lla[1]), S->ap3[bA], w);
        sl = sin(lla[0]); cl = cos(lla[0]); so = sin(lla[1]); co = cos(lla[1]);
        north[0] = -sl * co; north[1] = -sl * so; north[2] = cl;
        east[0] = -so; east[1] = co; east[2] = 0.0;
        for (i = 0; i < 3; i++) { wI[i] = w[0] * north[i] + w[1] * east[i]; }
        mul33(Ri, wI, wG);
        for (i = 0; i < 3; i++) { vrel[i] = vrel[i] - wG[i]; }
    }
    mul33t(Rb, vrel, vB);
    /* 3. density (env.density) */
    for (j = 0; j < NGF_SPECIES; j++) { sp[j] = 0.0; }
    if (S->density == NGF_DENS_EXPONENTIAL) {
        rho = S->expRho0 * exp(-(lla[2] - S->expH0) / S->expH);
    } else if (S->density == NGF_DENS_MSIS00) {
        struct nrlmsise_input in;
        struct nrlmsise_output out;
        struct nrlmsise_flags fl;
        struct ap_array aph;
        double latd = deg(lla[0]), lond = deg(lla[1]);
        fl.switches[0] = 0;
        for (i = 1; i < 24; i++) { fl.switches[i] = 1; }
        fl.switches[9] = -1;
        for (j = 0; j < 7; j++) { aph.a[j] = S->ap7[7 * bA + j]; }
        in.year = 0; in.doy = (int)S->doy[bT]; in.sec = utcSec; in.alt = lla[2] / 1000.0;
        in.g_lat = latd; in.g_long = lond; in.lst = utcSec / 3600.0 + lond / 15.0;
        in.f107A = S->f107a[bA]; in.f107 = S->f107[bA]; in.ap = aph.a[0]; in.ap_a = &aph;
        gtd7d(&in, &fl, &out);
        rho = out.d[5] * 1e3;
        sp[0] = out.d[0] * 1e6; sp[1] = out.d[1] * 1e6; sp[2] = out.d[2] * 1e6; sp[3] = out.d[3] * 1e6;
        sp[4] = out.d[4] * 1e6; sp[5] = out.d[6] * 1e6; sp[6] = out.d[7] * 1e6; sp[7] = out.d[8] * 1e6;
        Tg = out.t[1];
    } else {
        double dn[10], tn, tex, apk[7];
        if (!c->msis21Ready) {
            int st;
            ngf_msis21_init(c->msis21Dir, (int)strlen(c->msis21Dir), 1, 1, &st);
            if (st != 0) {
                return ngfFail(c, NGF_MODEL_FAIL, "NRLMSIS 2.1 refused its parameter folder");
            }
            c->msis21Ready = 1;
        }
        for (j = 0; j < 7; j++) { apk[j] = S->ap7[7 * bA + j]; }
        ngf_msis21_calc(S->doy[bT], utcSec, lla[2] / 1000.0, deg(lla[0]), deg(lla[1]), S->f107a[bA], S->f107[bA], apk, &tn, dn, &tex);
        for (j = 0; j < 10; j++) {
            if (dn[j] == 9.999e-38) { dn[j] = 0.0; }
        }
        rho = dn[0];
        sp[0] = dn[4]; sp[1] = dn[3]; sp[2] = dn[1]; sp[3] = dn[2];
        sp[4] = dn[6]; sp[5] = dn[5]; sp[6] = dn[7]; sp[7] = dn[8];
        Tg = tn;
    }
    /* 4. radiation sources (env.radiationSources), GCRS */
    {
        double d[3], nd, nrs = norm3(sun), sHat[3], Searth, zero[3] = {0.0, 0.0, 0.0};
        for (i = 0; i < 3; i++) { d[i] = sun[i] - r[i]; sHat[i] = sun[i] / nrs; }
        nd = norm3(d);
        nu = (S->eclipse == NGF_ECL_SOLAARS) ? ngfSolaars(r, sun, ecl, S->aEarth, S->bEarth)
                                             : ngfBodyShadow(S->eclipse, r, sun, zero, S->aEarth, S->Rsun);
        if (S->lunar) {
            nu = nu * ngfBodyShadow((S->eclipse == NGF_ECL_SOLAARS) ? NGF_ECL_CONICAL : S->eclipse, r, sun, moon, S->Rmoon, S->Rsun);
        }
        for (i = 0; i < 3; i++) { c->srcDir[i] = d[i] / nd; }
        c->srcP[0] = S->tsi / S->cLight * ((S->au / nd) * (S->au / nd)) * nu;
        c->band[0] = NGF_BAND_VIS;
        Searth = S->tsi * ((S->au / nrs) * (S->au / nrs));
        M = 1;
        if (S->erp == NGF_ERP_POINT) {
            double nr = norm3(r), rh[3], cosT, th, psi, geo;
            for (i = 0; i < 3; i++) { rh[i] = r[i] / nr; }
            cosT = fmax(-1.0, fmin(1.0, dot3(rh, sHat)));
            th = acos(cosT);
            psi = (sin(th) + (NGF_PI - th) * cos(th)) / NGF_PI;
            geo = (S->aEarth / nr) * (S->aEarth / nr);
            for (i = 0; i < 3; i++) { c->srcDir[3 + i] = -rh[i]; c->srcDir[6 + i] = -rh[i]; }
            c->srcP[1] = ((2.0 / 3.0) * 0.34 * Searth * geo * psi) / S->cLight;
            c->srcP[2] = (0.68 * (Searth / 4.0) * geo) / S->cLight;
            c->band[1] = NGF_BAND_VIS; c->band[2] = NGF_BAND_IR;
            M = 3;
        } else if (S->erp == NGF_ERP_KNOCKE || S->erp == NGF_ERP_CERES) {
            int Mt = ngfTiles(c, r, Ri, Searth, S->utcDay[bT], utcSec), ti;
            double Rt = S->aEarth + ((S->erp == NGF_ERP_CERES) ? S->ceresHeight : 0.0);
            if (Mt < 0) {
                return ngfFail(c, -Mt, "the Earth cap cannot be tiled here (inside the reference sphere, or outside the CERES months)");
            }
            for (ti = 0; ti < Mt; ti++) {
                const double *nt = c->tileN + 3 * ti;
                double rho3[3], dist2, dist, rhoHat[3], mu, muS, g;
                for (i = 0; i < 3; i++) { rho3[i] = r[i] - Rt * nt[i]; }
                dist2 = rho3[0] * rho3[0] + rho3[1] * rho3[1] + rho3[2] * rho3[2];
                dist = sqrt(dist2);
                for (i = 0; i < 3; i++) { rhoHat[i] = rho3[i] / dist; }
                mu = fmax(0.0, dot3(nt, rhoHat));
                muS = fmax(0.0, dot3(nt, sHat));
                g = mu * c->tileA[ti] / (NGF_PI * dist2);
                for (i = 0; i < 3; i++) {
                    c->srcDir[3 * (1 + ti) + i] = -rhoHat[i];
                    c->srcDir[3 * (1 + Mt + ti) + i] = -rhoHat[i];
                }
                c->srcP[1 + ti] = (c->tileAlb[ti] * Searth * muS * g) / S->cLight;
                c->srcP[1 + Mt + ti] = (c->tileLw[ti] * g) / S->cLight;
                c->band[1 + ti] = NGF_BAND_VIS;
                c->band[1 + Mt + ti] = NGF_BAND_IR;
            }
            M = 1 + 2 * Mt;
        }
        for (m = 0; m < M; m++) {                                /* to the body frame */
            double dg[3];
            for (i = 0; i < 3; i++) { dg[i] = c->srcDir[3 * m + i]; }
            mul33t(Rb, dg, c->srcDir + 3 * m);
        }
    }
    /* 5. occlusion (Model.occlusion) */
    {
        int kinetic = (S->gsi != NGF_GSI_CONSTANT_CD);
        double U = norm3(vB), vu[3];
        if (S->shadow == NGF_SHADOW_LUT) {
            for (i = 0; i < 3; i++) { vu[i] = vB[i] / U; }
            for (k = 0; k < N; k++) { c->beam[k] = ngfLut(c, k, vu); }
            if (kinetic && !(S->flowTable && ngfFlowTable(c, vB, U, Tg))) {
                ngfFlowShadow(c, vB, U, Tg);
            }
            for (m = 0; m < M; m++) {
                for (k = 0; k < N; k++) {
                    c->fRad[N * m + k] = ngfLut(c, k, c->srcDir + 3 * m);
                }
            }
        } else {
            for (m = 0; m < N * M; m++) { c->fRad[m] = 1.0; }
        }
    }
    /* 6. panel temperatures (geom.thermalState, steady state) */
    if (S->thermal == NGF_THERMAL_STEADY) {
        for (k = 0; k < N; k++) {
            const double *nk = S->normal + 3 * k;
            double q = 0.0;
            for (m = 0; m < M; m++) {
                double ab = (c->band[m] == NGF_BAND_VIS) ? S->optVis[3 * k] : S->optIR[3 * k];
                double cosI = fmax(0.0, dot3(nk, c->srcDir + 3 * m));
                q = q + ((ab * cosI) * (S->cLight * c->srcP[m])) * c->fRad[N * m + k];
            }
            if (S->internalPower > 0.0 && S->internalShare != NULL) {
                q = q + S->internalPower * S->internalShare[k] / S->area[k];
            }
            c->tPanel[k] = pow(q / (S->optIR[3 * k] * S->sigmaSB), 0.25);
        }
    }
    /* 7-8. per-panel forces and their reduction (Model.environment, reduceAccel) */
    {
        ngfSources src;
        src.m = M; src.dir = c->srcDir; src.pres = c->srcP; src.band = c->band;
        status = ngfDrag(c, rho, sp, Tg, vB, S->shadow == NGF_SHADOW_LUT, c->Fk);
        if (status != NGF_OK) { return ngfFail(c, status, "the gas-surface model refused its inputs"); }
        for (i = 0; i < 3 * N; i++) { c->Fk[i] = S->cdScale * c->Fk[i]; }
        status = ngfAssemble(&c->macro, c->Fk, S->rcm, FdSum, tau);
        if (status != NGF_OK) { return ngfFail(c, status, "assemble (drag)"); }
        status = ngfRadiationPanelForce(&c->macro, &c->con, &src, c->fRad, c->Fk, NULL);
        if (status != NGF_OK) { return ngfFail(c, status, "the radiation kernel refused its inputs"); }
        for (i = 0; i < 3 * N; i++) { c->Fk[i] = S->crScale * c->Fk[i]; }
        status = ngfAssemble(&c->macro, c->Fk, S->rcm, FrSum, tau);
        if (status != NGF_OK) { return ngfFail(c, status, "assemble (radiation)"); }
        if (S->thermal == NGF_THERMAL_STEADY) {
            status = ngfThermalRecoil(&c->macro, &c->con, c->tPanel, c->Fk, NULL);
            if (status != NGF_OK) { return ngfFail(c, status, "the thermal kernel refused its inputs"); }
            status = ngfAssemble(&c->macro, c->Fk, S->rcm, FtSum, tau);
            if (status != NGF_OK) { return ngfFail(c, status, "assemble (thermal)"); }
        } else {
            FtSum[0] = 0.0; FtSum[1] = 0.0; FtSum[2] = 0.0;
        }
    }
    for (i = 0; i < 3; i++) {
        aB[i] = FdSum[i] / mass + FrSum[i] / mass + FtSum[i] / mass;
    }
    mul33(Rb, aB, aG);
    if (pt != NULL) {
        for (i = 0; i < 3; i++) {
            pt->aD[i] = FdSum[i] / mass; pt->aR[i] = FrSum[i] / mass; pt->aT[i] = FtSum[i] / mass;
            pt->vB[i] = vB[i]; pt->lla[i] = lla[i];
        }
        pt->rho = rho; pt->Tgas = Tg; pt->nu = nu; pt->mass = mass;
        for (j = 0; j < NGF_SPECIES; j++) { pt->species[j] = sp[j]; }
        for (i = 0; i < 9; i++) { pt->Rb[i] = Rb[i]; }
    }
    return NGF_OK;
}

/* ========================================================================= public API */
static int ngfEntryAccel(double t, const double *R, const double *V, double *A, int n, int m, void *user)
{
    ngfRhsContext *c = (ngfRhsContext *)user;
    int status;
    if (n != 3 || m < 1 || m > c->s.maxColumns) {
        return GJ_NATIVE_ESHAPE;
    }
    status = ngfRhsAccel(c, t, R, V, A, m);
    if (status == NGF_OK) {
        return GJ_NATIVE_OK;
    }
    return (status == NGF_DOMAIN) ? GJ_NATIVE_EDOMAIN : GJ_NATIVE_EFAIL;
}

int ngfRhsCreate(const ngfRhsSetup *setup, ngfRhsContext **out)
{
    ngfRhsContext *c;
    ngfRhsSetup *S;
    int st = NGF_OK, N, Mt, k, cells;
    size_t nT, nb, grid;
    if (setup == NULL || out == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    *out = NULL;
    c = (ngfRhsContext *)calloc(1, sizeof(ngfRhsContext));
    if (c == NULL) {
        return NGF_ALLOC;
    }
    c->s = *setup;
    S = &c->s;
    N = S->n;
    if (N < 1 || S->nT < 4 || S->nb < 1 || !(S->dt > 0.0) || S->maxColumns < 1) {
        free(c);
        return NGF_BAD_RANGE;
    }
    nT = (size_t)S->nT;
    nb = (size_t)S->nb;
    st |= ngfOwnCopy(c, &S->area, (size_t)N, 1);
    st |= ngfOwnCopy(c, &S->normal, 3 * (size_t)N, 1);
    st |= ngfOwnCopy(c, &S->centroid, 3 * (size_t)N, 1);
    st |= ngfOwnCopy(c, &S->optVis, 3 * (size_t)N, 1);
    st |= ngfOwnCopy(c, &S->optIR, 3 * (size_t)N, 1);
    st |= ngfOwnCopy(c, &S->internalShare, (size_t)N, 0);
    cells = 6 * (S->lutResolution + 1) * (S->lutResolution + 1);
    st |= ngfOwnCopy(c, &S->lut, (S->shadow == NGF_SHADOW_LUT) ? (size_t)N * (size_t)cells : 0, S->shadow == NGF_SHADOW_LUT);
    st |= ngfOwnCopy(c, &S->xChi, (size_t)S->nChi, S->shadow == NGF_SHADOW_LUT);
    st |= ngfOwnCopy(c, &S->wChi, (size_t)S->nChi, S->shadow == NGF_SHADOW_LUT);
    st |= ngfOwnCopy(c, &S->xArc, (size_t)S->nArc, S->shadow == NGF_SHADOW_LUT);
    st |= ngfOwnCopy(c, &S->wArc, (size_t)S->nArc, S->shadow == NGF_SHADOW_LUT);
    {
        size_t ftRows = 24 * (size_t)S->ftNShaded, ftNodes = (size_t)S->ftN * (size_t)S->ftN * (size_t)S->ftNs;
        int ft = (S->flowTable != 0);
        st |= ngfOwnCopy(c, &S->ftShaded, ft ? (size_t)S->ftNShaded : 0, ft);
        st |= ngfOwnCopy(c, &S->ftPlus, ft ? ftRows * ftNodes : 0, ft);
        st |= ngfOwnCopy(c, &S->ftMinus, ft ? ftRows * ftNodes : 0, ft);
    }
    st |= ngfOwnCopy(c, &S->R, 9 * nT, 1);
    st |= ngfOwnCopy(c, &S->Rdot, 9 * nT, 1);
    st |= ngfOwnCopy(c, &S->sun, 3 * nT, 1);
    st |= ngfOwnCopy(c, &S->moon, S->lunar ? 3 * nT : 0, S->lunar);
    st |= ngfOwnCopy(c, &S->eclipticNormal, (S->eclipse == NGF_ECL_SOLAARS) ? 3 * nT : 0, S->eclipse == NGF_ECL_SOLAARS);
    st |= ngfOwnCopy(c, &S->mass, nT, 1);
    st |= ngfOwnCopy(c, &S->tq, (S->attitude == NGF_ATT_PROVIDER) ? (size_t)S->nq : 0, S->attitude == NGF_ATT_PROVIDER);
    st |= ngfOwnCopy(c, &S->q, (S->attitude == NGF_ATT_PROVIDER) ? 4 * (size_t)S->nq : 0, S->attitude == NGF_ATT_PROVIDER);
    st |= ngfOwnCopy(c, &S->alignedEdges, nb + 1, 1);
    st |= ngfOwnCopy(c, &S->f107, nb, 1);
    st |= ngfOwnCopy(c, &S->f107a, nb, 1);
    st |= ngfOwnCopy(c, &S->ap7, 7 * nb, 1);
    st |= ngfOwnCopy(c, &S->ap3, nb, 1);
    st |= ngfOwnCopy(c, &S->doy, nb, 1);
    st |= ngfOwnCopy(c, &S->yy, nb, 1);
    st |= ngfOwnCopy(c, &S->utcSecAtBin, nb, 1);
    st |= ngfOwnCopy(c, &S->utcDay, nb, 1);
    grid = (S->erp == NGF_ERP_CERES) ? (size_t)S->nLat * (size_t)S->nLon * (size_t)S->nMonth : 0;
    st |= ngfOwnCopy(c, &S->latC, (S->erp == NGF_ERP_CERES) ? (size_t)S->nLat : 0, S->erp == NGF_ERP_CERES);
    st |= ngfOwnCopy(c, &S->lonC, (S->erp == NGF_ERP_CERES) ? (size_t)S->nLon : 0, S->erp == NGF_ERP_CERES);
    st |= ngfOwnCopy(c, &S->tMonth, (S->erp == NGF_ERP_CERES) ? (size_t)S->nMonth : 0, S->erp == NGF_ERP_CERES);
    st |= ngfOwnCopy(c, &S->rsut, grid, S->erp == NGF_ERP_CERES);
    st |= ngfOwnCopy(c, &S->rlut, grid, S->erp == NGF_ERP_CERES);
    st |= ngfOwnCopy(c, &S->rsdt, grid, S->erp == NGF_ERP_CERES);
    if (st != NGF_OK) {
        ngfRhsDestroy(c);
        return NGF_NULL_ARGUMENT;
    }
    snprintf(c->msis21Dir, sizeof c->msis21Dir, "%s", S->msis21Dir ? S->msis21Dir : "");
    snprintf(c->hwm14Dir, sizeof c->hwm14Dir, "%s", S->hwm14Dir ? S->hwm14Dir : "");
    S->msis21Dir = c->msis21Dir;
    S->hwm14Dir = c->hwm14Dir;
    Mt = (S->erp == NGF_ERP_KNOCKE || S->erp == NGF_ERP_CERES) ? S->nRings * S->nSectors : 1;
    c->mTiles = Mt;
    c->mMax = 1 + 2 * Mt + 2;
    c->ringMax = 6 * ((S->nChi > 0) ? S->nChi : 1);
    c->srcDir = ngfOwnZeros(c, 3 * (size_t)c->mMax);
    c->srcP = ngfOwnZeros(c, (size_t)c->mMax);
    c->fRad = ngfOwnZeros(c, (size_t)N * (size_t)c->mMax);
    c->tileN = ngfOwnZeros(c, 3 * (size_t)Mt);
    c->tileA = ngfOwnZeros(c, (size_t)Mt);
    c->tileAlb = ngfOwnZeros(c, (size_t)Mt);
    c->tileLw = ngfOwnZeros(c, (size_t)Mt);
    c->fN = ngfOwnZeros(c, (size_t)N * NGF_SPECIES);
    c->fP = ngfOwnZeros(c, (size_t)N * NGF_SPECIES);
    c->fT = ngfOwnZeros(c, (size_t)N * NGF_SPECIES);
    c->ring = ngfOwnZeros(c, 8 * (size_t)c->ringMax);
    c->beam = ngfOwnZeros(c, (size_t)N);
    c->tPanel = ngfOwnZeros(c, (size_t)N);
    c->Fk = ngfOwnZeros(c, 3 * (size_t)N);
    c->baseA = ngfOwnZeros(c, 3 * (size_t)S->maxColumns);
    c->band = (int *)calloc((size_t)c->mMax, sizeof(int));
    c->lutAllOnes = (int *)calloc((size_t)N, sizeof(int));
    if (c->srcDir == NULL || c->srcP == NULL || c->fRad == NULL || c->tileN == NULL || c->tileA == NULL ||
        c->tileAlb == NULL || c->tileLw == NULL || c->fN == NULL || c->fP == NULL || c->fT == NULL ||
        c->ring == NULL || c->beam == NULL || c->tPanel == NULL || c->Fk == NULL || c->baseA == NULL ||
        c->band == NULL || c->lutAllOnes == NULL) {
        ngfRhsDestroy(c);
        return NGF_ALLOC;
    }
    if (S->shadow == NGF_SHADOW_LUT) {
        for (k = 0; k < N; k++) {
            int cc, all = 1;
            for (cc = 0; cc < cells && all; cc++) {
                if (S->lut[N * cc + k] != 1.0) { all = 0; }
            }
            c->lutAllOnes[k] = all;
        }
    }
    c->macro.n = N; c->macro.area = S->area; c->macro.normal = S->normal; c->macro.centroid = S->centroid;
    c->macro.optVis = S->optVis; c->macro.optIR = S->optIR;
    c->con.cLight = S->cLight; c->con.sigmaSB = S->sigmaSB; c->con.kB = S->kB;
    for (k = 0; k < NGF_SPECIES; k++) { c->con.speciesMass[k] = S->speciesMass[k]; }
    c->entry.abi_version = GJ_NATIVE_ABI_VERSION;
    c->entry.n = 3;
    c->entry.m_max = S->maxColumns;
    c->entry.accel = ngfEntryAccel;
    c->entry.user = c;
    *out = c;
    return NGF_OK;
}

void ngfRhsDestroy(ngfRhsContext *c)
{
    int i;
    if (c == NULL) {
        return;
    }
    for (i = 0; i < c->nOwn; i++) {
        free(c->own[i]);
    }
    free(c->band);
    free(c->lutAllOnes);
    free(c);
}

int ngfRhsAccel(ngfRhsContext *c, double t, const double *R, const double *V, double *A, int m)
{
    int j, i, status;
    double a[3];
    if (c == NULL || R == NULL || V == NULL || A == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    if (c->base != NULL) {
        status = gj_native_call(c->base, t, R, V, c->baseA, 3, m);
        if (status != GJ_NATIVE_OK) {
            return ngfFail(c, (status == GJ_NATIVE_EDOMAIN) ? NGF_DOMAIN : NGF_MODEL_FAIL, "the base (gravity) entry refused the evaluation");
        }
    }
    for (j = 0; j < m; j++) {
        status = ngfEvalOne(c, t, R + 3 * j, V + 3 * j, a, NULL);
        if (status != NGF_OK) {
            return status;
        }
        for (i = 0; i < 3; i++) {
            A[3 * j + i] = (c->base != NULL) ? c->baseA[3 * j + i] + a[i] : a[i];
        }
    }
    return NGF_OK;
}

int ngfRhsDiagnose(ngfRhsContext *c, double t, const double *r, const double *v, double *diag)
{
    ngfParts p;
    double a[3];
    int status, i;
    if (c == NULL || r == NULL || v == NULL || diag == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    status = ngfEvalOne(c, t, r, v, a, &p);
    if (status != NGF_OK) {
        return status;
    }
    for (i = 0; i < 3; i++) {
        diag[i] = p.aD[i]; diag[3 + i] = p.aR[i]; diag[6 + i] = p.aT[i];
        diag[19 + i] = p.vB[i]; diag[22 + i] = p.lla[i];
    }
    diag[9] = p.rho;
    for (i = 0; i < NGF_SPECIES; i++) { diag[10 + i] = p.species[i]; }
    diag[18] = p.Tgas;
    diag[25] = p.nu;
    diag[26] = p.mass;
    for (i = 0; i < 9; i++) { diag[27 + i] = p.Rb[i]; }
    return NGF_OK;
}

void *ngfRhsEntry(ngfRhsContext *c)
{
    return (c == NULL) ? NULL : (void *)&c->entry;
}

void ngfRhsSetBase(ngfRhsContext *c, void *baseEntry)
{
    if (c != NULL) {
        c->base = (const gj_native_entry *)baseEntry;
    }
}

const char *ngfRhsLastError(const ngfRhsContext *c)
{
    return (c == NULL) ? "no context" : c->lastError;
}
