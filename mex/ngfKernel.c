/*
 * ngfKernel.c
 * -----------
 * Implementation of the single-epoch NonGravForces kernel (LP-6b, WP-17b).
 *
 * HOW TO READ THIS FILE
 *   Every block below is a transcription of one MATLAB unit, named in its
 *   comment, and the arithmetic is written in the SAME ASSOCIATION ORDER as
 *   the MATLAB. That is deliberate and it is the reason the two agree to
 *   rounding rather than to "about 1e-12": MATLAB's `a .* b / c` is (a*b)/c,
 *   `2 * x * y` is (2*x)*y, and a re-association here would move the last two
 *   digits for no gain. Where the MATLAB writes a form chosen for
 *   conditioning -- erfc(-x) instead of 1 + erf(x), (1 - alpha) + alpha*b
 *   instead of 1 + alpha*(b - 1) -- the same form is written here, with the
 *   same reason.
 *
 *   The better-conditioned leeward form recorded in HANDOVER O-20 (continued
 *   fraction for erfc, Veltkamp/Dekker split) is deliberately NOT used: MATLAB
 *   is normative, and a C file that is more accurate than the unit it ports
 *   fails the parity test it exists to pass. O-20 stays where it is.
 *
 * NO MATLAB HEADERS. This file compiles with any C99 compiler and is exercised
 * without MATLAB by mex/test/ngfKernelParity.c.
 *
 * -------------------------------------------------------------------------
 * NonGravForces 0.6.1
 * Created 2026-09-08 by Claude Opus 5 (Anthropic) for Matthias Weigelt.
 * -------------------------------------------------------------------------
 */
#include "ngfKernel.h"

#include <math.h>
#include <stddef.h>

#define NGF_VERSION_STRING "0.6.1"

/* pi to double precision. NGF_PI is a POSIX extension: it is absent under a strict
 * -std=c99 and under MSVC without _USE_MATH_DEFINES, so it is not relied on. */
#define NGF_PI 3.14159265358979323846

/* MATLAB's realmin, the divisor guard in dragPanelForce's info.Cp. */
#define NGF_REALMIN 2.2250738585072014e-308

static int ngfFinite(double x)
{
    return (x == x) && (x > -HUGE_VAL) && (x < HUGE_VAL);
}

static int ngfFiniteArray(const double *x, int n)
{
    int i;
    if (x == NULL) {
        return 1;                       /* an absent optional array is finite */
    }
    for (i = 0; i < n; i++) {
        if (!ngfFinite(x[i])) {
            return 0;
        }
    }
    return 1;
}

static int ngfMacroOk(const ngfMacro *macro)
{
    if (macro == NULL || macro->area == NULL || macro->normal == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    if (macro->n < 1) {
        return NGF_BAD_PLATE_COUNT;
    }
    return NGF_OK;
}

const char *ngfKernelVersion(void)
{
    return NGF_VERSION_STRING;
}

const char *ngfStatusText(int status)
{
    switch (status) {
        case NGF_OK:               return "ok";
        case NGF_NULL_ARGUMENT:    return "a required pointer was NULL";
        case NGF_BAD_PLATE_COUNT:  return "the plate count must be at least 1";
        case NGF_BAD_MODEL:        return "the gas-surface model is not sentman, maxwell or constantCd";
        case NGF_BAD_RANGE:        return "an input is outside its permitted range";
        case NGF_NOT_FINITE:       return "an input is NaN or Inf";
        case NGF_BAD_SOURCE_COUNT: return "the source count is negative or a band is not vis|ir";
        default:                   return "unknown status";
    }
}

/* =========================================================================
 * gsiCoefficients  --  port of +NonGravForces/+kernel/gsiCoefficients.m
 * ========================================================================= */
int ngfGsiParts(ngfGsiModel model, double s, double theta, double alpha,
                double tw, double tinf, double cd, double specular,
                double *cpIncOut, double *cpReOut, double *ctau)
{
    double g, l, P, Z, cpInc, ctauInc, vRatio, cpRe, gmax;

    if (cpIncOut == NULL || cpReOut == NULL || ctau == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    if (model != NGF_GSI_SENTMAN && model != NGF_GSI_MAXWELL &&
        model != NGF_GSI_CONSTANT_CD) {
        return NGF_BAD_MODEL;
    }
    if (!ngfFinite(theta)) {
        return NGF_NOT_FINITE;
    }
    if (theta < 0.0 || theta > NGF_PI) {
        return NGF_BAD_RANGE;
    }

    g = cos(theta);
    l = sin(theta);

    if (model == NGF_GSI_CONSTANT_CD) {
        if (!ngfFinite(cd)) {
            return NGF_NOT_FINITE;
        }
        if (!(cd > 0.0)) {
            return NGF_BAD_RANGE;
        }
        gmax = (g > 0.0) ? g : 0.0;                 /* MATLAB max(g, 0)       */
        *cpIncOut = (cd * gmax) * g;
        *cpReOut  = 0.0;                            /* info.CpReemitted = zeros */
        *ctau     = (cd * gmax) * l;
        return NGF_OK;
    }

    if (!ngfFinite(s) || !ngfFinite(alpha) || !ngfFinite(tw) || !ngfFinite(tinf)) {
        return NGF_NOT_FINITE;
    }
    if (!(s > 0.0) || alpha < 0.0 || alpha > 1.0 || !(tw > 0.0) || !(tinf > 0.0)) {
        return NGF_BAD_RANGE;
    }
    if (model == NGF_GSI_MAXWELL) {
        if (!ngfFinite(specular)) {
            return NGF_NOT_FINITE;
        }
        if (specular < 0.0 || specular > 1.0) {
            return NGF_BAD_RANGE;
        }
    }

    /* P = exp(-(g s)^2) / s ; Z = erfc(-g s) -- erfc, not 1 + erf, because the
     * sum cancels to zero or below on a leeward plate (g s << 0). */
    P = exp(-((g * s) * (g * s))) / s;
    Z = erfc(-(g * s));

    cpInc   = (g * P) / sqrt(NGF_PI) + ((g * g) + 1.0 / (2.0 * (s * s))) * Z;
    ctauInc = l * (P / sqrt(NGF_PI) + g * Z);

    /* vRatio = sqrt((1/2)[(1 - alpha) + alpha 2 Tw/(Tinf s^2)]) -- written with
     * (1 - alpha) rather than 1 + alpha(b - 1), which loses 1e-8 at alpha = 1. */
    vRatio = sqrt(0.5 * ((1.0 - alpha) + ((alpha * 2.0) * tw) / (tinf * (s * s))));
    cpRe   = (0.5 * vRatio) * ((g * sqrt(NGF_PI)) * Z + P);

    if (model == NGF_GSI_SENTMAN) {
        *cpIncOut = cpInc;
        *cpReOut  = cpRe;
        *ctau     = ctauInc;
    } else {
        /* LP-C WP-C1: the diffuse part (1 - sig) carries incident and re-emitted
         * momentum, the specular part doubles the incident normal momentum. */
        *cpIncOut = (1.0 + specular) * cpInc;
        *cpReOut  = (1.0 - specular) * cpRe;
        *ctau     = (1.0 - specular) * ctauInc;
    }
    return NGF_OK;
}

int ngfGsiCoefficients(ngfGsiModel model, double s, double theta, double alpha,
                       double tw, double tinf, double cd, double specular,
                       double *cp, double *ctau)
{
    double inc, re;
    int status;
    if (cp == NULL || ctau == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    status = ngfGsiParts(model, s, theta, alpha, tw, tinf, cd, specular, &inc, &re, ctau);
    if (status != NGF_OK) {
        return status;
    }
    *cp = inc + re;                                 /* Cp = CpIncident + CpReemitted */
    return NGF_OK;
}

/* =========================================================================
 * dragPanelForce  --  port of +NonGravForces/+kernel/dragPanelForce.m
 *
 * The species loop is Doornbos (2011) eq. 3.57: the mixture coefficient is the
 * mass-weighted mean of the per-species coefficients, each with its own speed
 * ratio. Treating the gas as air biases s by sqrt(29/16) = 1.35 at 500 km
 * (HANDOVER defect D1b), so there is no mean-molar-mass shortcut here.
 * ========================================================================= */
int ngfDragPanelForce(const ngfMacro *macro, const ngfConstants *con,
                      double rho, const double *species, double tGas,
                      const double *vrel, ngfGsiModel model, double alpha,
                      double tw, double cd, double specular,
                      const double *shadow, double *F,
                      double *gammaOut, double *thetaOut, double *sRatioOut,
                      double *cpMixOut, double *ctauMixOut)
{
    int status, k, j, n, i;
    double U, v[3], nk[3], cr[3], t[3], gk, lk, thetaK;
    double sRatio[NGF_SPECIES], rhoJ, rhoSum, cpJ, ctJ, cpMix, ctMix;
    double scale, fn, ft, f;

    status = ngfMacroOk(macro);
    if (status != NGF_OK) {
        return status;
    }
    if (con == NULL || vrel == NULL || F == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    if (model != NGF_GSI_CONSTANT_CD && species == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    n = macro->n;
    if (!ngfFiniteArray(vrel, 3) || !ngfFinite(rho) || !ngfFiniteArray(shadow, n)) {
        return NGF_NOT_FINITE;
    }
    if (model != NGF_GSI_CONSTANT_CD &&
        (!ngfFiniteArray(species, NGF_SPECIES) || !ngfFinite(tGas))) {
        return NGF_NOT_FINITE;
    }
    if (rho < 0.0) {
        return NGF_BAD_RANGE;
    }
    if (model != NGF_GSI_CONSTANT_CD && !(tGas > 0.0)) {
        return NGF_BAD_RANGE;
    }
    if (shadow != NULL) {
        for (k = 0; k < n; k++) {
            if (shadow[k] < 0.0 || shadow[k] > 1.0) {
                return NGF_BAD_RANGE;
            }
        }
    }

    U = sqrt(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]);
    if (!(U > 0.0)) {
        return NGF_BAD_RANGE;
    }
    for (i = 0; i < 3; i++) {
        v[i] = vrel[i] / U;
    }

    /* Speed ratio per species, s_j = U / sqrt(2 k T / m_j). NaN for constantCd,
     * which the MATLAB reports the same way. */
    for (j = 0; j < NGF_SPECIES; j++) {
        if (model == NGF_GSI_CONSTANT_CD) {
            sRatio[j] = NAN;
        } else {
            sRatio[j] = U / sqrt(((2.0 * con->kB) * tGas) / con->speciesMass[j]);
        }
    }
    if (sRatioOut != NULL) {
        for (j = 0; j < NGF_SPECIES; j++) {
            sRatioOut[j] = sRatio[j];
        }
    }

    /* Mass density re-formed from the species, exactly as the MATLAB does; it
     * differs from the backend's own total by 4e-4 relative (measured
     * 2026-09-07), which is a statement about the reported total, not the force. */
    rhoSum = 0.0;
    if (model != NGF_GSI_CONSTANT_CD) {
        for (j = 0; j < NGF_SPECIES; j++) {
            rhoSum = rhoSum + species[j] * con->speciesMass[j];
        }
    } else {
        rhoSum = rho;
    }

    for (k = 0; k < n; k++) {
        for (i = 0; i < 3; i++) {
            nk[i] = macro->normal[3 * k + i];
        }
        /* gamma = n.v ; l = |n x v| ; theta = atan2(l, gamma). The atan2 form is
         * used because acos of a dot product loses sqrt(eps) near 0 and 180 deg
         * (a defect this project met twice). */
        gk = nk[0] * v[0] + nk[1] * v[1] + nk[2] * v[2];
        cr[0] = nk[1] * v[2] - nk[2] * v[1];
        cr[1] = nk[2] * v[0] - nk[0] * v[2];
        cr[2] = nk[0] * v[1] - nk[1] * v[0];
        lk = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
        thetaK = atan2(lk, gk);
        for (i = 0; i < 3; i++) {
            t[i] = -v[i] + gk * nk[i];
        }
        if (lk > 0.0) {
            for (i = 0; i < 3; i++) {
                t[i] = t[i] / lk;
            }
        } else {
            /* Edge-on or head-on: no tangential direction exists and Ctau -> 0
             * there anyway, so the direction is set to zero rather than left 0/0. */
            t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;
        }
        if (gammaOut != NULL) { gammaOut[k] = gk; }
        if (thetaOut != NULL) { thetaOut[k] = thetaK; }

        cpMix = 0.0;
        ctMix = 0.0;
        if (model == NGF_GSI_CONSTANT_CD) {
            status = ngfGsiCoefficients(model, 0.0, thetaK, 0.0, 0.0, 0.0,
                                        cd, specular, &cpJ, &ctJ);
            if (status != NGF_OK) {
                return status;
            }
            cpMix = rho * cpJ;
            ctMix = rho * ctJ;
        } else {
            for (j = 0; j < NGF_SPECIES; j++) {
                rhoJ = species[j] * con->speciesMass[j];
                if (!(rhoJ > 0.0)) {
                    continue;               /* MATLAB: if ~any(rhoj > 0), continue */
                }
                status = ngfGsiCoefficients(model, sRatio[j], thetaK, alpha,
                                            tw, tGas, cd, specular, &cpJ, &ctJ);
                if (status != NGF_OK) {
                    return status;
                }
                cpMix = cpMix + rhoJ * cpJ;
                ctMix = ctMix + rhoJ * ctJ;
            }
        }
        if (cpMixOut != NULL) {
            cpMixOut[k] = cpMix / ((rhoSum > NGF_REALMIN) ? rhoSum : NGF_REALMIN);
        }
        if (ctauMixOut != NULL) {
            ctauMixOut[k] = ctMix / ((rhoSum > NGF_REALMIN) ? rhoSum : NGF_REALMIN);
        }

        f = (shadow == NULL) ? 1.0 : shadow[k];
        scale = ((0.5 * (U * U)) * macro->area[k]) * f;
        fn = -(scale * cpMix);                       /* along -n              */
        ft = scale * ctMix;                          /* along t               */
        for (i = 0; i < 3; i++) {
            F[3 * k + i] = nk[i] * fn + t[i] * ft;
        }
    }
    return NGF_OK;
}

/* =========================================================================
 * radiationPanelForce  --  port of +kernel/radiationPanelForce.m
 *
 *   F_km = - P_m A_k f_km gamma [ (a + delta) d + (2 s gamma + (2/3) delta) n ]
 *
 * One expression serves sun, albedo and infrared through the band: sources, not
 * effects (principle P3). The eclipse factor is already inside P.
 * ========================================================================= */
int ngfRadiationPanelForce(const ngfMacro *macro, const ngfConstants *con,
                           const ngfSources *src, const double *shadow,
                           double *F, double *absorbed)
{
    int status, k, m, i, n, M;
    const double *opt;
    double dm[3], nk[3], gamma, gp, w, ad, cn, f, pm;

    status = ngfMacroOk(macro);
    if (status != NGF_OK) {
        return status;
    }
    if (con == NULL || F == NULL || macro->optVis == NULL || macro->optIR == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    n = macro->n;
    M = (src == NULL) ? 0 : src->m;
    if (M < 0) {
        return NGF_BAD_SOURCE_COUNT;
    }
    if (M > 0 && (src->dir == NULL || src->pres == NULL || src->band == NULL)) {
        return NGF_NULL_ARGUMENT;
    }
    if (M > 0 && (!ngfFiniteArray(src->dir, 3 * M) || !ngfFiniteArray(src->pres, M))) {
        return NGF_NOT_FINITE;
    }
    if (!ngfFiniteArray(shadow, n * M)) {
        return NGF_NOT_FINITE;
    }
    for (m = 0; m < M; m++) {
        if (src->band[m] != NGF_BAND_VIS && src->band[m] != NGF_BAND_IR) {
            return NGF_BAD_SOURCE_COUNT;
        }
        if (src->pres[m] < 0.0) {
            return NGF_BAD_RANGE;
        }
    }
    if (shadow != NULL) {
        for (i = 0; i < n * M; i++) {
            if (shadow[i] < 0.0 || shadow[i] > 1.0) {
                return NGF_BAD_RANGE;
            }
        }
    }

    for (i = 0; i < 3 * n; i++) {
        F[i] = 0.0;
    }
    if (absorbed != NULL) {
        for (k = 0; k < n; k++) {
            absorbed[k] = 0.0;
        }
    }

    for (m = 0; m < M; m++) {
        dm[0] = src->dir[3 * m];
        dm[1] = src->dir[3 * m + 1];
        dm[2] = src->dir[3 * m + 2];
        opt = (src->band[m] == NGF_BAND_VIS) ? macro->optVis : macro->optIR;
        pm  = src->pres[m];
        for (k = 0; k < n; k++) {
            nk[0] = macro->normal[3 * k];
            nk[1] = macro->normal[3 * k + 1];
            nk[2] = macro->normal[3 * k + 2];
            gamma = nk[0] * dm[0] + nk[1] * dm[1] + nk[2] * dm[2];
            gp = (gamma > 0.0) ? gamma : 0.0;
            f = (shadow == NULL) ? 1.0 : shadow[n * m + k];
            w = ((pm * macro->area[k]) * f) * gp;
            if (absorbed != NULL) {
                absorbed[k] = absorbed[k] + ((opt[3 * k] * w) * con->cLight) / macro->area[k];
            }
            ad = opt[3 * k] + opt[3 * k + 2];                    /* a + delta  */
            cn = (2.0 * opt[3 * k + 1]) * gp + (2.0 / 3.0) * opt[3 * k + 2];
            for (i = 0; i < 3; i++) {
                F[3 * k + i] = F[3 * k + i] - w * (ad * dm[i] + cn * nk[i]);
            }
        }
    }
    return NGF_OK;
}

/* =========================================================================
 * thermalRecoil  --  port of +kernel/thermalRecoil.m
 *   F_k = -(2/3) eps_k sigma T_k^4 A_k / c  n_k
 * ========================================================================= */
int ngfThermalRecoil(const ngfMacro *macro, const ngfConstants *con,
                     const double *tPanel, double *F, double *emitted)
{
    int status, k, i, n;
    double t2, t4, em, mag;

    status = ngfMacroOk(macro);
    if (status != NGF_OK) {
        return status;
    }
    if (con == NULL || tPanel == NULL || F == NULL || macro->optIR == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    n = macro->n;
    if (!ngfFiniteArray(tPanel, n)) {
        return NGF_NOT_FINITE;
    }
    for (k = 0; k < n; k++) {
        if (tPanel[k] < 0.0) {
            return NGF_BAD_RANGE;
        }
    }
    for (k = 0; k < n; k++) {
        t2 = tPanel[k] * tPanel[k];
        t4 = t2 * t2;                        /* MATLAB .^4: repeated squaring */
        em = (macro->optIR[3 * k] * con->sigmaSB) * t4;
        mag = (((2.0 / 3.0) * em) * macro->area[k]) / con->cLight;
        if (emitted != NULL) {
            emitted[k] = em;
        }
        for (i = 0; i < 3; i++) {
            F[3 * k + i] = -macro->normal[3 * k + i] * mag;
        }
    }
    return NGF_OK;
}

/* =========================================================================
 * assemble  --  port of +kernel/assemble.m
 * ========================================================================= */
int ngfAssemble(const ngfMacro *macro, const double *Fk, const double *rcm,
                double *F, double *tau)
{
    int status, k, i, n;
    double lever[3], f[3];

    status = ngfMacroOk(macro);
    if (status != NGF_OK) {
        return status;
    }
    if (Fk == NULL || rcm == NULL || F == NULL || tau == NULL ||
        macro->centroid == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    n = macro->n;
    if (!ngfFiniteArray(Fk, 3 * n) || !ngfFiniteArray(rcm, 3) ||
        !ngfFiniteArray(macro->centroid, 3 * n)) {
        return NGF_NOT_FINITE;
    }
    for (i = 0; i < 3; i++) {
        F[i] = 0.0;
        tau[i] = 0.0;
    }
    for (k = 0; k < n; k++) {
        for (i = 0; i < 3; i++) {
            lever[i] = macro->centroid[3 * k + i] - rcm[i];
            f[i] = Fk[3 * k + i];
            F[i] = F[i] + f[i];
        }
        tau[0] = tau[0] + (lever[1] * f[2] - lever[2] * f[1]);
        tau[1] = tau[1] + (lever[2] * f[0] - lever[0] * f[2]);
        tau[2] = tau[2] + (lever[0] * f[1] - lever[1] * f[0]);
    }
    return NGF_OK;
}

/* =========================================================================
 * singleEpoch  --  the composition Model.evaluate performs, for one epoch.
 *
 * cdScale and crScale are applied to the per-panel arrays at exactly the point
 * Model.evaluate applies them, so a scaled run and an unscaled run differ by
 * the factor and by nothing else.
 * ========================================================================= */
int ngfSingleEpoch(const ngfEpochInput *in, ngfEpochOutput *out)
{
    /* Per-panel scratch for at most this many plates without allocating; a
     * larger model writes straight into the caller's Fk, which is why Fk is
     * required when n exceeds the scratch. Nothing here calls malloc, so the
     * function is usable inside a right-hand side without a heap. */
    enum { NGF_SCRATCH_PLATES = 64 };
    double scratch[3 * NGF_SCRATCH_PLATES];
    double *panel;
    int status, n, i, e, useScratch;

    if (in == NULL || out == NULL || out->F == NULL || out->tau == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    status = ngfMacroOk(in->macro);
    if (status != NGF_OK) {
        return status;
    }
    if (in->con == NULL || in->rcm == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    n = in->macro->n;
    useScratch = (n <= NGF_SCRATCH_PLATES);
    if (!useScratch && out->Fk == NULL) {
        return NGF_NULL_ARGUMENT;
    }
    if (!ngfFinite(in->cdScale) || !ngfFinite(in->crScale)) {
        return NGF_NOT_FINITE;
    }

    for (e = 0; e < 3; e++) {
        panel = useScratch ? scratch : (out->Fk + 3 * n * e);
        switch (e) {
            case 0:
                status = ngfDragPanelForce(in->macro, in->con, in->rho, in->species,
                                           in->tGas, in->vrel, in->gsiModel,
                                           in->alpha, in->tw, in->cd, in->specular,
                                           in->shadowGas, panel,
                                           out->gamma, out->theta, out->sRatio,
                                           out->cpMix, out->ctauMix);
                if (status != NGF_OK) {
                    return status;
                }
                for (i = 0; i < 3 * n; i++) {
                    panel[i] = in->cdScale * panel[i];
                }
                break;
            case 1:
                status = ngfRadiationPanelForce(in->macro, in->con, in->src,
                                                in->shadowRad, panel, out->absorbed);
                if (status != NGF_OK) {
                    return status;
                }
                for (i = 0; i < 3 * n; i++) {
                    panel[i] = in->crScale * panel[i];
                }
                break;
            default:
                if (in->tPanel == NULL) {
                    for (i = 0; i < 3 * n; i++) {
                        panel[i] = 0.0;
                    }
                    if (out->emitted != NULL) {
                        for (i = 0; i < n; i++) {
                            out->emitted[i] = 0.0;
                        }
                    }
                } else {
                    status = ngfThermalRecoil(in->macro, in->con, in->tPanel,
                                              panel, out->emitted);
                    if (status != NGF_OK) {
                        return status;
                    }
                }
                break;
        }
        status = ngfAssemble(in->macro, panel, in->rcm, out->F + 3 * e, out->tau + 3 * e);
        if (status != NGF_OK) {
            return status;
        }
        if (useScratch && out->Fk != NULL) {
            for (i = 0; i < 3 * n; i++) {
                out->Fk[3 * n * e + i] = panel[i];
            }
        }
    }
    return NGF_OK;
}
