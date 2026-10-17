/*
 * ngfRhs.h
 * --------
 * The whole non-gravitational acceleration of one satellite state, in C, for an
 * orbit integrator (HANDOVER LP-C, WP-C9).
 *
 * WHAT THIS IS
 *   The composition NonGravForces.Model.accel performs -- attitude, relative
 *   velocity (co-rotation, HWM14), geodetic position, density (exponential,
 *   NRLMSISE-00, NRLMSIS 2.1), accommodation, the radiation source list (Sun with
 *   its shadow, Earth albedo and longwave from Knocke, a point source or CERES),
 *   self-shadowing of the gas stream and of every source (cube-map table), panel
 *   temperatures (steady state), the per-panel forces and their reduction -- for
 *   an n-by-3 block of GCRS states at one time, published as a gj_native_entry
 *   (odeSuite csrc/gj_native.h, ABI 1) so gjToolbox's native stepper can call it.
 *
 *   EVERYTHING THAT DEPENDS ON TIME ALONE IS A TABLE built once by
 *   NonGravForces.app.DynamicsSource.prepare in MATLAB, from the same units
 *   Model.accel calls: the ITRS->GCRS rotation and its rate, Sun and Moon, the
 *   attitude nodes, the space-weather drivers, the calendar of the UTC day, the
 *   mass, the CERES months and the self-shadow table. Everything that depends on
 *   the state is computed here. MATLAB is normative (CLAUDE.md sect. 2): every
 *   port below names the unit it follows and keeps its association order where
 *   that decides the last digits.
 *
 * NO ALLOCATION AFTER ngfRhsCreate, NO PRINTING, NO EXIT. The context owns its
 * tables and its scratch; ngfRhsDestroy frees them. The evaluation is not
 * re-entrant (it writes the context's scratch, and the Fortran models it calls
 * keep their own state), which is what the gj_native.h contract permits: the
 * stepper calls one entry from one thread.
 *
 * Created 2026-09-16 by Claude Opus 5 (Anthropic).
 */
#ifndef NGF_RHS_H
#define NGF_RHS_H

#include "ngfKernel.h"

#ifdef __cplusplus
extern "C" {
#endif

/* option enumerations, each the position of the value in Options' member list */
enum { NGF_DENS_EXPONENTIAL = 0, NGF_DENS_MSIS00 = 1, NGF_DENS_MSIS21 = 2 };
enum { NGF_WIND_NONE = 0, NGF_WIND_COROTATION = 1, NGF_WIND_HWM14 = 2 };
/* The order is Options.validateValue's `sets.drag_accommodation`, which app.DynamicsSource
 * indexes into. A value inserted in MATLAB and not here evaluates a DIFFERENT model on the
 * native path with nothing to say so, which is why the two move in one commit. */
enum { NGF_ACC_FIXED = 0, NGF_ACC_PILINSKI2010 = 1, NGF_ACC_PILINSKI2011 = 2,
       NGF_ACC_WALKER2014 = 3, NGF_ACC_SESAM = 4 };
enum { NGF_ECL_CYLINDRICAL = 0, NGF_ECL_CONICAL = 1, NGF_ECL_SOLAARS = 2 };
enum { NGF_ERP_NONE = 0, NGF_ERP_POINT = 1, NGF_ERP_KNOCKE = 2, NGF_ERP_CERES = 3 };
enum { NGF_THERMAL_NONE = 0, NGF_THERMAL_STEADY = 1 };
enum { NGF_SHADOW_NONE = 0, NGF_SHADOW_LUT = 1 };
enum { NGF_ATT_PROVIDER = 0, NGF_ATT_NOMINAL = 1 };

/* status codes beyond ngfKernel.h's */
#define NGF_DOMAIN        20   /* t outside the tables, or a state the models do not serve */
#define NGF_ALLOC         21
#define NGF_MODEL_FAIL    22   /* a density or wind model refused (its gateway checks failed) */

typedef struct ngfRhsContext ngfRhsContext;

/*
 * Everything the context copies at creation. Arrays are column-major as MATLAB
 * stores them; every pointer is read during ngfRhsCreate only and may be freed
 * by the caller afterwards.
 */
typedef struct {
    /* constants (NonGravForces.const, referenceSystem.iersConstants, WGS84) */
    double cLight, sigmaSB, kB, au, Rsun, Rmoon, aEarth, bEarth, fEarth, e2Earth, tsi, amu, eV, torr;
    double speciesMass[NGF_SPECIES];                  /* kg per molecule */
    /* macro model */
    int n;
    const double *area, *normal, *centroid, *optVis, *optIR, *internalShare;   /* internalShare may be NULL */
    double rcm[3];
    /* options */
    int gsi; double cd, specular, twK, cdScale, crScale;
    int accommodation; double alphaFixed, surfaceMassAmu;
    int density; double expRho0, expH0, expH;
    int wind;
    int eclipse, lunar;
    int erp, nRings, nSectors;
    int thermal; double internalPower;
    int shadow, lutResolution; const double *lut;     /* n x 6(res+1)^2 */
    int nChi, nArc; const double *xChi, *wChi, *xArc, *wArc;   /* Gauss-Legendre on [-1, 1] */
    /* geom.flowShadow Method "table" (Options geom.flowShadow): trilinear over (v_y/|v_x|, v_z/|v_x|, s_O),
       one table per cone (+x, -x), P = 24 * ftNShaded rows (plate, species, ratio) by ftN^2 ftNs nodes */
    int flowTable, ftN, ftNs, ftNShaded; double ftHalfWidthTan, ftT0, ftH, ftS0, ftHs;
    const double *ftShaded, *ftPlus, *ftMinus;       /* 1-based plate numbers; P x nodes, column-major */
    /* tables over [t0, t0 + (nT-1) dt], seconds from the run epoch (TAI) */
    double t0, dt; int nT;
    const double *R, *Rdot;                           /* 9 x nT, ITRS -> GCRS, column-major 3x3 */
    const double *sun, *moon, *eclipticNormal;        /* 3 x nT (moon, eclipticNormal may be NULL) */
    const double *mass;                               /* nT */
    int attitude, alongVelocity, senseTranspose;
    int nq; const double *tq, *q;                     /* nq, 4 x nq (sign-continued, normalised) */
    /* 3-hour UTC bins: true edges tb[i] = tb0 + i*10800 for the calendar, aligned edges for the drivers */
    int nb; double tb0; const double *alignedEdges;   /* nb + 1 */
    const double *f107, *f107a, *ap7, *ap3;           /* nb, nb, 7 x nb, nb */
    const double *doy, *yy, *utcSecAtBin, *utcDay;    /* nb each: day of year, year mod 100, UTC seconds and UTC MJD day at the bin start */
    /* CERES */
    int nLat, nLon, nMonth; const double *latC, *lonC, *tMonth;   /* tMonth: the month nodes as MJD UTC (io.importCERES .mjd) */
    const double *rsut, *rlut, *rsdt;                 /* nLat x nLon x nMonth */
    double ceresHeight;
    /* model data folders (msis21: with its trailing separator; hwm14: without) */
    const char *msis21Dir, *hwm14Dir;
    /* the widest n-by-m block the entry accepts */
    int maxColumns;
} ngfRhsSetup;

int  ngfRhsCreate(const ngfRhsSetup *setup, ngfRhsContext **ctx);
void ngfRhsDestroy(ngfRhsContext *ctx);

/*
 * The acceleration of m GCRS states at time t (seconds from the run epoch, TAI):
 * R, V and A are 3 x m, column-major. Returns NGF_OK or a status; A is written
 * only on success.
 */
int ngfRhsAccel(ngfRhsContext *ctx, double t, const double *R, const double *V, double *A, int m);

/*
 * The same for one state, with the intermediate quantities Model.accel reports,
 * for parity diagnostics: diag must hold NGF_RHS_DIAG doubles, laid out as
 *   [0..2] a_drag body, [3..5] a_rad body, [6..8] a_thermal body, [9] rho, [10..17] species,
 *   [18] T gas, [19..21] vrel body, [22..24] lat lon h, [25] nu, [26] mass, [27..35] R body->GCRS.
 */
#define NGF_RHS_DIAG 36
int ngfRhsDiagnose(ngfRhsContext *ctx, double t, const double *r, const double *v, double *diag);

/* The published entry (a gj_native_entry *), valid until ngfRhsDestroy. */
void *ngfRhsEntry(ngfRhsContext *ctx);

/* An optional second entry whose acceleration is added first (gravity). May be NULL. */
void ngfRhsSetBase(ngfRhsContext *ctx, void *baseEntry);

/* The status of the last failed evaluation, in words (which model or table refused). */
const char *ngfRhsLastError(const ngfRhsContext *ctx);

#ifdef __cplusplus
}
#endif

#endif /* NGF_RHS_H */
