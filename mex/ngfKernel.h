/*
 * ngfKernel.h
 * -----------
 * The single-epoch NonGravForces kernel in C (HANDOVER LP-6b, WP-17b).
 *
 * WHAT THIS IS
 *   A transcription of the FROZEN kernel layer of NonGravForces -- the
 *   gas-surface coefficients, the per-panel gas force, the per-panel radiation
 *   force, the thermal recoil and the force/torque reduction -- for ONE epoch.
 *   It exists so that an integrator written in C, or a MATLAB right-hand side
 *   called once per step, can evaluate the non-gravitational force without
 *   MATLAB in the inner loop.
 *
 * WHAT IT IS NOT
 *   It is not a second authority. MATLAB is normative (CLAUDE.md sect. 2): the
 *   expressions below are written in the same order as
  *   the +NonGravForces/+kernel units so that the two agree to rounding, and the
 *   agreement is asserted, not assumed -- see mex/test/ngfKernelParity.c, which
 *   runs the same frozen fixtures the MATLAB kernel is checked against
 *   (tools/pyxcheck/mirror_*.csv, WP-33).
 *
 *   It computes NO environment: no density, no radiation source list, no
 *   attitude, no self-shadowing, no panel temperature. Those are env/geom and
 *   they stay in MATLAB. Everything here is body frame, SI, radians.
 *
 * NO ALLOCATION, NO GLOBALS, NO PRINTING
 *   Every function writes only through its output pointers and is re-entrant,
 *   so a parallel integration may share one ngfMacro across threads without a
 *   lock. Errors are returned as codes; nothing is printed and nothing exits.
 *
 * LAYOUT
 *   Vectors of length 3 and 3-by-N arrays are column-major, exactly as MATLAB
 *   stores them, so a MEX gateway hands mxGetPr() straight through:
 *   normal[3*k + i] is component i of plate k.
 *
 * -------------------------------------------------------------------------
 * NonGravForces 0.6.1
 * Created 2026-09-08 by Claude Opus 5 (Anthropic) for Matthias Weigelt.
 * -------------------------------------------------------------------------
 */
#ifndef NGF_KERNEL_H
#define NGF_KERNEL_H

#ifdef __cplusplus
extern "C" {
#endif

#define NGF_SPECIES 8            /* He O N2 O2 Ar H N O_anomalous            */

/* Return codes. Anything but NGF_OK leaves every output untouched, so a failed
 * call can never be mistaken for a computed force. */
#define NGF_OK                 0
#define NGF_NULL_ARGUMENT      1
#define NGF_BAD_PLATE_COUNT    2   /* n < 1                                   */
#define NGF_BAD_MODEL          3   /* gsi model outside the enum              */
#define NGF_BAD_RANGE          4   /* s <= 0, theta outside [0,pi], alpha
                                    * outside [0,1], Tw/Tinf <= 0, |v| = 0,
                                    * rho < 0, P < 0, T < 0, shadow outside
                                    * [0,1], Cd <= 0, specular outside [0,1]  */
#define NGF_NOT_FINITE         5   /* NaN or Inf reached an input             */
#define NGF_BAD_SOURCE_COUNT   6   /* m < 0, or a band outside {0,1}          */

/* Gas-surface interaction model. cll and schamberg are refused in MATLAB
 * (HANDOVER O-17) and have no code here either. */
typedef enum {
    NGF_GSI_SENTMAN     = 0,
    NGF_GSI_MAXWELL     = 1,
    NGF_GSI_CONSTANT_CD = 2
} ngfGsiModel;

/* Radiation band of a source. */
#define NGF_BAND_VIS 0
#define NGF_BAND_IR  1

/*
 * The macro model. Pointers are BORROWED: nothing here copies or frees them.
 *
 *   n         plate count, >= 1
 *   area      [n]    m^2, > 0
 *   normal    [3*n]  outward UNIT normals, body frame. MacroModel renormalises
 *                    them (HANDOVER O-2); nothing here renormalises again.
 *   centroid  [3*n]  m, body frame; only ngfAssemble reads them
 *   optVis    [3*n]  per plate: absorbed, specular, diffuse (sum 1)
 *   optIR     [3*n]  per plate: emissivity, specular, diffuse (sum 1)
 */
typedef struct {
    int           n;
    const double *area;
    const double *normal;
    const double *centroid;
    const double *optVis;
    const double *optIR;
} ngfMacro;

/*
 * A radiation source list for one epoch (env.radiationSources, one page).
 *
 *   m      source count, >= 0 (0 is legal and gives zero force)
 *   dir    [3*m]  UNIT vectors satellite -> source, body frame
 *   pres   [m]    Pa, >= 0; the eclipse factor is already inside the Sun's
 *   band   [m]    NGF_BAND_VIS or NGF_BAND_IR
 */
typedef struct {
    int           m;
    const double *dir;
    const double *pres;
    const int    *band;
} ngfSources;

/*
 * Physical constants the caller supplies rather than this file inventing them:
 * one authority per fact (rule R1). NonGravForces.const and
 * referenceSystem.iersConstants own these numbers in MATLAB.
 *
 *   cLight       m/s      referenceSystem.iersConstants().speedOfLight
 *   sigmaSB      W/m^2/K^4 NonGravForces.const().sigmaSB
 *   kB           J/K      NonGravForces.const().kB
 *   speciesMass  [8] kg   PER MOLECULE, i.e. const().speciesMass / const().NA
 */
typedef struct {
    double cLight;
    double sigmaSB;
    double kB;
    double speciesMass[NGF_SPECIES];
} ngfConstants;

/*
 * The normal and tangential free-molecular coefficients of one flat plate, in
 * the (1/2) rho U^2 convention. Port of +kernel/gsiCoefficients.m.
 *
 *   model     which coefficient model
 *   s         speed ratio U / sqrt(2 k T / m), > 0 (ignored by constantCd)
 *   theta     angle between the outward normal and the relative velocity,
 *             [0, pi] rad
 *   alpha     energy accommodation coefficient, [0, 1] (ignored by constantCd)
 *   tw        wall temperature, K, > 0 (ignored by constantCd)
 *   tinf      gas kinetic temperature, K, > 0 (ignored by constantCd)
 *   cd        drag coefficient, > 0; read only by constantCd
 *   specular  specular fraction, [0, 1]; read only by maxwell
 *   cp, ctau  outputs, may not be NULL
 */
int ngfGsiCoefficients(ngfGsiModel model, double s, double theta, double alpha,
                       double tw, double tinf, double cd, double specular,
                       double *cp, double *ctau);

/*
 * Per-panel gas force for one epoch. Port of +kernel/dragPanelForce.m.
 *
 *   macro     panel model
 *   con       constants (kB and speciesMass are read; cLight and sigmaSB are not)
 *   rho       kg/m^3 total mass density; read by constantCd only
 *   species   [8] m^-3 number densities in const.speciesName order; may be
 *             NULL for constantCd
 *   tGas      K, > 0; ignored by constantCd
 *   vrel      [3] m/s satellite velocity relative to the gas, body frame
 *   model, alpha, tw, cd, specular   as in ngfGsiCoefficients
 *   shadow    [n] self-shadow fraction per plate in [0,1], or NULL for all 1
 *   F         [3*n] output force per plate, N, body frame
 *   gamma     [n]  output n.v per plate, or NULL
 *   theta     [n]  output incidence angle per plate, rad, or NULL
 *   sRatio    [8]  output speed ratio per species, or NULL (NaN for constantCd)
 *   cpMix     [n]  output mixture normal coefficient per unit rho, or NULL
 *   ctauMix   [n]  output mixture tangential coefficient per unit rho, or NULL
 */
int ngfDragPanelForce(const ngfMacro *macro, const ngfConstants *con,
                      double rho, const double *species, double tGas,
                      const double *vrel, ngfGsiModel model, double alpha,
                      double tw, double cd, double specular,
                      const double *shadow, double *F,
                      double *gamma, double *theta, double *sRatio,
                      double *cpMix, double *ctauMix);

/*
 * Per-panel radiation force for one epoch, summed over the sources. Port of
 * +kernel/radiationPanelForce.m.
 *
 *   shadow    [n*m] fraction per plate and source, plate-major
 *             (shadow[n*j + k] for plate k and source j), or NULL for all 1.
 *             That is the k-fastest layout of MATLAB's Shadow(:, 1, m).
 *   F         [3*n] output, N
 *   absorbed  [n]   output absorbed flux, W/m^2, or NULL
 */
int ngfRadiationPanelForce(const ngfMacro *macro, const ngfConstants *con,
                           const ngfSources *src, const double *shadow,
                           double *F, double *absorbed);

/*
 * Per-panel Lambertian thermal recoil. Port of +kernel/thermalRecoil.m.
 *
 *   tPanel    [n] K, >= 0
 *   F         [3*n] output, N
 *   emitted   [n] output eps sigma T^4, W/m^2, or NULL
 */
int ngfThermalRecoil(const ngfMacro *macro, const ngfConstants *con,
                     const double *tPanel, double *F, double *emitted);

/*
 * Reduce per-panel forces to the total force and the torque about rcm. Port of
 * +kernel/assemble.m. Sums in plate order, never through a BLAS product, so K
 * single calls and one batched call agree bitwise (HANDOVER pitfall list).
 *
 *   Fk    [3*n] per-panel forces
 *   rcm   [3]   reference point of the torque, body frame
 *   F     [3]   output total force, may not be NULL
 *   tau   [3]   output torque, may not be NULL
 */
int ngfAssemble(const ngfMacro *macro, const double *Fk, const double *rcm,
                double *F, double *tau);

/*
 * Everything a single-epoch right-hand side needs, in one call.
 *
 * The three effects stay SEPARATE all the way out, because Model.accel divides
 * each by the mass and Model.budget reports them one by one; summing them here
 * would force the caller to re-split what it needs.
 *
 *   cdScale, crScale   Options.drag.cdScale and Options.srp.crScale, applied to
 *                      the gas and radiation per-panel forces exactly where
 *                      Model.evaluate applies them
 *   tPanel             [n] K panel temperatures, or NULL for no thermal recoil
 *                      (Options.thermal.model == "none")
 */
typedef struct {
    const ngfMacro     *macro;
    const ngfConstants *con;
    /* gas */
    double              rho;
    const double       *species;      /* [8] or NULL for constantCd          */
    double              tGas;
    const double       *vrel;         /* [3]                                 */
    ngfGsiModel         gsiModel;
    double              alpha;
    double              tw;
    double              cd;
    double              specular;
    const double       *shadowGas;    /* [n] or NULL                         */
    double              cdScale;
    /* radiation */
    const ngfSources   *src;          /* may be NULL: no radiation           */
    const double       *shadowRad;    /* [n*m] or NULL                       */
    double              crScale;
    /* thermal */
    const double       *tPanel;       /* [n] or NULL: no thermal recoil      */
    /* reduction */
    const double       *rcm;          /* [3]                                 */
} ngfEpochInput;

/*
 * Where ngfSingleEpoch writes. F and tau are required; every other member may
 * be NULL, and a NULL member costs nothing -- the diagnostics are computed by
 * the same pass that produces the force, never by a second evaluation.
 *
 *   F         [9]    columns 0,1,2 = drag, radiation, thermal; body frame, N
 *   tau       [9]    same order, about rcm, N m
 *   Fk        [9*n]  per-panel forces, page 0 drag, 1 radiation, 2 thermal
 *   gamma     [n]    n.v per plate for the gas stream
 *   theta     [n]    incidence angle of the gas stream, rad
 *   sRatio    [8]    speed ratio per species (NaN for constantCd)
 *   cpMix     [n]    mixture normal coefficient per unit rho
 *   ctauMix   [n]    mixture tangential coefficient per unit rho
 *   absorbed  [n]    absorbed radiative flux, W/m^2 (before crScale)
 *   emitted   [n]    emitted thermal flux eps sigma T^4, W/m^2
 */
typedef struct {
    double *F;
    double *tau;
    double *Fk;
    double *gamma;
    double *theta;
    double *sRatio;
    double *cpMix;
    double *ctauMix;
    double *absorbed;
    double *emitted;
} ngfEpochOutput;

int ngfSingleEpoch(const ngfEpochInput *in, ngfEpochOutput *out);

/* A human-readable name for a return code, for a gateway's error message. */
const char *ngfStatusText(int status);

/* The version this file was transcribed from, for provenance. */
const char *ngfKernelVersion(void);

#ifdef __cplusplus
}
#endif

#endif /* NGF_KERNEL_H */
