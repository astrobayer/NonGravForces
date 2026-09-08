/*
 * ngfKernelParity.c
 * -----------------
 * Runs mex/ngfKernel.c against the FROZEN kernel fixtures in tools/pyxcheck and
 * reports the worst disagreement per unit. No MATLAB, no MEX, no numpy: a C99
 * compiler and the CSV files, so the CI `audit` job on ubuntu-latest can gate
 * the C exactly as `tests/tNgfMirror.m` gates the MATLAB.
 *
 * WHY THIS IS THE RIGHT ORACLE
 *   tools/pyxcheck/ngfkernel/ is an independently derived mirror of the kernel,
 *   written from the equations and frozen when LP-5 went green (WP-33). The
 *   MATLAB kernel is measured against it by tests/tNgfMirror.m; this program
 *   measures the C against the SAME files with the SAME bounds and the same
 *   big/negligible policy. Two ports, one authority -- which is the reason
 *   HANDOVER gives for freezing the mirror before LP-6b was written.
 *
 *   It does NOT replace the direct C-vs-MATLAB parity suite (tests/tNgfMex.m):
 *   that one needs MATLAB and is what the 1e-14 claim of the level prompt rests
 *   on. This program is what can be run where MATLAB cannot.
 *
 * BUILD AND RUN
 *   cc -std=c99 -O2 -I mex -o ngfKernelParity mex/test/ngfKernelParity.c \
 *      mex/ngfKernel.c -lm
 *   ./ngfKernelParity tools/pyxcheck
 *
 * Exit status 0 when every bound holds, 1 otherwise; every measured number is
 * printed whether it passes or not (rule T3: a passing bound leaves its number).
 *
 * -------------------------------------------------------------------------
 * NonGravForces 0.6.1
 * Created 2026-09-08 by Claude Opus 5 (Anthropic) for Matthias Weigelt.
 * -------------------------------------------------------------------------
 */
#include "ngfKernel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PLATES 16
#define MAX_FIELDS 64
#define MAX_LINE   4096
#define REALMIN    2.2250738585072014e-308

/* Constants: the same numbers NonGravForces.const and
 * referenceSystem.iersConstants publish, and the same the mirror uses. */
static const double MOLAR_G[NGF_SPECIES] = {
    4.002602, 15.9994, 28.0134, 31.9988, 39.948, 1.00794, 14.0067, 15.9994
};
#define N_AVOGADRO 6.02214076e23

static ngfConstants makeConstants(void)
{
    ngfConstants c;
    int j;
    c.cLight  = 299792458.0;
    c.sigmaSB = 5.670374419e-8;
    c.kB      = 1.380649e-23;
    for (j = 0; j < NGF_SPECIES; j++) {
        c.speciesMass[j] = (MOLAR_G[j] * 1e-3) / N_AVOGADRO;
    }
    return c;
}

/* ------------------------------------------------------------------ bodies */
typedef struct {
    char     name[32];
    int      n;
    double   area[MAX_PLATES];
    double   normal[3 * MAX_PLATES];
    double   centroid[3 * MAX_PLATES];
    double   optVis[3 * MAX_PLATES];
    double   optIR[3 * MAX_PLATES];
    ngfMacro macro;
} bodyModel;

static bodyModel bodies[4];
static int       bodyCount = 0;

static bodyModel *findBody(const char *name)
{
    int i;
    for (i = 0; i < bodyCount; i++) {
        if (strcmp(bodies[i].name, name) == 0) {
            return &bodies[i];
        }
    }
    return NULL;
}

/* --------------------------------------------------------------- CSV reader
 * The fixtures are written by csv.writer with repr() floats and bare words; no
 * field is quoted and none contains a comma (auditPackage enforces the same for
 * the catalogue). A minimal splitter is therefore correct and is easier to
 * trust than a general parser. */
typedef struct {
    FILE  *fp;
    char   line[MAX_LINE];
    char  *field[MAX_FIELDS];
    int    nField;
    char  *header[MAX_FIELDS];
    char   headerBuf[MAX_LINE];
    int    nHeader;
} csvReader;

static int csvSplit(char *s, char **out, int maxOut)
{
    int n = 0;
    char *p = s;
    size_t len = strlen(s);
    while (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r')) {
        s[--len] = '\0';
    }
    out[n++] = p;
    while (*p != '\0' && n < maxOut) {
        if (*p == ',') {
            *p = '\0';
            out[n++] = p + 1;
        }
        p++;
    }
    return n;
}

static int csvOpen(csvReader *r, const char *path)
{
    r->fp = fopen(path, "r");
    if (r->fp == NULL) {
        fprintf(stderr, "cannot open %s\n", path);
        return 0;
    }
    if (fgets(r->headerBuf, MAX_LINE, r->fp) == NULL) {
        fprintf(stderr, "%s is empty\n", path);
        return 0;
    }
    r->nHeader = csvSplit(r->headerBuf, r->header, MAX_FIELDS);
    return 1;
}

static int csvNext(csvReader *r)
{
    if (fgets(r->line, MAX_LINE, r->fp) == NULL) {
        return 0;
    }
    if (r->line[0] == '\n' || r->line[0] == '\r' || r->line[0] == '\0') {
        return 0;
    }
    r->nField = csvSplit(r->line, r->field, MAX_FIELDS);
    return 1;
}

static int csvIndex(const csvReader *r, const char *name)
{
    int i;
    for (i = 0; i < r->nHeader; i++) {
        if (strcmp(r->header[i], name) == 0) {
            return i;
        }
    }
    fprintf(stderr, "column '%s' missing\n", name);
    exit(2);
}

static const char *csvStr(const csvReader *r, int i)
{
    return (i < r->nField) ? r->field[i] : "";
}

static double csvNum(const csvReader *r, int i)
{
    const char *s = csvStr(r, i);
    if (s[0] == '\0') {
        return NAN;                        /* an empty field: not applicable */
    }
    return strtod(s, NULL);
}

/* ------------------------------------------------------------- book-keeping */
typedef struct {
    const char *name;
    double      measured;
    double      bound;
    const char *note;
} boundRecord;

static boundRecord records[32];
static int         recordCount = 0;
static int         failures = 0;

static void record(const char *name, double measured, double bound, const char *note)
{
    records[recordCount].name = name;
    records[recordCount].measured = measured;
    records[recordCount].bound = bound;
    records[recordCount].note = note;
    recordCount++;
    if (!(measured == measured) || measured > bound) {
        failures++;
    }
}

static double dmax(double a, double b) { return (a > b) ? a : b; }

/* ------------------------------------------------------------- load bodies */
static void loadBodies(const char *dir)
{
    csvReader r;
    char path[512];
    int iBody, iArea, inx, iny, inz, icx, icy, icz, ia, is, id, ie, isi, idi;
    bodyModel *b;
    int k, j;

    snprintf(path, sizeof path, "%s/kernel_mirror_plates.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iBody = csvIndex(&r, "body");
    iArea = csvIndex(&r, "area_m2");
    inx = csvIndex(&r, "nx"); iny = csvIndex(&r, "ny"); inz = csvIndex(&r, "nz");
    icx = csvIndex(&r, "cx"); icy = csvIndex(&r, "cy"); icz = csvIndex(&r, "cz");
    ia = csvIndex(&r, "absVis"); is = csvIndex(&r, "specVis"); id = csvIndex(&r, "diffVis");
    ie = csvIndex(&r, "emisIR"); isi = csvIndex(&r, "specIR"); idi = csvIndex(&r, "diffIR");
    while (csvNext(&r)) {
        b = findBody(csvStr(&r, iBody));
        if (b == NULL) {
            b = &bodies[bodyCount++];
            snprintf(b->name, sizeof b->name, "%s", csvStr(&r, iBody));
            b->n = 0;
        }
        k = b->n;
        b->area[k] = csvNum(&r, iArea);
        b->normal[3 * k] = csvNum(&r, inx);
        b->normal[3 * k + 1] = csvNum(&r, iny);
        b->normal[3 * k + 2] = csvNum(&r, inz);
        b->centroid[3 * k] = csvNum(&r, icx);
        b->centroid[3 * k + 1] = csvNum(&r, icy);
        b->centroid[3 * k + 2] = csvNum(&r, icz);
        b->optVis[3 * k] = csvNum(&r, ia);
        b->optVis[3 * k + 1] = csvNum(&r, is);
        b->optVis[3 * k + 2] = csvNum(&r, id);
        b->optIR[3 * k] = csvNum(&r, ie);
        b->optIR[3 * k + 1] = csvNum(&r, isi);
        b->optIR[3 * k + 2] = csvNum(&r, idi);
        b->n = k + 1;
    }
    fclose(r.fp);
    for (j = 0; j < bodyCount; j++) {
        bodies[j].macro.n = bodies[j].n;
        bodies[j].macro.area = bodies[j].area;
        bodies[j].macro.normal = bodies[j].normal;
        bodies[j].macro.centroid = bodies[j].centroid;
        bodies[j].macro.optVis = bodies[j].optVis;
        bodies[j].macro.optIR = bodies[j].optIR;
        printf("  body %-12s %d plates\n", bodies[j].name, bodies[j].n);
    }
}

/* ------------------------------------------------------------------- gsi */
static void checkGsi(const char *dir)
{
    csvReader r;
    char path[512];
    int iModel, iS, iTheta, iAlpha, iTw, iTinf, iCd, iSpec, iCp, iCtau;
    double cp, ct, wantCp, wantCt, s, theta, alpha, tw, tinf, cd, spec;
    double bigCp = 0.0, bigCt = 0.0, negCp = 0.0, negCt = 0.0;
    int nBigCp = 0, nBigCt = 0, n = 0, status;
    ngfGsiModel model;

    snprintf(path, sizeof path, "%s/mirror_gsi.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iModel = csvIndex(&r, "model"); iS = csvIndex(&r, "s"); iTheta = csvIndex(&r, "theta");
    iAlpha = csvIndex(&r, "alpha"); iTw = csvIndex(&r, "tw_k"); iTinf = csvIndex(&r, "tinf_k");
    iCd = csvIndex(&r, "cd"); iSpec = csvIndex(&r, "specular");
    iCp = csvIndex(&r, "cp"); iCtau = csvIndex(&r, "ctau");
    while (csvNext(&r)) {
        n++;
        if (strcmp(csvStr(&r, iModel), "constantCd") == 0) {
            model = NGF_GSI_CONSTANT_CD;
        } else if (strcmp(csvStr(&r, iModel), "maxwell") == 0) {
            model = NGF_GSI_MAXWELL;
        } else {
            model = NGF_GSI_SENTMAN;
        }
        s = csvNum(&r, iS); theta = csvNum(&r, iTheta); alpha = csvNum(&r, iAlpha);
        tw = csvNum(&r, iTw); tinf = csvNum(&r, iTinf);
        cd = csvNum(&r, iCd); spec = csvNum(&r, iSpec);
        status = ngfGsiCoefficients(model, s, theta, alpha, tw, tinf, cd, spec, &cp, &ct);
        if (status != NGF_OK) {
            fprintf(stderr, "gsi row %d refused: %s\n", n, ngfStatusText(status));
            failures++;
            continue;
        }
        wantCp = csvNum(&r, iCp);
        wantCt = csvNum(&r, iCtau);
        if (fabs(wantCp) >= 1e-12) {
            bigCp = dmax(bigCp, fabs(cp - wantCp) / fabs(wantCp)); nBigCp++;
        } else {
            negCp = dmax(negCp, dmax(fabs(cp), fabs(wantCp)));
        }
        if (fabs(wantCt) >= 1e-12) {
            bigCt = dmax(bigCt, fabs(ct - wantCt) / fabs(wantCt)); nBigCt++;
        } else {
            negCt = dmax(negCt, dmax(fabs(ct), fabs(wantCt)));
        }
    }
    fclose(r.fp);
    printf("  gsiCoefficients: %d rows, %d/%d cp and ctau above 1e-12\n", n, nBigCp, nBigCt);
    record("gsiCoefficients cp vs the mirror (relative)", bigCp, 1e-12,
           "same closed form, C and MATLAB rounding");
    record("gsiCoefficients cp, negligible cases", negCp, 1e-11,
           "leeward plates: both sides negligible");
    record("gsiCoefficients ctau vs the mirror (relative)", bigCt, 1e-12,
           "same closed form, C and MATLAB rounding");
    record("gsiCoefficients ctau, negligible cases", negCt, 1e-11,
           "leeward plates: both sides negligible");
}

/* ------------------------------------------------------------------ drag */
static void checkDrag(const char *dir)
{
    csvReader r;
    char path[512];
    int iCase, iBody, iModel, iAlpha, iTw, iCd, iSpec, iT, iVx, iVy, iVz;
    int iD[NGF_SPECIES], iShadow, iPlate, iFx, iFy, iFz;
    int j, k, n, status, caseId, prevCase = -1, nCase = 0;
    char body[32];
    double species[NGF_SPECIES], shadow[MAX_PLATES], want[3 * MAX_PLATES];
    double got[3 * MAX_PLATES];
    double alpha, tw, cd, spec, tGas, vrel[3], rho;
    double worst = 0.0, worstNeg = 0.0, scale, per, plateMax;
    ngfGsiModel model = NGF_GSI_SENTMAN;
    bodyModel *b = NULL;
    ngfConstants con = makeConstants();
    char nameBuf[16];
    int have = 0;

    snprintf(path, sizeof path, "%s/mirror_drag.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iCase = csvIndex(&r, "caseId"); iBody = csvIndex(&r, "body"); iModel = csvIndex(&r, "model");
    iAlpha = csvIndex(&r, "alpha"); iTw = csvIndex(&r, "tw_k"); iCd = csvIndex(&r, "cd");
    iSpec = csvIndex(&r, "specular"); iT = csvIndex(&r, "t_k");
    iVx = csvIndex(&r, "vx"); iVy = csvIndex(&r, "vy"); iVz = csvIndex(&r, "vz");
    for (j = 0; j < NGF_SPECIES; j++) {
        snprintf(nameBuf, sizeof nameBuf, "d%d", j + 1);
        iD[j] = csvIndex(&r, nameBuf);
    }
    iShadow = csvIndex(&r, "shadow"); iPlate = csvIndex(&r, "plate");
    iFx = csvIndex(&r, "fx"); iFy = csvIndex(&r, "fy"); iFz = csvIndex(&r, "fz");

    /* Rows of one case are contiguous and ordered by plate; a case is evaluated
     * when the caseId changes and once more at end of file. */
    while (1) {
        have = csvNext(&r);
        caseId = have ? (int)csvNum(&r, iCase) : -1;
        if ((!have || caseId != prevCase) && prevCase != -1) {
            n = b->n;
            rho = 0.0;
            for (j = 0; j < NGF_SPECIES; j++) {
                rho += species[j] * con.speciesMass[j];
            }
            status = ngfDragPanelForce(&b->macro, &con, rho,
                                       (model == NGF_GSI_CONSTANT_CD) ? NULL : species,
                                       tGas, vrel, model, alpha, tw, cd, spec,
                                       shadow, got, NULL, NULL, NULL, NULL, NULL);
            if (status != NGF_OK) {
                fprintf(stderr, "drag case %d refused: %s\n", prevCase, ngfStatusText(status));
                failures++;
            } else {
                scale = 0.0;
                for (k = 0; k < 3 * n; k++) {
                    scale = dmax(scale, fabs(want[k]));
                }
                if (scale < REALMIN) {
                    scale = REALMIN;
                }
                for (k = 0; k < n; k++) {
                    plateMax = dmax(dmax(fabs(want[3 * k]), fabs(want[3 * k + 1])),
                                    fabs(want[3 * k + 2]));
                    if (plateMax >= 1e-12 * scale) {
                        for (j = 0; j < 3; j++) {
                            per = fabs(got[3 * k + j] - want[3 * k + j]) / scale;
                            worst = dmax(worst, per);
                        }
                    } else {
                        for (j = 0; j < 3; j++) {
                            worstNeg = dmax(worstNeg, fabs(got[3 * k + j]) / scale);
                        }
                    }
                }
            }
            nCase++;
        }
        if (!have) {
            break;
        }
        if (caseId != prevCase) {
            snprintf(body, sizeof body, "%s", csvStr(&r, iBody));
            b = findBody(body);
            if (b == NULL) {
                fprintf(stderr, "unknown body %s\n", body);
                exit(2);
            }
            if (strcmp(csvStr(&r, iModel), "constantCd") == 0) {
                model = NGF_GSI_CONSTANT_CD;
            } else if (strcmp(csvStr(&r, iModel), "maxwell") == 0) {
                model = NGF_GSI_MAXWELL;
            } else {
                model = NGF_GSI_SENTMAN;
            }
            alpha = csvNum(&r, iAlpha); tw = csvNum(&r, iTw);
            cd = csvNum(&r, iCd); spec = csvNum(&r, iSpec); tGas = csvNum(&r, iT);
            vrel[0] = csvNum(&r, iVx); vrel[1] = csvNum(&r, iVy); vrel[2] = csvNum(&r, iVz);
            for (j = 0; j < NGF_SPECIES; j++) {
                species[j] = csvNum(&r, iD[j]);
            }
            prevCase = caseId;
        }
        k = (int)csvNum(&r, iPlate) - 1;
        shadow[k] = csvNum(&r, iShadow);
        want[3 * k] = csvNum(&r, iFx);
        want[3 * k + 1] = csvNum(&r, iFy);
        want[3 * k + 2] = csvNum(&r, iFz);
    }
    fclose(r.fp);
    printf("  dragPanelForce: %d cases\n", nCase);
    record("dragPanelForce vs the mirror (relative to the case scale)", worst, 1e-12,
           "per-species mixture rule, C and MATLAB rounding");
    record("dragPanelForce, plates negligible in their case", worstNeg, 1e-11,
           "leeward plates: negligible on both sides");
}

/* ------------------------------------------------------------- radiation */
static void checkRadiation(const char *dir)
{
    csvReader r;
    char path[512];
    int iCase, iBody, iPlate, iFx, iFy, iFz, iAbs;
    int iDx[3], iDy[3], iDz[3], iP[3], iBand[3], iSh[3];
    int m, k, j, n, status, caseId, prevCase = -1, nCase = 0, have;
    char nameBuf[16], body[32];
    double dir3[9], pres[3], want[3 * MAX_PLATES], wantAbs[MAX_PLATES];
    double got[3 * MAX_PLATES], absorbed[MAX_PLATES], shadow[3 * MAX_PLATES];
    int band[3];
    double worst = 0.0, worstAbs = 0.0, scale, absScale;
    bodyModel *b = NULL;
    ngfConstants con = makeConstants();
    ngfSources src;

    snprintf(path, sizeof path, "%s/mirror_radiation.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iCase = csvIndex(&r, "caseId"); iBody = csvIndex(&r, "body");
    for (m = 0; m < 3; m++) {
        snprintf(nameBuf, sizeof nameBuf, "s%ddx", m + 1); iDx[m] = csvIndex(&r, nameBuf);
        snprintf(nameBuf, sizeof nameBuf, "s%ddy", m + 1); iDy[m] = csvIndex(&r, nameBuf);
        snprintf(nameBuf, sizeof nameBuf, "s%ddz", m + 1); iDz[m] = csvIndex(&r, nameBuf);
        snprintf(nameBuf, sizeof nameBuf, "s%dp", m + 1);  iP[m]  = csvIndex(&r, nameBuf);
        snprintf(nameBuf, sizeof nameBuf, "s%dband", m + 1); iBand[m] = csvIndex(&r, nameBuf);
        snprintf(nameBuf, sizeof nameBuf, "s%dshadow", m + 1); iSh[m] = csvIndex(&r, nameBuf);
    }
    iPlate = csvIndex(&r, "plate");
    iFx = csvIndex(&r, "fx"); iFy = csvIndex(&r, "fy"); iFz = csvIndex(&r, "fz");
    iAbs = csvIndex(&r, "absorbed");

    while (1) {
        have = csvNext(&r);
        caseId = have ? (int)csvNum(&r, iCase) : -1;
        if ((!have || caseId != prevCase) && prevCase != -1) {
            n = b->n;
            src.m = 3; src.dir = dir3; src.pres = pres; src.band = band;
            status = ngfRadiationPanelForce(&b->macro, &con, &src, shadow, got, absorbed);
            if (status != NGF_OK) {
                fprintf(stderr, "radiation case %d refused: %s\n", prevCase, ngfStatusText(status));
                failures++;
            } else {
                scale = REALMIN; absScale = REALMIN;
                for (k = 0; k < 3 * n; k++) {
                    scale = dmax(scale, fabs(want[k]));
                }
                for (k = 0; k < n; k++) {
                    absScale = dmax(absScale, fabs(wantAbs[k]));
                }
                for (k = 0; k < 3 * n; k++) {
                    worst = dmax(worst, fabs(got[k] - want[k]) / scale);
                }
                for (k = 0; k < n; k++) {
                    worstAbs = dmax(worstAbs, fabs(absorbed[k] - wantAbs[k]) / absScale);
                }
            }
            nCase++;
        }
        if (!have) {
            break;
        }
        if (caseId != prevCase) {
            snprintf(body, sizeof body, "%s", csvStr(&r, iBody));
            b = findBody(body);
            if (b == NULL) {
                fprintf(stderr, "unknown body %s\n", body);
                exit(2);
            }
            for (m = 0; m < 3; m++) {
                dir3[3 * m] = csvNum(&r, iDx[m]);
                dir3[3 * m + 1] = csvNum(&r, iDy[m]);
                dir3[3 * m + 2] = csvNum(&r, iDz[m]);
                pres[m] = csvNum(&r, iP[m]);
                band[m] = (strcmp(csvStr(&r, iBand[m]), "ir") == 0) ? NGF_BAND_IR : NGF_BAND_VIS;
            }
            prevCase = caseId;
        }
        k = (int)csvNum(&r, iPlate) - 1;
        for (m = 0; m < 3; m++) {
            shadow[b->n * m + k] = csvNum(&r, iSh[m]);
        }
        want[3 * k] = csvNum(&r, iFx);
        want[3 * k + 1] = csvNum(&r, iFy);
        want[3 * k + 2] = csvNum(&r, iFz);
        wantAbs[k] = csvNum(&r, iAbs);
        (void)j;
    }
    fclose(r.fp);
    printf("  radiationPanelForce: %d cases\n", nCase);
    record("radiationPanelForce vs the mirror (relative to the case scale)", worst, 1e-12,
           "flat-plate optical model, C and MATLAB rounding");
    record("radiationPanelForce absorbed flux vs the mirror", worstAbs, 1e-12,
           "a E gamma summed over sources");
}

/* --------------------------------------------------------------- thermal */
static void checkThermal(const char *dir)
{
    csvReader r;
    char path[512];
    int iCase, iBody, iPlate, iT, iFx, iFy, iFz;
    int k, n, status, caseId, prevCase = -1, nCase = 0, have;
    char body[32];
    double temp[MAX_PLATES], want[3 * MAX_PLATES], got[3 * MAX_PLATES];
    double worst = 0.0, scale;
    bodyModel *b = NULL;
    ngfConstants con = makeConstants();

    snprintf(path, sizeof path, "%s/mirror_thermal.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iCase = csvIndex(&r, "caseId"); iBody = csvIndex(&r, "body");
    iPlate = csvIndex(&r, "plate"); iT = csvIndex(&r, "temp_k");
    iFx = csvIndex(&r, "fx"); iFy = csvIndex(&r, "fy"); iFz = csvIndex(&r, "fz");
    while (1) {
        have = csvNext(&r);
        caseId = have ? (int)csvNum(&r, iCase) : -1;
        if ((!have || caseId != prevCase) && prevCase != -1) {
            n = b->n;
            status = ngfThermalRecoil(&b->macro, &con, temp, got, NULL);
            if (status != NGF_OK) {
                fprintf(stderr, "thermal case %d refused: %s\n", prevCase, ngfStatusText(status));
                failures++;
            } else {
                scale = REALMIN;
                for (k = 0; k < 3 * n; k++) {
                    scale = dmax(scale, fabs(want[k]));
                }
                for (k = 0; k < 3 * n; k++) {
                    worst = dmax(worst, fabs(got[k] - want[k]) / scale);
                }
            }
            nCase++;
        }
        if (!have) {
            break;
        }
        if (caseId != prevCase) {
            snprintf(body, sizeof body, "%s", csvStr(&r, iBody));
            b = findBody(body);
            if (b == NULL) {
                fprintf(stderr, "unknown body %s\n", body);
                exit(2);
            }
            prevCase = caseId;
        }
        k = (int)csvNum(&r, iPlate) - 1;
        temp[k] = csvNum(&r, iT);
        want[3 * k] = csvNum(&r, iFx);
        want[3 * k + 1] = csvNum(&r, iFy);
        want[3 * k + 2] = csvNum(&r, iFz);
    }
    fclose(r.fp);
    printf("  thermalRecoil: %d cases\n", nCase);
    record("thermalRecoil vs the mirror", worst, 1e-13,
           "Stefan-Boltzmann and one hemispherical factor: only rounding may differ");
}

/* -------------------------------------------------------------- assemble */
static void checkAssemble(const char *dir)
{
    csvReader r;
    char path[512];
    int iCase, iBody, iPlate, iIx, iIy, iIz, iRx, iRy, iRz;
    int iFx, iFy, iFz, iTx, iTy, iTz;
    int k, status, caseId, prevCase = -1, nCase = 0, have;
    char body[32];
    double fk[3 * MAX_PLATES], rcm[3], F[3], tau[3], wantF[3], wantT[3];
    double wF = 0.0, wT = 0.0, sF, sT;
    bodyModel *b = NULL;

    snprintf(path, sizeof path, "%s/mirror_assemble.csv", dir);
    if (!csvOpen(&r, path)) {
        exit(2);
    }
    iCase = csvIndex(&r, "caseId"); iBody = csvIndex(&r, "body"); iPlate = csvIndex(&r, "plate");
    iIx = csvIndex(&r, "inx"); iIy = csvIndex(&r, "iny"); iIz = csvIndex(&r, "inz");
    iRx = csvIndex(&r, "rcmx"); iRy = csvIndex(&r, "rcmy"); iRz = csvIndex(&r, "rcmz");
    iFx = csvIndex(&r, "Fx"); iFy = csvIndex(&r, "Fy"); iFz = csvIndex(&r, "Fz");
    iTx = csvIndex(&r, "Tx"); iTy = csvIndex(&r, "Ty"); iTz = csvIndex(&r, "Tz");
    while (1) {
        have = csvNext(&r);
        caseId = have ? (int)csvNum(&r, iCase) : -1;
        if ((!have || caseId != prevCase) && prevCase != -1) {
            status = ngfAssemble(&b->macro, fk, rcm, F, tau);
            if (status != NGF_OK) {
                fprintf(stderr, "assemble case %d refused: %s\n", prevCase, ngfStatusText(status));
                failures++;
            } else {
                sF = REALMIN; sT = REALMIN;
                for (k = 0; k < 3; k++) {
                    sF = dmax(sF, fabs(wantF[k]));
                    sT = dmax(sT, fabs(wantT[k]));
                }
                for (k = 0; k < 3; k++) {
                    wF = dmax(wF, fabs(F[k] - wantF[k]) / sF);
                    wT = dmax(wT, fabs(tau[k] - wantT[k]) / sT);
                }
            }
            nCase++;
        }
        if (!have) {
            break;
        }
        if (caseId != prevCase) {
            snprintf(body, sizeof body, "%s", csvStr(&r, iBody));
            b = findBody(body);
            if (b == NULL) {
                fprintf(stderr, "unknown body %s\n", body);
                exit(2);
            }
            rcm[0] = csvNum(&r, iRx); rcm[1] = csvNum(&r, iRy); rcm[2] = csvNum(&r, iRz);
            wantF[0] = csvNum(&r, iFx); wantF[1] = csvNum(&r, iFy); wantF[2] = csvNum(&r, iFz);
            wantT[0] = csvNum(&r, iTx); wantT[1] = csvNum(&r, iTy); wantT[2] = csvNum(&r, iTz);
            prevCase = caseId;
        }
        k = (int)csvNum(&r, iPlate) - 1;
        fk[3 * k] = csvNum(&r, iIx);
        fk[3 * k + 1] = csvNum(&r, iIy);
        fk[3 * k + 2] = csvNum(&r, iIz);
    }
    fclose(r.fp);
    printf("  assemble: %d cases\n", nCase);
    record("assemble force vs the mirror", wF, 1e-13, "sums of 6 to 8 terms in two orders");
    record("assemble torque vs the mirror", wT, 1e-12, "cross products then sums, two orders");
}

/* ----------------------------------------------------- composition check
 * ngfSingleEpoch must equal the three units called separately, exactly: it is a
 * composition, not a re-derivation. Anything but zero here means the composite
 * path took a different route through the same physics. */
static void checkComposite(void)
{
    ngfConstants con = makeConstants();
    bodyModel *b = findBody("graceFO_C");
    ngfEpochInput in;
    ngfSources src;
    double species[NGF_SPECIES] = {2.0e13, 8.0e13, 4.0e11, 1.0e10, 2.0e9, 5.0e10, 1.0e12, 1.0e12};
    double vrel[3] = {7590.0, 120.0, -45.0};
    double dir3[9] = {1.0, 0.0, 0.0,
                      -0.19518001458970663, 0.09759000729485331, -0.9759000729485331,
                      0.0984135662610246, -0.14762034939153687, -0.9841356626102459};
    double pres[3] = {4.5606e-6, 1.83e-7, 7.4e-7};
    int band[3] = {NGF_BAND_VIS, NGF_BAND_VIS, NGF_BAND_IR};
    double tPanel[MAX_PLATES], shadowGas[MAX_PLATES], shadowRad[3 * MAX_PLATES];
    double Fk[3 * 3 * MAX_PLATES], F[9], tau[9];
    double refPanel[3 * MAX_PLATES], refF[3], refTau[3];
    double rcm[3], rho = 0.0, worst = 0.0, worstT = 0.0;
    double gammaOut[MAX_PLATES], thetaOut[MAX_PLATES], sRatioOut[NGF_SPECIES];
    double gammaRef[MAX_PLATES], thetaRef[MAX_PLATES], sRatioRef[NGF_SPECIES];
    ngfEpochOutput out;
    int k, m, j, n, status, e;

    if (b == NULL) {
        fprintf(stderr, "graceFO_C missing\n");
        exit(2);
    }
    n = b->n;
    rcm[0] = -0.000249; rcm[1] = -0.000403; rcm[2] = 0.001266;
    for (k = 0; k < n; k++) {
        tPanel[k] = 250.0 + 7.0 * k;
        shadowGas[k] = 0.6 + 0.05 * k;
        for (m = 0; m < 3; m++) {
            shadowRad[n * m + k] = 0.25 + 0.25 * ((k + m) % 3);
        }
    }
    for (j = 0; j < NGF_SPECIES; j++) {
        rho += species[j] * con.speciesMass[j];
    }
    src.m = 3; src.dir = dir3; src.pres = pres; src.band = band;

    in.macro = &b->macro; in.con = &con;
    in.rho = rho; in.species = species; in.tGas = 950.0; in.vrel = vrel;
    in.gsiModel = NGF_GSI_SENTMAN; in.alpha = 0.85; in.tw = 300.0;
    in.cd = 2.2; in.specular = 0.0; in.shadowGas = shadowGas; in.cdScale = 1.03;
    in.src = &src; in.shadowRad = shadowRad; in.crScale = 0.97;
    in.tPanel = tPanel; in.rcm = rcm;

    memset(&out, 0, sizeof out);
    out.F = F; out.tau = tau; out.Fk = Fk;
    out.gamma = gammaOut; out.theta = thetaOut; out.sRatio = sRatioOut;
    status = ngfSingleEpoch(&in, &out);
    if (status != NGF_OK) {
        fprintf(stderr, "singleEpoch refused: %s\n", ngfStatusText(status));
        failures++;
        return;
    }
    for (e = 0; e < 3; e++) {
        if (e == 0) {
            ngfDragPanelForce(&b->macro, &con, rho, species, 950.0, vrel,
                              NGF_GSI_SENTMAN, 0.85, 300.0, 2.2, 0.0, shadowGas,
                              refPanel, gammaRef, thetaRef, sRatioRef, NULL, NULL);
            for (k = 0; k < n; k++) {
                worst = dmax(worst, fabs(gammaOut[k] - gammaRef[k]));
                worst = dmax(worst, fabs(thetaOut[k] - thetaRef[k]));
            }
            for (k = 0; k < NGF_SPECIES; k++) {
                worst = dmax(worst, fabs(sRatioOut[k] - sRatioRef[k]));
            }
            for (k = 0; k < 3 * n; k++) { refPanel[k] = 1.03 * refPanel[k]; }
        } else if (e == 1) {
            ngfRadiationPanelForce(&b->macro, &con, &src, shadowRad, refPanel, NULL);
            for (k = 0; k < 3 * n; k++) { refPanel[k] = 0.97 * refPanel[k]; }
        } else {
            ngfThermalRecoil(&b->macro, &con, tPanel, refPanel, NULL);
        }
        ngfAssemble(&b->macro, refPanel, rcm, refF, refTau);
        for (k = 0; k < 3 * n; k++) {
            worst = dmax(worst, fabs(Fk[3 * n * e + k] - refPanel[k]));
        }
        for (k = 0; k < 3; k++) {
            worst = dmax(worst, fabs(F[3 * e + k] - refF[k]));
            worstT = dmax(worstT, fabs(tau[3 * e + k] - refTau[k]));
        }
    }
    printf("  singleEpoch: composition of 3 effects over %d plates\n", n);
    record("singleEpoch equals the units called separately (force)", worst, 0.0,
           "a composition must be bitwise, not merely close");
    record("singleEpoch equals the units called separately (torque)", worstT, 0.0,
           "a composition must be bitwise, not merely close");
}

/* ------------------------------------------------------ refusal behaviour
 * A refusal that never fires is a refusal that is not there (rule T6): each
 * guard is provoked once and paired with the case that must succeed. */
static void checkRefusals(void)
{
    ngfConstants con = makeConstants();
    bodyModel *b = findBody("graceFO_C");
    double cp, ct, F[3 * MAX_PLATES];
    double vrel[3] = {7600.0, 0.0, 0.0}, zero[3] = {0.0, 0.0, 0.0};
    double species[NGF_SPECIES] = {2.0e13, 8.0e13, 4.0e11, 1.0e10, 2.0e9, 5.0e10, 1.0e12, 1.0e12};
    int bad = 0;

    if (ngfGsiCoefficients(NGF_GSI_SENTMAN, 8.0, 0.3, 0.85, 300.0, 1000.0, 0.0, 0.0, &cp, &ct) != NGF_OK) { bad++; }
    if (ngfGsiCoefficients(NGF_GSI_SENTMAN, -1.0, 0.3, 0.85, 300.0, 1000.0, 0.0, 0.0, &cp, &ct) != NGF_BAD_RANGE) { bad++; }
    if (ngfGsiCoefficients(NGF_GSI_SENTMAN, 8.0, 4.0, 0.85, 300.0, 1000.0, 0.0, 0.0, &cp, &ct) != NGF_BAD_RANGE) { bad++; }
    if (ngfGsiCoefficients(NGF_GSI_SENTMAN, 8.0, NAN, 0.85, 300.0, 1000.0, 0.0, 0.0, &cp, &ct) != NGF_NOT_FINITE) { bad++; }
    if (ngfGsiCoefficients((ngfGsiModel)7, 8.0, 0.3, 0.85, 300.0, 1000.0, 0.0, 0.0, &cp, &ct) != NGF_BAD_MODEL) { bad++; }
    if (ngfGsiCoefficients(NGF_GSI_CONSTANT_CD, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, &cp, &ct) != NGF_BAD_RANGE) { bad++; }
    if (ngfGsiCoefficients(NGF_GSI_SENTMAN, 8.0, 0.3, 0.85, 300.0, 1000.0, 0.0, 0.0, NULL, &ct) != NGF_NULL_ARGUMENT) { bad++; }

    if (ngfDragPanelForce(&b->macro, &con, 1e-13, species, 950.0, vrel,
                          NGF_GSI_SENTMAN, 0.85, 300.0, 2.2, 0.0, NULL, F,
                          NULL, NULL, NULL, NULL, NULL) != NGF_OK) { bad++; }
    if (ngfDragPanelForce(&b->macro, &con, 1e-13, species, 950.0, zero,
                          NGF_GSI_SENTMAN, 0.85, 300.0, 2.2, 0.0, NULL, F,
                          NULL, NULL, NULL, NULL, NULL) != NGF_BAD_RANGE) { bad++; }
    if (ngfDragPanelForce(&b->macro, &con, -1.0, species, 950.0, vrel,
                          NGF_GSI_SENTMAN, 0.85, 300.0, 2.2, 0.0, NULL, F,
                          NULL, NULL, NULL, NULL, NULL) != NGF_BAD_RANGE) { bad++; }
    if (ngfThermalRecoil(&b->macro, &con, NULL, F, NULL) != NGF_NULL_ARGUMENT) { bad++; }

    printf("  refusals: 11 guards provoked\n");
    record("guards that behaved as specified (count of misses)", (double)bad, 0.0,
           "each refusal paired with the call that must succeed");
}

int main(int argc, char **argv)
{
    const char *dir = (argc > 1) ? argv[1] : "tools/pyxcheck";
    int i;

    printf("ngfKernelParity: NonGravForces C kernel %s against the frozen mirror in %s\n",
           ngfKernelVersion(), dir);
    loadBodies(dir);
    checkGsi(dir);
    checkDrag(dir);
    checkRadiation(dir);
    checkThermal(dir);
    checkAssemble(dir);
    checkComposite();
    checkRefusals();

    printf("\n%-58s %12s %12s\n", "measurement", "measured", "bound");
    for (i = 0; i < recordCount; i++) {
        printf("%-58s %12.3e %12.3e  %s\n", records[i].name, records[i].measured,
               records[i].bound, (records[i].measured <= records[i].bound) ? "ok" : "EXCEEDED");
    }
    if (failures > 0) {
        printf("\nngfKernelParity: %d bound(s) exceeded\n", failures);
        return 1;
    }
    printf("\nngfKernelParity: clean (%d bounds)\n", recordCount);
    return 0;
}
