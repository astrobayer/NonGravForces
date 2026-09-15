"""tpmc_check.py -- independent test-particle Monte Carlo oracle for self-shadowing
of a flat plate in free-molecular flow (pure standard-library Python).

Written from the physics statement only; the MATLAB implementation was not opened.

PHYSICS
-------
A flat target plate (outward unit normal n) sits in a rarefied gas. Relative to
the body the gas is a drifting Maxwellian:

  bulk velocity     W = -U * vhat, vhat = unit vector of the satellite velocity
                    relative to the gas. Molecules stream along -vhat; the
                    upstream direction (where most of them come from) is +vhat.
  thermal velocity  each Cartesian component ~ Normal(0, sigma^2),
                    sigma = U / (s * sqrt(2)),  s = U / sqrt(2 k T / m).

U = 1 throughout: every coefficient below is dimensionless, so U cancels.

A molecule with velocity c strikes the plate only if n.c < 0. It arrived along a
straight line, i.e. from the direction -c/|c| as seen from its impact point p. It
is OCCLUDED if the ray p - lambda*c, lambda > 0, meets an occluder rectangle.

With H = 1 if n.c < 0 else 0, and O = 1 if the molecule is NOT occluded else 0,
expectations taken over the normalised Maxwellian pdf f(c):

  number-flux coefficient      CN = E[ (-n.c) H ] / U
  normal-momentum coefficient  CP = 2 E[ (n.c)^2 H ] / U^2
  tangential-momentum coeff.   CT = 2 E[ (-n.c) (c.t) H ] / U^2

t is the unit vector of the tangential part of the stream direction g = -vhat:
t = (g - (g.n) n) / |g - (g.n) n|. If |g - (g.n) n| < 1e-12 (flow along the
normal), t is DEFINED as (1,0,0) projected onto the plate plane and normalised,
which for n = (0,0,1) is (1,0,0). This fallback applies to the elevation 90 deg
rows; the CSV header is fixed, so those rows are recognisable by t = (1,0,0)
where every other row has t = (-1,0,0), and the script says so on stdout.
Each coefficient is reported "free" (O ignored) and "occ" (O inside the expectation).

GEOMETRY (every case)
---------------------
Target plate : square of side 1e-4 m centred at the origin in the plane z = 0,
               n = (0,0,1); p is sampled uniformly on it (the plate is so small
               that occlusion is effectively a function of direction only).
Occluder     : one rectangle in the plane x = 0.5 m, y in [-2, 2] m,
               z in [0, 0.5] m, blocking rays from either side (two-sided).
Ray test     : lambda = (p_x - 0.5) / c_x (c_x != 0), occluded iff lambda > 0,
               |p_y - lambda c_y| <= 2 and 0 <= p_z - lambda c_z <= 0.5.

CASES
-----
vhat = (cos e, 0, sin e); s in {2, 5, 8}; e in {90, 60, 30, 10, 0, -10} deg;
caseId = 1..18, s ascending, then e in that order. cos and sin are taken exactly
(0 or 1) at the multiples of 90 deg, so vhat at e = 90 is exactly (0,0,1).

ESTIMATORS
----------
Plain Monte Carlo: c = W + sigma*(z1, z2, z3), z_i independent standard normal,
drawn with random.Random(seed), seed = SEED_BASE + caseId. For each sample the
integrand x_j (one per coefficient, free and occ) is accumulated; the estimate is
the sample mean m_j = sum(x_j)/N and its standard error is
se_j = sqrt( (sum(x_j^2) - N m_j^2) / (N - 1) / N ).
Component z3 is drawn first; a molecule with n.c >= 0 contributes 0 to every
integrand, so z1, z2 and p are drawn only for molecules that hit. They are
independent of z3, so this is the same plain estimator at lower cost, and the
stream stays reproducible for a given seed.

SAMPLE SIZE
-----------
With Y = -n.c ~ Normal(q, sigma^2), q = -n.W = sin e, k = q/sigma, and
Phi/phi the standard normal cdf/pdf, the truncated moments are exact:
  E[Y   H] = sigma   ( k Phi(k) + phi(k) )
  E[Y^2 H] = sigma^2 ( (k^2 + 1) Phi(k) + k phi(k) )
  E[Y^4 H] = sigma^4 ( (k^4 + 6 k^2 + 3) Phi(k) + (k^3 + 5 k) phi(k) )
so the per-sample relative standard deviation of the CP_free integrand is
r = sqrt( E[Y^4 H] / E[Y^2 H]^2 - 1 ) and N = ceil(r^2 / target^2), rounded up
to a multiple of 100000 and clipped to [nmin, nmax] (defaults: target 2e-3,
nmin 1e6, nmax 1e8). The same formulas, with CT_free = 2 E[Y H] (W.t) (c.t is
independent of n.c because t lies in the plate plane), give an analytic
free-stream check (d) of the sampler; they say nothing about occlusion.

SELF-CHECKS (stdout; exit status 1 if any fails)
------------------------------------------------
(a) s = 8, e = 90: |CP_free - 2| <= 5 se + 0.02.
(b) s = 8, e = 90: the occluder only blocks arrival directions below ~45 deg
    elevation toward +x, ~11 sigma away, so CP_occ / CP_free ~ 1 (1e-3 allowed).
(c) every case: CN_occ <= CN_free + 3 se(CN_free).
(d) every case: MC free coefficients vs the closed forms above, |z| <= 5.

USAGE
-----
  python tpmc_check.py [--out tpmc_oracle.csv] [--target-rse 2e-3]
                       [--nmin 1000000] [--nmax 100000000]
Progress goes to stderr, self-check results to stdout.

Created 2026-09-15 by Claude Opus 5 (Anthropic).
"""

import argparse
import csv
import math
import random
import sys
import time

SEED_BASE = 20260915
S_LIST = (2.0, 5.0, 8.0)
E_LIST = (90.0, 60.0, 30.0, 10.0, 0.0, -10.0)

PLATE_SIDE = 1e-4          # m, target square centred at the origin, z = 0
OCC_X = 0.5                # m, occluder plane x = OCC_X
OCC_YMIN, OCC_YMAX = -2.0, 2.0
OCC_ZMIN, OCC_ZMAX = 0.0, 0.5
NORMAL = (0.0, 0.0, 1.0)

HEADER = ["caseId", "s", "elevation_deg", "vhat_x", "vhat_y", "vhat_z",
          "t_x", "t_y", "t_z", "N",
          "CN_free", "CN_free_se", "CP_free", "CP_free_se", "CT_free", "CT_free_se",
          "CN_occ", "CN_occ_se", "CP_occ", "CP_occ_se", "CT_occ", "CT_occ_se",
          "seed"]


def cos_sin_deg(e_deg):
    """cos and sin of an angle in degrees, exact at multiples of 90 deg."""
    if e_deg % 90.0 == 0.0:
        quadrant = int(round(e_deg / 90.0)) % 4
        return ((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))[quadrant]
    r = math.radians(e_deg)
    return math.cos(r), math.sin(r)


def tangent(vhat, n):
    """Unit tangential direction of g = -vhat; (fallback, flag)."""
    g = (-vhat[0], -vhat[1], -vhat[2])
    gn = g[0] * n[0] + g[1] * n[1] + g[2] * n[2]
    d = (g[0] - gn * n[0], g[1] - gn * n[1], g[2] - gn * n[2])
    norm = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
    if norm < 1e-12:
        e = (1.0, 0.0, 0.0)
        en = e[0] * n[0] + e[1] * n[1] + e[2] * n[2]
        d = (e[0] - en * n[0], e[1] - en * n[1], e[2] - en * n[2])
        norm = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        return (d[0] / norm, d[1] / norm, d[2] / norm), True
    return (d[0] / norm, d[1] / norm, d[2] / norm), False


def free_moments(q, sigma):
    """E[Y H], E[Y^2 H], E[Y^4 H] for Y ~ Normal(q, sigma^2), H = [Y > 0]."""
    k = q / sigma
    Phi = 0.5 * math.erfc(-k / math.sqrt(2.0))
    phi = math.exp(-0.5 * k * k) / math.sqrt(2.0 * math.pi)
    m1 = sigma * (k * Phi + phi)
    m2 = sigma ** 2 * ((k * k + 1.0) * Phi + k * phi)
    m4 = sigma ** 4 * ((k ** 4 + 6.0 * k * k + 3.0) * Phi + (k ** 3 + 5.0 * k) * phi)
    return m1, m2, m4


def choose_n(q, sigma, target, nmin, nmax):
    _, m2, m4 = free_moments(q, sigma)
    r2 = m4 / (m2 * m2) - 1.0
    n = int(math.ceil(r2 / (target * target)))
    n = ((n + 99999) // 100000) * 100000
    return max(nmin, min(nmax, n)), math.sqrt(r2)


def run_case(vhat, t, sigma, N, seed):
    """Plain Monte Carlo of the six coefficients; returns (means, ses)."""
    rng = random.Random(seed)
    gauss = rng.gauss
    rnd = rng.random
    wx, wy, wz = -vhat[0], -vhat[1], -vhat[2]
    tx, ty, tz = t
    side = PLATE_SIDE
    x0, ymin, ymax, zmin, zmax = OCC_X, OCC_YMIN, OCC_YMAX, OCC_ZMIN, OCC_ZMAX
    # free sums: a = -n.c, a^2, a^4, a*(c.t), (a*(c.t))^2 ; same for occ
    f1 = f2 = f4 = g1 = g2 = 0.0
    o1 = o2 = o4 = h1 = h2 = 0.0
    for _ in range(N):
        cz = gauss(wz, sigma)          # n = (0,0,1): n.c = c_z
        if cz >= 0.0:
            continue                   # H = 0: zero in every integrand
        cx = gauss(wx, sigma)
        cy = gauss(wy, sigma)
        px = (rnd() - 0.5) * side
        py = (rnd() - 0.5) * side
        a = -cz
        a2 = a * a
        act = a * (cx * tx + cy * ty + cz * tz)
        f1 += a
        f2 += a2
        f4 += a2 * a2
        g1 += act
        g2 += act * act
        if cx != 0.0:
            lam = (px - x0) / cx       # ray p - lam*c reaches the plane x = x0
            if lam > 0.0:
                zh = -lam * cz         # p_z = 0
                if zmin <= zh <= zmax:
                    yh = py - lam * cy
                    if ymin <= yh <= ymax:
                        continue       # occluded: O = 0
        o1 += a
        o2 += a2
        o4 += a2 * a2
        h1 += act
        h2 += act * act

    def mean_se(sum1, sum2, scale):
        m = sum1 / N
        var = (sum2 - sum1 * m) / (N - 1)
        return scale * m, scale * math.sqrt(max(var, 0.0) / N)

    return (mean_se(f1, f2, 1.0), mean_se(f2, f4, 2.0), mean_se(g1, g2, 2.0),
            mean_se(o1, o2, 1.0), mean_se(o2, o4, 2.0), mean_se(h1, h2, 2.0))


def fmt(x):
    return "%.17g" % (x + 0.0)       # + 0.0 turns -0.0 into 0.0


def default_out():
    here = __file__
    cut = max(here.rfind("/"), here.rfind("\\"))
    return (here[:cut + 1] if cut >= 0 else "") + "tpmc_oracle.csv"


def main():
    ap = argparse.ArgumentParser(description="TPMC self-shadowing oracle (stdlib only).")
    ap.add_argument("--out", default=default_out())
    ap.add_argument("--target-rse", type=float, default=2e-3)
    ap.add_argument("--nmin", type=int, default=1000000)
    ap.add_argument("--nmax", type=int, default=100000000)
    args = ap.parse_args()

    t_start = time.perf_counter()
    rows, results = [], []
    case_id = 0
    for s in S_LIST:
        sigma = 1.0 / (s * math.sqrt(2.0))
        for e in E_LIST:
            case_id += 1
            ce, se_ = cos_sin_deg(e)
            vhat = (ce, 0.0, se_)
            t, fallback = tangent(vhat, NORMAL)
            # mean of Y = -n.c is -n.W = n.vhat
            q = vhat[0] * NORMAL[0] + vhat[1] * NORMAL[1] + vhat[2] * NORMAL[2]
            N, r = choose_n(q, sigma, args.target_rse, args.nmin, args.nmax)
            seed = SEED_BASE + case_id
            tc = time.perf_counter()
            res = run_case(vhat, t, sigma, N, seed)
            sys.stderr.write("case %2d  s=%g e=%+5.1f  N=%9d  r=%7.3f  %6.1f s%s\n"
                             % (case_id, s, e, N, r, time.perf_counter() - tc,
                                "  (t fallback)" if fallback else ""))
            sys.stderr.flush()
            (cnf, cnf_se), (cpf, cpf_se), (ctf, ctf_se), \
                (cno, cno_se), (cpo, cpo_se), (cto, cto_se) = res
            rows.append([str(case_id), "%d" % s, "%d" % e,
                         fmt(vhat[0]), fmt(vhat[1]), fmt(vhat[2]),
                         fmt(t[0]), fmt(t[1]), fmt(t[2]), str(N),
                         fmt(cnf), fmt(cnf_se), fmt(cpf), fmt(cpf_se), fmt(ctf), fmt(ctf_se),
                         fmt(cno), fmt(cno_se), fmt(cpo), fmt(cpo_se), fmt(cto), fmt(cto_se),
                         str(seed)])
            wt = -vhat[0] * t[0] - vhat[1] * t[1] - vhat[2] * t[2]
            results.append(dict(case=case_id, s=s, e=e, N=N, fallback=fallback,
                                sigma=sigma, q=q, wt=wt, res=res))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    runtime = time.perf_counter() - t_start

    ok = True
    print("tpmc_check: %d cases written to %s, runtime %.1f s" % (len(rows), args.out, runtime))
    for r in results:
        if r["fallback"]:
            print("note: case %d (s=%g, e=%g): flow along the normal, t set to (1,0,0) "
                  "projected onto the plate plane" % (r["case"], r["s"], r["e"]))

    ref = [r for r in results if r["s"] == 8.0 and r["e"] == 90.0][0]
    cpf, cpf_se = ref["res"][1]
    cpo = ref["res"][4][0]
    pass_a = abs(cpf - 2.0) <= 5.0 * cpf_se + 0.02
    print("(a) s=8 e=90: CP_free = %.6f +- %.2e, |CP_free-2| = %.5f, bound %.5f -> %s"
          % (cpf, cpf_se, abs(cpf - 2.0), 5.0 * cpf_se + 0.02, "PASS" if pass_a else "FAIL"))
    ratio = cpo / cpf
    pass_b = abs(ratio - 1.0) <= 1e-3
    print("(b) s=8 e=90: CP_occ/CP_free = %.12f -> %s" % (ratio, "PASS" if pass_b else "FAIL"))
    pass_c = True
    worst_z = 0.0
    for r in results:
        (cnf, cnf_se) = r["res"][0]
        cno = r["res"][3][0]
        good = cno <= cnf + 3.0 * cnf_se
        pass_c = pass_c and good
        m1, m2, _ = free_moments(r["q"], r["sigma"])
        exact = (m1, 2.0 * m2, 2.0 * m1 * r["wt"])
        zs = []
        for j in range(3):
            est, se = r["res"][j]
            zs.append((est - exact[j]) / se if se > 0.0 else 0.0)
        worst_z = max(worst_z, max(abs(z) for z in zs))
        print("(c) case %2d s=%g e=%+4g: CN_occ/CN_free = %.6f %s | "
              "(d) z(CN,CP,CT) = %+.2f %+.2f %+.2f"
              % (r["case"], r["s"], r["e"], cno / cnf, "ok" if good else "FAIL",
                 zs[0], zs[1], zs[2]))
    pass_d = worst_z <= 5.0
    print("(c) CN_occ <= CN_free + 3 se in every case -> %s" % ("PASS" if pass_c else "FAIL"))
    print("(d) max |z| of MC free coefficients vs closed form = %.2f -> %s"
          % (worst_z, "PASS" if pass_d else "FAIL"))
    ok = pass_a and pass_b and pass_c and pass_d
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
