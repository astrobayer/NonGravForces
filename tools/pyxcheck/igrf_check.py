#!/usr/bin/env python3
"""Independent oracle for the IGRF-14 geomagnetic main field and its secular variation, pure Python.

Written 2026-09-15 by Claude Opus 5 (Anthropic): written from a specification by an agent that
did not open the MATLAB toolbox (T9). No MATLAB file, no test and no toolbox source was read; the
coefficient reader, the Legendre functions and the synthesis are this file's own code.
Standard library only (argparse, csv, decimal, fractions, math, os, sys, time) - no numpy, no
scipy, no pyIGRF - so it runs on Python 3.12 (Ubuntu) and 3.13 (Windows).

INPUT
  data/igrf/igrf14coeffs.txt, the IAGA V-MOD IGRF-14 coefficient table. Format as read here:
  lines starting with '#' are comments; a label line starting 'c/s'; one header line
      g/h n m 1900.0 1905.0 ... 2020.0 2025.0 2025-30
  whose tokens after 'g/h n m' are the tabulated epochs and, last, the secular-variation label
  'Y0-YY' (valid from Y0 to the year ending in YY of the same century, here 2025 to 2030); then one
  row per coefficient 'g n m v_1 ... v_E sv' or 'h n m v_1 ... v_E sv' in nT (sv in nT/yr).
  Every row must carry E + 1 numbers; every g_n^m (0 <= m <= n) and h_n^m (1 <= m <= n) up to the
  file's maximum degree must be present exactly once, and no h_n^0; the script stops otherwise.
  The values are taken as the EXACT decimals printed (fractions.Fraction of the token).
  Observed in this file: all epochs are tabulated to n = 13, but every coefficient above n = 10 is
  0 in the epochs 1900-1995 (the script prints the highest non-zero degree per epoch), and the SV
  column is printed as explicit 0.0 for 9 <= n <= 13 - used as given.

MODEL (Alken et al. 2021, Earth Planets Space 73:49, eqs. 1-3; IAGA conventions)
  V(r, theta, phi, t) = a sum_{n=1}^{N} (a/r)^(n+1) sum_{m=0}^{n}
                          [g_n^m(t) cos(m phi) + h_n^m(t) sin(m phi)] P_n^m(cos theta)
  a = 6371.2 km; r geocentric radius [km]; theta geocentric colatitude; phi east longitude.
  P_n^m Schmidt quasi-normalised, no Condon-Shortley phase:
      P_n^m = sqrt(2 (n-m)! / (n+m)!) P_{n,m}  (m > 0),   P_n^0 = P_{n,0},
      P_{n,m}(x) = (1 - x^2)^(m/2) d^m/dx^m P_n(x).
  B = -grad V:
      B_r     = -dV/dr               =  sum (n+1) (a/r)^(n+2) [g cos m phi + h sin m phi] P_n^m
      B_theta = -(1/r) dV/dtheta     = -sum (a/r)^(n+2) [g cos m phi + h sin m phi] dP_n^m/dtheta
      B_phi   = -(1/(r sin theta)) dV/dphi
                                     =  sum (a/r)^(n+2) m [g sin m phi - h cos m phi] P_n^m / sin theta
  North X = -B_theta, East Y = B_phi, Down Z = -B_r (geocentric local frame);
  Earth-fixed B_xyz = B_r e_r + B_theta e_theta + B_phi e_phi,
      e_r = (sin t cos p, sin t sin p, cos t), e_theta = (cos t cos p, cos t sin p, -sin t),
      e_phi = (-sin p, cos p, 0);   F = |B|.
  Time (IAGA rule): for E_i <= y < E_{i+1} (two tabulated epochs)
      g(y) = g(E_i) + (g(E_{i+1}) - g(E_i)) (y - E_i) / (E_{i+1} - E_i);
  for E_last <= y <= 2030.0, g(y) = g(E_last) + SV (y - E_last). Outside [1900, 2030] the script
  stops (no extrapolation beyond the model's validity).

METHOD (why it is independent of the usual implementations)
  Legendre functions are NOT computed by a three-term (Gauss or Holmes-Featherstone) recursion.
  From Rodrigues' formula P_n(x) = 1/(2^n n!) d^n/dx^n (x^2 - 1)^n, expanding the binomial and
  differentiating n + m times term by term gives the closed polynomial
      Q_nm(x) = d^m P_n / dx^m = 1/(2^n n!) sum_{k=0}^{floor((n-m)/2)}
                                 (-1)^k C(n,k) (2n-2k)! / (n-m-2k)!  x^(n-m-2k),
  whose coefficients are built as exact integers/rationals (fractions.Fraction). Then, with
  s = sin theta and c = cos theta,
      P_n^m          = S_nm s^m Q_nm(c)
      dP_n^m/dtheta  = S_nm s^(m-1) [m c Q_nm(c) - s^2 Q'_nm(c)]   (m >= 1),   = -S_n0 s Q'_n0(c)
      P_n^m / s      = S_nm s^(m-1) Q_nm(c)                          (m >= 1; m = 0 drops out of B_phi)
  (only ds/dtheta = c and dc/dtheta = -s are used), so NO division by sin theta occurs anywhere:
  at a pole (s = 0 exactly) the m >= 2 terms vanish through s^(m-1) and m = 1 stays finite.
  cos(m phi), sin(m phi) come from de Moivre, (cos phi + i sin phi)^m, not from cos(m*phi).
  EVALUATION IS IN EXACT RATIONAL ARITHMETIC: the only inexact steps are
    * sin/cos of the input angles, by math.sin/math.cos after an EXACT reduction of the degree
      value to a multiple of 90 deg plus a remainder in [-45, 45] (so lat = +-90 gives cos = 0 and
      sin = +-1 exactly, like a degree-argument sine); cos theta = sin(lat), sin theta = cos(lat);
    * S_nm = sqrt(2 (n-m)!/(n+m)!), a 60-digit decimal square root (error ~1e-60);
    * the final conversion of each output to float (correctly rounded) and F = sqrt(|B|^2)
      through a 60-digit decimal.
  Everything else (coefficient interpolation, Legendre polynomials, powers of a/r, the sums, the
  frame rotation) is exact. The inputs year, r_km, lat_deg, lon_deg and a = 6371.2 are taken as
  the IEEE doubles a double-precision caller holds (so r = 6371.2 gives a/r = 1 exactly).
  The CSV therefore holds the model value for those doubles to within the rounding of the input
  trig (~1e-16 relative), independent of summation order.

SELF-CHECKS (printed, not written)
  1. Exact identities of the Legendre code, as Fraction equalities (True/False), for n = 1..13 at
     rational points on the unit circle s = 2t/(1+t^2), c = (1-t^2)/(1+t^2) and at both poles:
        sum_m (P_n^m)^2 = 1  (addition theorem),   sum_m (dP_n^m/dtheta)^2 = n(n+1)/2,
        sum_m m^2 (P_n^m / sin theta)^2 = n(n+1)/2
     (the last two follow from d^2/dtheta1 dtheta2 and d^2/dphi1 dphi2 of the addition theorem at
     coincidence, with P_n(cos gamma) = 1 - n(n+1) gamma^2 / 4 + O(gamma^4)); they use S_nm^2 as
     an exact rational, so they test the polynomials, the normalisation, the theta derivative and
     the pole-safe P/sin form exactly.
  2. Dipole closed form (plain float code) at case 1 (the equator) and at case 7's point:
        B_r = 2 (a/r)^3 [g10 cos t + (g11 cos p + h11 sin p) sin t]
        B_t = (a/r)^3 [g10 sin t - (g11 cos p + h11 sin p) cos t]
        B_p = (a/r)^3 [g11 sin p - h11 cos p]
     against this script's synthesis with N = 1.
  3. -grad V by central differences of this script's OWN potential function at case 3, in exact
     arithmetic: r +- h_r (h_r = 1e-6 km), and theta, phi rotated by +-h (h = 1e-9 rad, sin h and
     cos h by truncated series, error < 1e-50), so the difference quotient has no rounding error
     and only its O(h^2) truncation.
  4. At the poles (exact s = 0) B_xyz must not depend on longitude, and must be the limit of
     nearby latitudes.

OUTPUT: igrf_oracle.csv, one row per case (%.17g)
  caseId, year, r_km, lat_deg, lon_deg, N, X_nT, Y_nT, Z_nT, Bx_nT, By_nT, Bz_nT, F_nT
  lat_deg is GEOCENTRIC latitude (90 - colatitude). At the exact poles X and Y are the limits along
  the meridian lon_deg (they depend on lon_deg); Bx, By, Bz do not.

Run:  python igrf_check.py [--coeffs PATH] [--out PATH]
      (defaults: ../../data/igrf/igrf14coeffs.txt and igrf_oracle.csv next to this script)
"""
import argparse
import csv
import math
import os
import sys
import time
from decimal import Decimal, localcontext
from fractions import Fraction

DEC_PREC = 60
A_KM = Fraction(6371.2)            # the IEEE double nearest 6371.2 km, exactly

HEADER = ["caseId", "year", "r_km", "lat_deg", "lon_deg", "N",
          "X_nT", "Y_nT", "Z_nT", "Bx_nT", "By_nT", "Bz_nT", "F_nT"]

# caseId, decimal year, geocentric radius [km], geocentric latitude [deg], east longitude [deg], N
CASES = [
    (1, 2020.0, 6371.2, 0.0, 0.0, 13),
    (2, 2020.0, 6371.2, 45.0, 90.0, 13),
    (3, 2020.0, 6861.2, -60.5, 200.25, 13),
    (4, 2020.0, 6861.2, 89.999, 10.0, 13),
    (5, 2020.0, 6861.2, 90.0, 0.0, 13),
    (6, 2020.0, 6861.2, -90.0, 123.0, 13),
    (7, 2020.0, 6861.2, 30.0, -75.0, 1),
    (8, 2022.5, 6800.0, 10.0, 50.0, 13),
    (9, 2027.25, 6800.0, -20.0, 300.0, 13),
    (10, 1900.0, 6371.2, 51.5, -0.1, 10),
    (11, 2020.0, 42164.0, 0.0, 75.0, 13),
    (12, 2019.999, 6861.2, 5.0, 5.0, 13),
]


# ---------------------------------------------------------------------------------------------
# Coefficient file
# ---------------------------------------------------------------------------------------------
class CoefficientTable:
    """Parsed IGRF table: epochs (Fraction), SV validity, and per-coefficient value lists."""

    def __init__(self, path):
        self.path = path
        self.epochs = None
        self.sv_start = None
        self.sv_end = None
        self.rows = {}               # (kind, n, m) -> [v_1 .. v_E, sv] as Fraction
        self._read(path)
        self.nmax = max(n for (_, n, _) in self.rows)
        self._validate()

    def _read(self, path):
        with open(path, "r") as f:
            for lineno, line in enumerate(f, 1):
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                tok = s.split()
                if tok[0] == "c/s":
                    continue
                if tok[0] == "g/h":
                    if tok[1:3] != ["n", "m"]:
                        raise ValueError("line %d: header must start 'g/h n m'" % lineno)
                    labels = tok[3:]
                    if "-" not in labels[-1]:
                        raise ValueError("line %d: last header label must be the SV span 'Y0-YY'"
                                         % lineno)
                    self.epochs = [Fraction(x) for x in labels[:-1]]
                    y0, y1 = labels[-1].split("-")
                    self.sv_start = Fraction(y0)
                    if len(y1) == 2:
                        self.sv_end = Fraction((int(Fraction(y0)) // 100) * 100 + int(y1))
                    else:
                        self.sv_end = Fraction(y1)
                    continue
                if tok[0] in ("g", "h"):
                    if self.epochs is None:
                        raise ValueError("line %d: data row before the 'g/h n m' header" % lineno)
                    n, m = int(tok[1]), int(tok[2])
                    vals = [Fraction(x) for x in tok[3:]]
                    if len(vals) != len(self.epochs) + 1:
                        raise ValueError("line %d: %d values, expected %d epochs + SV"
                                         % (lineno, len(vals), len(self.epochs)))
                    key = (tok[0], n, m)
                    if key in self.rows:
                        raise ValueError("line %d: duplicate coefficient %s" % (lineno, key))
                    self.rows[key] = vals
                    continue
                raise ValueError("line %d: unrecognised line %r" % (lineno, s[:40]))

    def _validate(self):
        if self.epochs is None:
            raise ValueError("no 'g/h n m' header line")
        if any(b <= a for a, b in zip(self.epochs, self.epochs[1:])):
            raise ValueError("epochs must increase")
        if self.sv_start != self.epochs[-1] or self.sv_end <= self.sv_start:
            raise ValueError("SV span must start at the last epoch")
        expected = set()
        for n in range(1, self.nmax + 1):
            for m in range(0, n + 1):
                expected.add(("g", n, m))
                if m > 0:
                    expected.add(("h", n, m))
        if set(self.rows) != expected:
            raise ValueError("coefficient set incomplete or with extras: missing %s, extra %s"
                             % (sorted(expected - set(self.rows))[:5],
                                sorted(set(self.rows) - expected)[:5]))

    def highest_nonzero_degree(self, column):
        return max((n for (_, n, _), v in self.rows.items() if v[column] != 0), default=0)

    def at(self, year):
        """Coefficient dict {(kind, n, m): Fraction} at a decimal year (Fraction), IAGA rule."""
        ep = self.epochs
        if year < ep[0] or year > self.sv_end:
            raise ValueError("year %s outside [%s, %s]" % (float(year), float(ep[0]),
                                                           float(self.sv_end)))
        out = {}
        if year >= ep[-1]:
            dt = year - ep[-1]
            for key, v in self.rows.items():
                out[key] = v[-2] + v[-1] * dt
            return out
        i = max(k for k in range(len(ep) - 1) if ep[k] <= year)
        w = (year - ep[i]) / (ep[i + 1] - ep[i])
        for key, v in self.rows.items():
            out[key] = v[i] + (v[i + 1] - v[i]) * w
        return out


# ---------------------------------------------------------------------------------------------
# Legendre functions from Rodrigues' formula (exact coefficients)
# ---------------------------------------------------------------------------------------------
def rodrigues_derivative_poly(n, m):
    """Coefficients c[j] of x^j in d^m P_n/dx^m, exact, from Rodrigues' formula."""
    c = [Fraction(0)] * (n - m + 1)
    den0 = (2 ** n) * math.factorial(n)
    for k in range(0, (n - m) // 2 + 1):
        num = (-1) ** k * math.comb(n, k) * math.factorial(2 * n - 2 * k)
        c[n - m - 2 * k] += Fraction(num, den0 * math.factorial(n - m - 2 * k))
    return c


def poly_derivative(c):
    return [c[j] * j for j in range(1, len(c))] or [Fraction(0)]


def poly_eval(c, x):
    acc = Fraction(0)
    for coef in reversed(c):
        acc = acc * x + coef
    return acc


def decimal_sqrt_fraction(q):
    """sqrt of a non-negative Fraction as an exact Fraction of a DEC_PREC-digit decimal."""
    with localcontext() as ctx:
        ctx.prec = DEC_PREC
        d = (Decimal(q.numerator) / Decimal(q.denominator)).sqrt()
    return Fraction(d)


class Legendre:
    """Q_nm, Q'_nm (exact polynomials) and Schmidt factors S_nm for n <= nmax."""

    def __init__(self, nmax):
        self.nmax = nmax
        self.Q = {}
        self.dQ = {}
        self.S2 = {}
        self.S = {}
        for n in range(0, nmax + 1):
            for m in range(0, n + 1):
                q = rodrigues_derivative_poly(n, m)
                self.Q[n, m] = q
                self.dQ[n, m] = poly_derivative(q)
                if m == 0:
                    self.S2[n, m] = Fraction(1)
                    self.S[n, m] = Fraction(1)
                else:
                    self.S2[n, m] = Fraction(2 * math.factorial(n - m), math.factorial(n + m))
                    self.S[n, m] = decimal_sqrt_fraction(self.S2[n, m])

    def unnormalised(self, n, m, s, c):
        """(s^m Q, s^(m-1) [m c Q - s^2 Q'] or -s Q', s^(m-1) Q or None) without S_nm."""
        q = poly_eval(self.Q[n, m], c)
        dq = poly_eval(self.dQ[n, m], c)
        if m == 0:
            return q, -s * dq, None
        sm1 = s ** (m - 1)
        return sm1 * s * q, sm1 * (m * c * q - s * s * dq), sm1 * q


# ---------------------------------------------------------------------------------------------
# Trigonometry and synthesis
# ---------------------------------------------------------------------------------------------
def sin_cos_deg(x_deg):
    """(sin x, cos x) for x in degrees: exact reduction to 90 k + rem, rem in [-45, 45]."""
    k = int(round(x_deg / 90.0))
    rem = x_deg - 90.0 * k
    t = math.radians(rem)
    s, c = math.sin(t), math.cos(t)
    k %= 4
    if k == 0:
        return s, c
    if k == 1:
        return c, -s
    if k == 2:
        return -s, -c
    return -c, s


def multiple_angles(cp, sp, nmax):
    """cos(m phi), sin(m phi), m = 0..nmax, from (cos phi + i sin phi)^m, exact."""
    cosm, sinm = [Fraction(1)], [Fraction(0)]
    re, im = Fraction(1), Fraction(0)
    for _ in range(nmax):
        re, im = re * cp - im * sp, re * sp + im * cp
        cosm.append(re)
        sinm.append(im)
    return cosm, sinm


def field_spherical(gh, leg, N, r, st, ct, sp, cp):
    """(B_r, B_theta, B_phi) [nT], exact, for Fraction inputs; st >= 0 is sin(colatitude)."""
    rho = A_KM / r
    cosm, sinm = multiple_angles(cp, sp, N)
    br = bt = bp = Fraction(0)
    rho_n2 = rho * rho                  # becomes (a/r)^(n+2) inside the loop
    for n in range(1, N + 1):
        rho_n2 *= rho
        for m in range(0, n + 1):
            g = gh["g", n, m]
            h = gh["h", n, m] if m > 0 else Fraction(0)
            p, dp, p_over_s = leg.unnormalised(n, m, st, ct)
            S = leg.S[n, m]
            cs = g * cosm[m] + h * sinm[m]
            br += (n + 1) * rho_n2 * cs * S * p
            bt -= rho_n2 * cs * S * dp
            if m > 0:
                bp += rho_n2 * m * (g * sinm[m] - h * cosm[m]) * S * p_over_s
    return br, bt, bp


def potential(gh, leg, N, r, st, ct, sp, cp):
    """V [nT km], exact, straight from the definition (uses P only, no derivative)."""
    rho = A_KM / r
    cosm, sinm = multiple_angles(cp, sp, N)
    v = Fraction(0)
    rho_n1 = rho
    for n in range(1, N + 1):
        rho_n1 *= rho
        for m in range(0, n + 1):
            g = gh["g", n, m]
            h = gh["h", n, m] if m > 0 else Fraction(0)
            p = leg.S[n, m] * (st ** m) * poly_eval(leg.Q[n, m], ct)
            v += rho_n1 * (g * cosm[m] + h * sinm[m]) * p
    return A_KM * v


def to_cartesian(br, bt, bp, st, ct, sp, cp):
    bx = br * st * cp + bt * ct * cp - bp * sp
    by = br * st * sp + bt * ct * sp + bp * cp
    bz = br * ct - bt * st
    return bx, by, bz


def exact_norm_float(*comps):
    q = sum((c * c for c in comps), Fraction(0))
    with localcontext() as ctx:
        ctx.prec = DEC_PREC
        return float((Decimal(q.numerator) / Decimal(q.denominator)).sqrt())


def geometry(lat_deg, lon_deg):
    if not (-90.0 <= lat_deg <= 90.0):
        raise ValueError("latitude must lie in [-90, 90] deg")
    slat, clat = sin_cos_deg(lat_deg)
    slon, clon = sin_cos_deg(lon_deg)
    # colatitude theta = 90 - lat: sin theta = cos lat (>= 0), cos theta = sin lat
    return Fraction(clat), Fraction(slat), Fraction(slon), Fraction(clon)


def evaluate(table, leg, year, r_km, lat_deg, lon_deg, N):
    if not (1 <= N <= table.nmax):
        raise ValueError("N must lie in [1, %d]" % table.nmax)
    if not (r_km > 0.0):
        raise ValueError("r must be positive")
    gh = table.at(Fraction(year))
    st, ct, sp, cp = geometry(lat_deg, lon_deg)
    r = Fraction(r_km)
    br, bt, bp = field_spherical(gh, leg, N, r, st, ct, sp, cp)
    bx, by, bz = to_cartesian(br, bt, bp, st, ct, sp, cp)
    return {
        "Br": br, "Bt": bt, "Bp": bp,
        "X": -bt, "Y": bp, "Z": -br, "Bx": bx, "By": by, "Bz": bz,
        "F": exact_norm_float(br, bt, bp),
        "gh": gh, "geom": (st, ct, sp, cp), "r": r,
    }


# ---------------------------------------------------------------------------------------------
# Self-checks
# ---------------------------------------------------------------------------------------------
def check_identities(leg, nmax):
    points = []
    for t in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1), Fraction(3), Fraction(7)):
        points.append((2 * t / (1 + t * t), (1 - t * t) / (1 + t * t)))
    points.append((Fraction(0), Fraction(-1)))       # south pole (t -> infinity)
    ok = [True, True, True]
    for s, c in points:
        for n in range(1, nmax + 1):
            a1 = a2 = a3 = Fraction(0)
            for m in range(0, n + 1):
                p, dp, ps = leg.unnormalised(n, m, s, c)
                a1 += leg.S2[n, m] * p * p
                a2 += leg.S2[n, m] * dp * dp
                if m > 0:
                    a3 += leg.S2[n, m] * m * m * ps * ps
            half = Fraction(n * (n + 1), 2)
            ok[0] &= (a1 == 1)
            ok[1] &= (a2 == half)
            ok[2] &= (a3 == half)
    out = sys.stdout
    out.write("check 1  exact Legendre identities, n = 1..%d, %d points incl. both poles:\n"
              % (nmax, len(points)))
    out.write("         sum (P)^2 = 1: %s;  sum (dP/dtheta)^2 = n(n+1)/2: %s;  "
              "sum m^2 (P/sin)^2 = n(n+1)/2: %s\n" % tuple(ok))
    return all(ok)


def dipole_closed_form(gh, r_km, lat_deg, lon_deg):
    """Plain float closed-form dipole (B_r, B_theta, B_phi)."""
    a = 6371.2
    t = math.radians(90.0 - lat_deg)
    p = math.radians(lon_deg)
    g10, g11, h11 = float(gh["g", 1, 0]), float(gh["g", 1, 1]), float(gh["h", 1, 1])
    k = (a / r_km) ** 3
    e = g11 * math.cos(p) + h11 * math.sin(p)
    br = 2.0 * k * (g10 * math.cos(t) + e * math.sin(t))
    bt = k * (g10 * math.sin(t) - e * math.cos(t))
    bp = k * (g11 * math.sin(p) - h11 * math.cos(p))
    return br, bt, bp


def check_dipole(table, leg):
    worst = 0.0
    for (year, r_km, lat, lon, label) in ((2020.0, 6371.2, 0.0, 0.0, "case 1 point"),
                                          (2020.0, 6861.2, 30.0, -75.0, "case 7 point")):
        res = evaluate(table, leg, year, r_km, lat, lon, 1)
        ref = dipole_closed_form(res["gh"], r_km, lat, lon)
        got = (float(res["Br"]), float(res["Bt"]), float(res["Bp"]))
        d = [abs(x - y) for x, y in zip(got, ref)]
        scale = math.sqrt(sum(x * x for x in ref))
        worst = max(worst, max(d) / scale)
        sys.stdout.write("check 2  dipole closed form at %s: |dB_r|, |dB_theta|, |dB_phi| = "
                         "%.3g, %.3g, %.3g nT (max relative %.3g)\n"
                         % ((label,) + tuple(d) + (max(d) / scale,)))
    return worst


def check_gradient(table, leg, case):
    _, year, r_km, lat, lon, N = case
    res = evaluate(table, leg, year, r_km, lat, lon, N)
    gh, (st, ct, sp, cp), r = res["gh"], res["geom"], res["r"]
    hr = Fraction(1, 10 ** 6)
    h = Fraction(1, 10 ** 9)
    sh = h - h ** 3 / 6 + h ** 5 / 120
    ch = 1 - h ** 2 / 2 + h ** 4 / 24 - h ** 6 / 720

    def V(rr, s_t, c_t, s_p, c_p):
        return potential(gh, leg, N, rr, s_t, c_t, s_p, c_p)

    dv_dr = (V(r + hr, st, ct, sp, cp) - V(r - hr, st, ct, sp, cp)) / (2 * hr)
    dv_dt = (V(r, st * ch + ct * sh, ct * ch - st * sh, sp, cp)
             - V(r, st * ch - ct * sh, ct * ch + st * sh, sp, cp)) / (2 * h)
    dv_dp = (V(r, st, ct, sp * ch + cp * sh, cp * ch - sp * sh)
             - V(r, st, ct, sp * ch - cp * sh, cp * ch + sp * sh)) / (2 * h)
    num = (-dv_dr, -dv_dt / r, -dv_dp / (r * st))
    ana = (res["Br"], res["Bt"], res["Bp"])
    rel = [float(abs(a - b) / abs(a)) for a, b in zip(num, ana)]
    vec = float(sum(((a - b) ** 2 for a, b in zip(num, ana)), Fraction(0))
                / sum((a * a for a in ana), Fraction(0))) ** 0.5
    sys.stdout.write("check 3  -grad V (exact central differences, h_r = 1e-6 km, h = 1e-9 rad) "
                     "vs analytic B at case 3:\n         relative |dB_r|, |dB_theta|, |dB_phi| = "
                     "%.3g, %.3g, %.3g;  |dB| / |B| = %.3g\n" % (rel[0], rel[1], rel[2], vec))
    return vec


def rel_vec_diff(a, b):
    num = sum(((x - y) ** 2 for x, y in zip(a, b)), Fraction(0))
    den = sum((x * x for x in a), Fraction(0))
    return math.sqrt(float(num / den))


def check_poles(table, leg):
    def xyz(lat, lon):
        res = evaluate(table, leg, 2020.0, 6861.2, lat, lon, 13)
        return res["Bx"], res["By"], res["Bz"]

    north0, north10 = xyz(90.0, 0.0), xyz(90.0, 10.0)
    south0, south123 = xyz(-90.0, 0.0), xyz(-90.0, 123.0)
    near_n, near_s = xyz(90.0 - 1e-7, 0.0), xyz(-90.0 + 1e-7, 123.0)
    sys.stdout.write("check 4  poles: B_xyz relative change lon 0 -> 10 at N pole %.3g, "
                     "lon 0 -> 123 at S pole %.3g\n"
                     "         limit 1e-7 deg away (expect ~1e-9): N %.3g, S %.3g\n"
                     % (rel_vec_diff(north0, north10), rel_vec_diff(south0, south123),
                        rel_vec_diff(north0, near_n), rel_vec_diff(south123, near_s)))


# ---------------------------------------------------------------------------------------------
def fmt(x):
    return "%.17g" % x


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--coeffs", default=os.path.normpath(
        os.path.join(here, "..", "..", "data", "igrf", "igrf14coeffs.txt")),
        help="IGRF coefficient file")
    parser.add_argument("--out", default=os.path.join(here, "igrf_oracle.csv"),
                        help="output CSV path")
    args = parser.parse_args(argv)

    t0 = time.perf_counter()
    table = CoefficientTable(args.coeffs)
    leg = Legendre(table.nmax)
    sys.stdout.write("read %s: %d coefficients, degree %d, epochs %s..%s, SV %s-%s\n"
                     % (args.coeffs, len(table.rows), table.nmax, float(table.epochs[0]),
                        float(table.epochs[-1]), float(table.sv_start), float(table.sv_end)))
    degs = {}
    for i, e in enumerate(table.epochs):
        degs.setdefault(table.highest_nonzero_degree(i), []).append(float(e))
    for d in sorted(degs):
        sys.stdout.write("         highest non-zero degree %d: epochs %s..%s (%d)\n"
                         % (d, degs[d][0], degs[d][-1], len(degs[d])))
    sys.stdout.write("         highest non-zero degree in SV: %d\n"
                     % table.highest_nonzero_degree(len(table.epochs)))

    rows = []
    for case in CASES:
        cid, year, r_km, lat, lon, N = case
        res = evaluate(table, leg, year, r_km, lat, lon, N)
        rows.append([str(cid), fmt(year), fmt(r_km), fmt(lat), fmt(lon), str(N)]
                    + [fmt(float(res[k])) for k in ("X", "Y", "Z", "Bx", "By", "Bz")]
                    + [fmt(res["F"])])
    t_cases = time.perf_counter() - t0

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    sys.stdout.write("wrote %d cases to %s (%.3f s incl. reading)\n"
                     % (len(rows), args.out, t_cases))

    check_identities(leg, table.nmax)
    check_dipole(table, leg)
    check_gradient(table, leg, CASES[2])
    check_poles(table, leg)
    sys.stdout.write("total %.3f s\n" % (time.perf_counter() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
