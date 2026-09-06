#!/usr/bin/env python3
"""Independent oracle for NonGravForces.env.eclipseFactor (conical model), pure Python.

The conical shadow function is the unocculted fraction of the solar disc as seen
from the satellite, with the Earth a sphere of radius R_E. Montenbruck & Gill
(2000, 3.4.2) give it in closed form as the overlap area of two discs of angular
radii a (Sun) and b (Earth) at angular separation c. This script does NOT use
that formula: it integrates the overlap numerically (composite Simpson over the
chord lengths of the Sun's disc that lie inside the Earth's disc), an independent
code path. The integrand has a kink where the two circles cross and square-root
endpoints where a chord vanishes; the range is split at the kink and each piece
is integrated in the variable u with x = lo + (hi - lo)(1 - cos(pi u))/2, which
makes the integrand smooth, so Simpson with N = 20000 per piece is accurate to
about 1e-12.

It also emits the exact limiting cases (full light, umbra, annular) and a
cylindrical-model column for the same geometries.

Output: eclipse_oracle.csv with columns
    case, rx, ry, rz, sx, sy, sz, R_E, R_S, nu_conical, nu_cylindrical
positions in metres (Earth-centred, Sun vector Earth -> Sun).

Written 2026-09-07 by Claude Fable 5.1 (Anthropic) from the geometry, not from the MATLAB.
"""
import csv, math, os

R_E = 6378137.0          # WGS-84 semi-major axis, the sphere the toolbox uses
R_S = 6.957e8            # IAU 2015 B3 nominal solar radius
AU = 149597870700.0


def angles(r, s):
    """Angular radii a (Sun), b (Earth) and separation c, seen from the satellite at r."""
    d = [s[i] - r[i] for i in range(3)]              # satellite -> Sun
    e = [-r[i] for i in range(3)]                    # satellite -> Earth centre
    nd = math.sqrt(sum(x * x for x in d)); ne = math.sqrt(sum(x * x for x in e))
    a = math.asin(R_S / nd)
    b = math.asin(R_E / ne)
    cosc = sum(d[i] * e[i] for i in range(3)) / (nd * ne)
    c = math.acos(max(-1.0, min(1.0, cosc)))
    return a, b, c, nd, ne


def overlap_area_numeric(a, b, c, n=20000):
    """Area of the intersection of two discs (radii a, b, centre distance c), numerically."""
    if c >= a + b:
        return 0.0
    if c <= abs(a - b):
        return math.pi * min(a, b) ** 2
    # Disc A centred at 0, disc B centred at c on the x axis. Integrate the chord overlap over x.
    x0 = max(-a, c - b); x1 = min(a, c + b)

    def chord(x):
        ya = math.sqrt(max(0.0, a * a - x * x))
        yb = math.sqrt(max(0.0, b * b - (x - c) ** 2))
        return 2.0 * min(ya, yb)

    xk = (a * a - b * b + c * c) / (2.0 * c)         # where the circles cross: the kink of min(ya, yb)
    cuts = [x0] + ([xk] if x0 < xk < x1 else []) + [x1]
    total = 0.0
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        L = hi - lo

        def f(u):
            x = lo + L * (1.0 - math.cos(math.pi * u)) / 2.0
            return chord(x) * L * math.pi * math.sin(math.pi * u) / 2.0

        h = 1.0 / n
        acc = f(0.0) + f(1.0)
        for k in range(1, n):
            acc += (4 if k % 2 else 2) * f(k * h)
        total += acc * h / 3.0
    return total


def nu_conical(r, s):
    a, b, c, nd, ne = angles(r, s)
    if ne > nd:                       # Earth farther than the Sun: no shadow
        return 1.0
    A = overlap_area_numeric(a, b, c)
    return 1.0 - A / (math.pi * a * a)


def nu_cylindrical(r, s):
    ns = math.sqrt(sum(x * x for x in s)); u = [x / ns for x in s]
    along = sum(r[i] * u[i] for i in range(3))
    if along >= 0:
        return 1.0
    perp = math.sqrt(sum((r[i] - along * u[i]) ** 2 for i in range(3)))
    return 0.0 if perp < R_E else 1.0


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, 'eclipse_oracle.csv')
    s = [AU, 0.0, 0.0]                              # Sun on +x
    rows = []
    h = 490e3; rr = R_E + h
    # sweep the anti-sun side: angle psi between -x and the satellite position, in the x-z plane
    psis = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 65.0, 66.0, 67.0, 67.5, 68.0, 68.2, 68.4, 68.6, 68.8,
            69.0, 69.2, 69.4, 69.6, 69.8, 70.0, 70.5, 71.0, 72.0, 75.0, 80.0, 90.0, 120.0, 180.0]
    for k, psi in enumerate(psis):
        p = math.radians(psi)
        r = [-rr * math.cos(p), 0.0, rr * math.sin(p)]
        rows.append((k, r, s))
    # a higher orbit, a GEO-like distance, and an out-of-plane case
    rows.append((len(rows), [-(R_E + 1500e3) * math.cos(math.radians(60)), 0.0, (R_E + 1500e3) * math.sin(math.radians(60))], s))
    rows.append((len(rows), [-42164e3 * math.cos(math.radians(8.6)), 0.0, 42164e3 * math.sin(math.radians(8.6))], s))
    rows.append((len(rows), [-rr * math.cos(math.radians(68.5)), rr * math.sin(math.radians(68.5)) * 0.6, rr * math.sin(math.radians(68.5)) * 0.8], s))
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'rx', 'ry', 'rz', 'sx', 'sy', 'sz', 'R_E', 'R_S', 'nu_conical', 'nu_cylindrical'])
        for k, r, sv in rows:
            w.writerow([k] + ['%.17g' % v for v in r] + ['%.17g' % v for v in sv] + ['%.17g' % R_E, '%.17g' % R_S,
                       '%.17g' % nu_conical(r, sv), '%.17g' % nu_cylindrical(r, sv)])
    print('wrote', out, len(rows), 'cases')


if __name__ == '__main__':
    main()
