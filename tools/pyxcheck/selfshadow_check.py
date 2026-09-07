#!/usr/bin/env python3
"""Independent oracle for NonGravForces.geom.selfShadow: a pure-Python ray cast on the GRACE-FO plate polygons.

The plate polygons are rebuilt here from the same exterior dimensions the MATLAB factory
uses (JPL D-56935 Figs 2-3: length 3.1225 m, zenith width 0.6950 m, height 0.7210 m, array
flanks along the 50 deg slant from the zenith edge to the nadir plane and 0.2282913/3.1225 m
beyond it; the nadir half-width is where the flank crosses the nadir plane, 0.9525 m). The geometry is therefore a shared assumption; the ray casting -
sampling, ray-plane intersection, point-in-polygon - is an independent code path with
a finer grid (nS x nS = 160 x 160 points per plate).

Convention: dir is the unit vector from the body to the light source; a plate with
n . dir <= 0 is unlit (fraction 0); the exposed fraction is the share of sample points
whose ray towards the source meets no other plate in front of them.

Output: selfshadow_oracle.csv, rows (case, dx, dy, dz, f_1 .. f_8) in the panel order
Front Rear StarboardOuter StarboardInner PortOuter PortInner Nadir Zenith.

Written 2026-09-07 by Claude Fable 5.1 (Anthropic) from the geometry, not from the MATLAB.
"""
import csv, math, os

Lx, Hz, Wt = 3.1225, 0.7210, 0.6950
s50, c50 = 0.766044, 0.642787
wInner = 0.2282913 / Lx
tCross = Hz / s50                    # slant length from the zenith edge to the nadir plane
yNadir = Wt / 2 + tCross * c50       # nadir half-width where the flank crosses the nadir plane


def sub(a, b): return [a[i] - b[i] for i in range(3)]
def add(a, b): return [a[i] + b[i] for i in range(3)]
def mul(a, s): return [a[i] * s for i in range(3)]
def dot(a, b): return sum(a[i] * b[i] for i in range(3))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def norm(a): return math.sqrt(dot(a, a))
def unit(a): n = norm(a); return [x / n for x in a]


def plates():
    """(name, normal, corners) with corners as lists of 3-vectors."""
    P = []
    sec = [(-yNadir, Hz/2), (yNadir, Hz/2), (Wt/2, -Hz/2), (-Wt/2, -Hz/2)]
    P.append(("Front", [1, 0, 0], [[Lx/2, y, z] for y, z in sec]))
    P.append(("Rear", [-1, 0, 0], [[-Lx/2, y, z] for y, z in sec]))
    for sgn, name_o, name_i in [(1, "StarboardOuter", "StarboardInner"), (-1, "PortOuter", "PortInner")]:
        a = (sgn * Wt/2, -Hz/2)                                   # zenith edge
        b = (sgn * yNadir, Hz/2)                                  # nadir edge
        e = (b[0] + wInner * sgn * c50, b[1] + wInner * s50)      # end of the overhang
        outer = [[-Lx/2, a[0], a[1]], [Lx/2, a[0], a[1]], [Lx/2, e[0], e[1]], [-Lx/2, e[0], e[1]]]
        inner = [[-Lx/2, b[0], b[1]], [Lx/2, b[0], b[1]], [Lx/2, e[0], e[1]], [-Lx/2, e[0], e[1]]]
        P.append((name_o, [0, sgn * s50, -c50], outer))
        P.append((name_i, [0, -sgn * s50, c50], inner))
    P.append(("Nadir", [0, 0, 1], [[-Lx/2, -yNadir, Hz/2], [Lx/2, -yNadir, Hz/2], [Lx/2, yNadir, Hz/2], [-Lx/2, yNadir, Hz/2]]))
    P.append(("Zenith", [0, 0, -1], [[-Lx/2, Wt/2, -Hz/2], [Lx/2, Wt/2, -Hz/2], [Lx/2, -Wt/2, -Hz/2], [-Lx/2, -Wt/2, -Hz/2]]))
    return P


def frame(normal, corners):
    """In-plane orthonormal axes and the 2-D polygon, wound counter-clockwise about the normal."""
    o = corners[0]
    e1 = unit(sub(corners[1], o)); e2 = cross(normal, e1)
    poly = [(dot(sub(c, o), e1), dot(sub(c, o), e2)) for c in corners]
    area = 0.5 * sum(poly[i][0] * poly[(i+1) % len(poly)][1] - poly[(i+1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))
    if area < 0:
        poly = poly[::-1]; corners = corners[::-1]
    return o, e1, e2, poly, corners


def inside(poly, x, y):
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]; x1, y1 = poly[(i+1) % n]
        if (x1 - x0) * (y - y0) - (y1 - y0) * (x - x0) < -1e-12:
            return False
    return True


def exposed(P, d, nS=160):
    frames = [frame(nrm, cor) for _, nrm, cor in P]
    out = []
    for k, (name, nk, ck) in enumerate(P):
        if dot(nk, d) <= 0:
            out.append(0.0); continue
        o, e1, e2, poly, cor = frames[k]
        # bilinear sampling of the quadrilateral (all GRACE-FO plates are quadrilaterals), area-weighted
        hit = 0.0; wsum = 0.0
        for i in range(nS):
            u = (i + 0.5) / nS
            for j in range(nS):
                w = (j + 0.5) / nS
                p = [(1-u)*(1-w)*cor[0][t] + u*(1-w)*cor[1][t] + u*w*cor[2][t] + (1-u)*w*cor[3][t] for t in range(3)]
                du = [-(1-w)*cor[0][t] + (1-w)*cor[1][t] + w*cor[2][t] - w*cor[3][t] for t in range(3)]
                dw = [-(1-u)*cor[0][t] - u*cor[1][t] + u*cor[2][t] + (1-u)*cor[3][t] for t in range(3)]
                jac = norm(cross(du, dw))
                wsum += jac
                occ = False
                for jn, (name2, n2, cor2) in enumerate(P):
                    if jn == k: continue
                    den = dot(n2, d)
                    if abs(den) < 1e-12: continue
                    o2, f1, f2, poly2, _ = frames[jn]
                    t = dot(n2, sub(o2, p)) / den
                    if t <= 1e-9: continue
                    x = add(p, mul(d, t))
                    rel = sub(x, o2)
                    if inside(poly2, dot(rel, f1), dot(rel, f2)):
                        occ = True; break
                if occ:
                    hit += jac
        out.append(1.0 - hit / wsum)
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    P = plates()
    dirs = [unit([0, 0, -1]), unit([0, math.sin(math.radians(70)), math.cos(math.radians(70))]),
            unit([0, math.sin(math.radians(85)), math.cos(math.radians(85))]), unit([0, -math.sin(math.radians(85)), math.cos(math.radians(85))]),
            unit([1, 0, 0]), unit([0.3, 0.8, 0.52]), unit([-0.5, -0.7, 0.5])]
    with open(os.path.join(here, 'selfshadow_oracle.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'dx', 'dy', 'dz'] + ['f_' + name for name, _, _ in P])
        for k, d in enumerate(dirs):
            fr = exposed(P, d)
            w.writerow([k] + ['%.17g' % x for x in d] + ['%.17g' % x for x in fr])
            print('case', k, ' '.join('%.4f' % x for x in fr))
    print('wrote selfshadow_oracle.csv')


if __name__ == '__main__':
    main()
