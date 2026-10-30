#!/usr/bin/env python3
"""attitude_quality_check.py -- independent oracle for NonGravForces.geom.attitudeQuality (N-12).

    python3 tools/pyxcheck/attitude_quality_check.py [day]      # day: 2020-01-01 (default), -02, -03

Reads the committed 30 min excerpts data/samples/multiday/SCA1B_<day>_C_04.txt and IMU1B_<day>_C_04.txt
with plain string splitting (not gracefoReading) and recomputes what the MATLAB unit computes, from the
definitions and not from its code:

  - the SCA1B quaternion (file: inertial -> SRF, scalar first) conjugated once to BODY -> inertial and
    turned into the matrix R whose columns are the body axes in the inertial frame;
  - over each 1 s interval the rotation vector of R_k^T R_k+1 divided by dt: the mean body rate, in the
    convention skew(omega) = R^T dR/dt;
  - each gyro's accumulated FiltAng (deg -> rad) interpolated linearly to the SCA1B epochs and
    differenced over the same intervals;
  - per gyro the least-squares fit rate = a . omega + b and its residual, the residual as an exact
    projection (modified Gram-Schmidt on normalised columns): on an Earth-pointing arc omega is nearly
    constant, the fit's parameters are not separable and only the residual is a result.

THE QUATERNION CONVENTION IS MEASURED, NOT ASSUMED: both matrix conventions (R and R^T) are evaluated and
the one whose rates the gyros fit is reported (a body-fixed gyro cannot follow a rate expressed in the
inertial frame over a third of a revolution); the other is printed beside it.

Standard library only (CLAUDE.md section 7). Prints the numbers tests/tNgfAttitude.m compares with.
Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
"""
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def records(path):
    rows, body = [], False
    with open(path, encoding="latin-1") as f:
        for ln in f:
            if body:
                if ln.strip():
                    rows.append(ln.split())
            elif ln.startswith("# End of YAML header"):
                body = True
    return rows


def quat_to_matrix(w, x, y, z):
    """Hamilton, scalar first: the matrix of the active rotation by q."""
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]


def transpose(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def rotvec(m):
    w = [(m[2][1] - m[1][2]) / 2, (m[0][2] - m[2][0]) / 2, (m[1][0] - m[0][1]) / 2]
    s = math.sqrt(sum(c * c for c in w))
    th = math.atan2(s, (m[0][0] + m[1][1] + m[2][2] - 1) / 2)
    f = th / s if s > 0 else 1.0
    return [c * f for c in w]


def solve(a, b):
    """Gaussian elimination with partial pivoting, small dense systems."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][k] * x[k] for k in range(r + 1, n))) / m[r][r]
    return x


def lstsq(X, y):
    p = len(X[0])
    ata = [[sum(row[i] * row[j] for row in X) for j in range(p)] for i in range(p)]
    aty = [sum(row[i] * yy for row, yy in zip(X, y)) for i in range(p)]
    return solve(ata, aty)


def project_out(X, y):
    """Residual of y after removing its orthogonal projection on the columns of X (modified Gram-Schmidt
    on column-normalised X), and the 2-norm condition estimate is not needed: the residual is exact."""
    cols = [[row[j] for row in X] for j in range(len(X[0]))]
    cols = [[v / math.sqrt(sum(c * c for c in col)) for v in col] for col in cols]
    basis = []
    for col in cols:
        v = col[:]
        for b in basis:
            d = sum(x * z for x, z in zip(v, b))
            v = [x - d * z for x, z in zip(v, b)]
        n = math.sqrt(sum(x * x for x in v))
        basis.append([x / n for x in v])
    r = y[:]
    for b in basis:
        d = sum(x * z for x, z in zip(r, b))
        r = [x - d * z for x, z in zip(r, b)]
    return r


def interp(tq, t, v):
    out, j = [], 0
    for x in tq:
        if x < t[0] or x > t[-1]:
            out.append(None)
            continue
        while j + 1 < len(t) and t[j + 1] < x:
            j += 1
        if t[j] == x or j + 1 >= len(t):
            out.append(v[j])
        else:
            out.append(v[j] + (v[j + 1] - v[j]) * (x - t[j]) / (t[j + 1] - t[j]))
    return out


def main():
    day = sys.argv[1] if len(sys.argv) > 1 else "2020-01-01"
    d = os.path.join(ROOT, "data", "samples", "multiday")
    sca = records(os.path.join(d, f"SCA1B_{day}_C_04.txt"))
    imu = records(os.path.join(d, f"IMU1B_{day}_C_04.txt"))
    t0 = float(sca[0][0])
    ts = [float(r[0]) - t0 for r in sca]
    # file quaternion inertial -> SRF; conjugate once
    q = [(float(r[3]), -float(r[4]), -float(r[5]), -float(r[6])) for r in sca]
    gy = {}
    for r in imu:
        t = float(r[0]) + float(r[1]) * 1e-6 - t0
        gy.setdefault(int(r[4]), ([], []))
        gy[int(r[4])][0].append(t)
        gy[int(r[4])][1].append(math.radians(float(r[5])))
    results = {}
    for name, conv in (("R", lambda m: m), ("R^T", transpose)):
        R = [conv(quat_to_matrix(*qq)) for qq in q]
        om = []
        for k in range(len(R) - 1):
            dt = ts[k + 1] - ts[k]
            om.append([c / dt for c in rotvec(matmul(transpose(R[k]), R[k + 1]))])
        res = {}
        for g in sorted(gy):
            th = interp(ts, gy[g][0], gy[g][1])
            X, y = [], []
            for k in range(len(ts) - 1):
                if th[k] is None or th[k + 1] is None:
                    continue
                dt = ts[k + 1] - ts[k]
                X.append(om[k] + [1.0])
                y.append((th[k + 1] - th[k]) / dt)
            c = lstsq(X, y)
            r = project_out(X, y)
            a = c[:3]
            n = math.sqrt(sum(v * v for v in a))
            res[g] = dict(axis=[v / n for v in a], scale=n, bias=c[3],
                          rms=math.sqrt(sum(v * v for v in r) / len(r)),
                          rmsRate=math.sqrt(sum(v * v for v in y) / len(y)), n=len(y))
        mag = sorted(math.sqrt(sum(c * c for c in o)) for o in om)
        results[name] = (res, mag[len(mag) // 2], om)
    best = min(results, key=lambda k: sum(v["rms"] for v in results[k][0].values()))
    for name in ("R", "R^T"):
        res, med, om = results[name]
        print(f"convention {name}{'  <- the gyros fit this one' if name == best else ''}: median |omega| {med:.6e} rad/s")
        for g, v in res.items():
            print(f"  gyro {g}: n {v['n']}  axis [{v['axis'][0]:+.9f} {v['axis'][1]:+.9f} {v['axis'][2]:+.9f}]"
                  f"  scale {v['scale']:.9f}  bias {v['bias']:+.6e}  rms {v['rms']:.12e}  rmsRate {v['rmsRate']:.12e}")
    om = results[best][2]
    k = 100
    print(f"omega[{k}] ({best}) = [{om[k][0]:+.15e} {om[k][1]:+.15e} {om[k][2]:+.15e}] rad/s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
