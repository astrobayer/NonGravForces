"""Independent oracle for NonGravForces.qc.detectJumps (pure Python, no numpy).

Exhaustive optimal partitioning (dynamic programming over EVERY admissible last
change point, no pruning) of a fixed-seed series with planted steps, for the
penalised least-squares objective the MATLAB unit minimises with PELT:

    sum over segments of RSS(segment) + beta * (number of change points)

RSS about a constant ("level") or about a straight line in the sample index
("slope"). PELT is exact, so the MATLAB change points must EQUAL these indices
(bitwise) when called with the same explicit 'Penalty' and 'MinSegment'.

Writes pelt_series.csv (i, x) and pelt_oracle.csv (model, beta, minSeg, cost,
changePoints as ';'-joined 1-based indices of the first sample of each new segment).
Run:  python pelt_check.py
"""
import csv
import os
import random


def rss_level(x, a, b):
    seg = x[a:b + 1]
    n = len(seg)
    m = sum(seg) / n
    return sum((v - m) ** 2 for v in seg)


def rss_slope(x, a, b):
    seg = x[a:b + 1]
    n = len(seg)
    ts = list(range(a + 1, b + 2))                     # 1-based sample index as regressor
    tm = sum(ts) / n
    xm = sum(seg) / n
    sxt = sum((t - tm) * (v - xm) for t, v in zip(ts, seg))
    stt = sum((t - tm) ** 2 for t in ts)
    slope = sxt / stt
    return sum((v - xm - slope * (t - tm)) ** 2 for t, v in zip(ts, seg))


def optimal_partition(x, beta, min_seg, cost):
    n = len(x)
    F = [float("inf")] * (n + 1)
    F[0] = -beta
    prev = [0] * (n + 1)
    for t in range(min_seg, n + 1):
        best = float("inf"); arg = 0
        for tau in range(0, t - min_seg + 1):
            if F[tau] == float("inf"):
                continue
            if tau > 0 and tau < min_seg:
                continue
            val = F[tau] + cost(x, tau, t - 1) + beta
            if val < best:
                best = val; arg = tau
        F[t] = best; prev[t] = arg
    cps = []
    tau = prev[n]
    while tau > 0:
        cps.append(tau + 1)                            # 1-based first sample of the new segment
        tau = prev[tau]
    return F[n], sorted(cps)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    rng = random.Random(20260907)
    K = 160
    x = [rng.gauss(0.0, 1.0) for _ in range(K)]
    for i in range(40, 90):
        x[i] += 4.0                                    # level step up at 41, down at 91
    for i in range(120, K):
        x[i] += 0.12 * (i - 119)                       # slope change at 121
    with open(os.path.join(here, "pelt_series.csv"), "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(["i", "x"])
        for i, v in enumerate(x):
            wr.writerow([i + 1, repr(v)])
    rows = []
    for model, cost, min_seg in (("level", rss_level, 2), ("slope", rss_slope, 3)):
        for beta in (10.0, 25.0, 60.0):
            c, cps = optimal_partition(x, beta, min_seg, cost)
            rows.append([model, repr(beta), min_seg, repr(c), ";".join(str(v) for v in cps)])
            print(model, beta, "cost %.6f" % c, "changepoints", cps)
    with open(os.path.join(here, "pelt_oracle.csv"), "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(["model", "beta", "minSeg", "cost", "changePoints"])
        wr.writerows(rows)


if __name__ == "__main__":
    main()
