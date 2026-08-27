"""Independent oracle for NonGravForces.qc.hampel (numpy only).

Writes tools/pyxcheck/hampel_oracle.csv with columns i, x, flagged for a
fixed-seed series with planted outliers; the MATLAB reference test reads it.
Run:  python hampel_check.py
"""
import numpy as np, csv, os
rng = np.random.default_rng(20260827)
K, W, t = 400, 21, 4.0
x = rng.standard_normal(K)
plant = {50: 9.0, 51: -8.0, 200: 15.0, 333: 7.5}
for i, v in plant.items():
    x[i] += v
h = (W - 1) // 2
flag = np.zeros(K, bool)
for i in range(K):
    w = x[max(0, i - h):min(K, i + h + 1)]
    m = np.median(w); s = 1.4826 * np.median(np.abs(w - m))
    flag[i] = abs(x[i] - m) > t * s if s > 0 else False
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hampel_oracle.csv")
with open(out, "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["i", "x", "flagged"])
    for i in range(K):
        wr.writerow([i + 1, repr(float(x[i])), int(flag[i])])
print(out, int(flag.sum()), "flagged")
