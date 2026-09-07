"""Independent oracle for NonGravForces.qc.hampel (standard library only).

WHAT IS FROZEN AND WHAT IS RE-RUN
  hampel_oracle.csv holds a 400-sample series with four planted outliers (columns
  i, x) and the flags a correct Hampel filter produces on it (column flagged).
  The SERIES is committed data: it was drawn once on 2026-08-27 from
  numpy.random.default_rng(20260827).standard_normal(400) with +9.0, -8.0, +15.0
  and +7.5 added at indices 50, 51, 200 and 333 (0-based). It is not redrawn,
  exactly as the real L1B excerpts under data/samples are not re-fetched.
  The FLAGS are re-run here, from the committed series, by an implementation that
  shares no code with the MATLAB unit, and compared with the committed column.

  This matters because until 2026-09-07 this file imported numpy, which is
  installed on neither tier (the owner's host has Python 3.14.7 without it, and
  the ubuntu audit job installs nothing). The one oracle in the tree could
  therefore not be regenerated on any machine that runs the gate: it had become
  a claim rather than an instrument. The rewrite below is stdlib-only, so the
  flags column is now genuinely reproducible.

ALGORITHM (Hampel 1974; Pearson 2002)
  Over a centred window of W samples, truncated at the edges,
      med  = median of the window
      sigma = 1.4826 * median(|window - med|)
      flag  = |x[i] - med| > t * sigma       (never flagged when sigma == 0)
  with W = 21 and t = 4.0 for this fixture.

Run:  python hampel_check.py        (verifies, and rewrites the file identically)
Exit code 1 and a diff report if the recomputed flags disagree with the committed ones.
"""
import csv
import os
import sys

W = 21
THRESHOLD = 4.0
HALF = (W - 1) // 2


def median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        raise ValueError("median of an empty window")
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def hampel_flags(x, window=W, threshold=THRESHOLD):
    """Flags of a centred, edge-truncated Hampel filter. Pure Python, no dependencies."""
    half = (window - 1) // 2
    n = len(x)
    flags = []
    for i in range(n):
        w = x[max(0, i - half):min(n, i + half + 1)]
        med = median(w)
        sigma = 1.4826 * median([abs(v - med) for v in w])
        flags.append(abs(x[i] - med) > threshold * sigma if sigma > 0 else False)
    return flags


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "hampel_oracle.csv")
    if not os.path.isfile(path):
        sys.exit(
            "hampel_oracle.csv is missing. The series is FROZEN COMMITTED DATA drawn once with "
            "numpy on 2026-08-27 (see this file's docstring) and is not redrawn here; restore it "
            "from git rather than regenerating it."
        )
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    x = [float(r["x"]) for r in rows]
    committed = [int(r["flagged"]) for r in rows]

    recomputed = [int(v) for v in hampel_flags(x)]

    differing = [i for i, (a, b) in enumerate(zip(committed, recomputed)) if a != b]
    if differing:
        for i in differing[:20]:
            print(f"  row {i + 1}: committed {committed[i]}, recomputed {recomputed[i]}, x = {x[i]!r}")
        sys.exit(f"{len(differing)} of {len(rows)} flags disagree with the committed oracle")

    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["i", "x", "flagged"])
        for i, (xi, fi) in enumerate(zip(x, recomputed)):
            wr.writerow([i + 1, repr(xi), fi])
    print(f"{path} {len(rows)} rows, {sum(recomputed)} flagged; recomputed flags match the committed column")


if __name__ == "__main__":
    main()
