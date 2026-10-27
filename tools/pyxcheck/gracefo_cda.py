#!/usr/bin/env python3
"""gracefo_cda.py -- the effective drag area Cd*A of the GRACE-FO macro model under Sentman.

    python3 -I tools/pyxcheck/gracefo_cda.py

Splits the free-molecular drag of geom.graceFO's eight plates (areas and normals of the
GRACE-FO handbook, flight along -x, no shadowing) into the frontal plate and the edge-on
shear of the 16.7 m^2 of side, nadir and zenith plates, through the independent kernel
mirror tools/pyxcheck/ngfkernel, for three compositions and two accommodations. Written for
HANDOVER O-19 (2026-10-27): how much of the "drag twice too large" against ACT1B and GROOPS
is the gas-surface model rather than the density. Standard library only.

Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ngfkernel import drag  # noqa: E402

S50, C50 = math.sin(math.radians(50)), math.cos(math.radians(50))
PLATES = [("Front", 0.9551567, (1, 0, 0)), ("Rear", 0.9551567, (-1, 0, 0)),
          ("StarboardOuter", 3.1554792, (0, S50, -C50)), ("StarboardInner", 0.2282913, (0, -S50, C50)),
          ("PortOuter", 3.1554792, (0, -S50, -C50)), ("PortInner", 0.2282913, (0, S50, C50)),
          ("Nadir", 6.0711120, (0, 0, 1)), ("Zenith", 2.1673620, (0, 0, -1))]
U = 7600.0
FRONT = 0.9551567
CASES = [("O only, 750 K", [0, 3e13, 0, 0, 0, 0, 0, 0], 750.0),
         ("O + He, 700 K", [1.5e13, 2.5e13, 1e11, 0, 0, 0, 0, 0], 700.0),
         ("He-rich, 650 K", [3e13, 1.5e13, 0, 0, 0, 0, 0, 0], 650.0)]


def main():
    plates = [{"area": a, "normal": n, "centroid": (0, 0, 0), "vis": (0.3, 0.2, 0.5), "ir": (0.8, 0.1, 0.1)}
              for _, a, n in PLATES]
    print("%-16s %5s %9s %9s %9s %12s" % ("composition", "alpha", "CdA m^2", "frontal", "edge-on", "Cd on 0.955"))
    for label, species, t_k in CASES:
        q = 0.5 * sum(drag.species_mass_densities(species)) * U * U
        for alpha in (0.85, 1.0):
            f = drag.drag_panel_force(plates, species, t_k, (-U, 0.0, 0.0), "sentman", alpha=alpha, tw_k=300.0)
            fx = abs(sum(v[0] for v in f))
            fr = abs(f[0][0] + f[1][0])
            print("%-16s %5.2f %9.3f %9.3f %9.3f %12.2f" % (label, alpha, fx / q, fr / q, (fx - fr) / q, fx / q / FRONT))
    print("constantCd 2.4 on the frontal plate: CdA %.3f m^2" % (2.4 * FRONT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
