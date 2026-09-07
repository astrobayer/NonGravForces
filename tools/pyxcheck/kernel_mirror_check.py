"""Parity driver for the independent kernel mirror (standard library only).

Evaluates tools/pyxcheck/ngfkernel/ -- a mirror of the NonGravForces kernel derived
from the physics, not transcribed from the MATLAB -- over a fixed case list and writes
one CSV per kernel unit. tests/tNgfMirror.m reads those files, feeds the SAME inputs to
the MATLAB kernel and compares. Inputs and expected outputs travel in the same file
(rule R2), so a case can never be read with the wrong expectation.

Plate geometry and optics are INPUT DATA here: kernel_mirror_plates.csv and
kernel_mirror_bodies.csv were exported once from NonGravForces.geom (2026-09-07,
GRACE-C from the JPL handbook and a unit cube with a diagonal inertia). The ruling of
2026-08-27 mirrors the KERNEL only, so geom is not under test and its output is data.

No randomness: every case value is a fixed number or a deterministic function of the
case index, so the files regenerate byte for byte on any host.

Run:  python kernel_mirror_check.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "ngfkernel"))

from gsi import gsi_coefficients          # noqa: E402
from drag import drag_panel_force         # noqa: E402
from radiation import radiation_panel_force  # noqa: E402
from thermal import thermal_recoil        # noqa: E402
from assemble import assemble             # noqa: E402
from torques import body_torques          # noqa: E402

# MSIS species order and molar masses [g/mol]: He O N2 O2 Ar H N anomalous-O
MOLAR_G = [4.002602, 15.9994, 28.0134, 31.9988, 39.948, 1.00794, 14.0067, 15.9994]
N_A = 6.02214076e23
GM_EARTH = 3.986004418e14


def load_bodies():
    plates = {}
    with open(os.path.join(HERE, "kernel_mirror_plates.csv"), newline="") as f:
        for r in csv.DictReader(f):
            plates.setdefault(r["body"], []).append({
                "area": float(r["area_m2"]),
                "normal": (float(r["nx"]), float(r["ny"]), float(r["nz"])),
                "centroid": (float(r["cx"]), float(r["cy"]), float(r["cz"])),
                "vis": (float(r["absVis"]), float(r["specVis"]), float(r["diffVis"])),
                "ir": (float(r["emisIR"]), float(r["specIR"]), float(r["diffIR"])),
            })
    bodies = {}
    with open(os.path.join(HERE, "kernel_mirror_bodies.csv"), newline="") as f:
        for r in csv.DictReader(f):
            bodies[r["body"]] = {
                "cm": (float(r["cmx"]), float(r["cmy"]), float(r["cmz"])),
                "inertia": [[float(r["I%d%d" % (i, j)]) for j in (1, 2, 3)] for i in (1, 2, 3)],
                "plates": plates[r["body"]],
            }
    return bodies


def unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


def w(path, header, rows):
    with open(os.path.join(HERE, path), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(header)
        for row in rows:
            wr.writerow([repr(v) if isinstance(v, float) else v for v in row])
    print("%-28s %5d rows" % (path, len(rows)))


# ----------------------------------------------------------------- gsiCoefficients
def cases_gsi():
    rows = []
    ss = [0.05, 0.5, 2.0, 5.0, 8.0, 12.0, 25.0]
    thetas = [0.0, 0.3, math.pi / 4, 1.2, math.pi / 2, 2.0, 2.8, math.pi]
    alphas = [0.0, 0.5, 0.85, 1.0]
    temps = [(300.0, 1000.0), (250.0, 700.0), (600.0, 1200.0)]
    for s in ss:
        for th in thetas:
            for a in alphas:
                for tw, ti in temps:
                    cp, ct = gsi_coefficients("sentman", s, th, a, tw, ti)
                    rows.append(["sentman", s, th, a, tw, ti, "", "", cp, ct])
    for s in [0.5, 8.0, 25.0]:
        for th in thetas:
            for sp in [0.0, 0.25, 1.0]:
                cp, ct = gsi_coefficients("maxwell", s, th, 0.85, 300.0, 1000.0, None, sp)
                rows.append(["maxwell", s, th, 0.85, 300.0, 1000.0, "", sp, cp, ct])
    for cd in [2.2, 3.0]:
        for th in thetas:
            cp, ct = gsi_coefficients("constantCd", None, th, None, None, None, cd)
            rows.append(["constantCd", "", th, "", "", "", cd, "", cp, ct])
    w("mirror_gsi.csv",
      ["model", "s", "theta", "alpha", "tw_k", "tinf_k", "cd", "specular", "cp", "ctau"], rows)


# ----------------------------------------------------------------- dragPanelForce
def gas(kind, scale):
    """Two thermospheres: oxygen-dominated (490 km) and a denser nitrogen-rich one."""
    if kind == "O":
        n = [2.0e13, 8.0e13, 4.0e11, 1.0e10, 2.0e9, 5.0e10, 1.0e12, 1.0e12]
    else:
        n = [1.0e12, 6.0e13, 9.0e14, 4.0e13, 6.0e12, 1.0e10, 2.0e11, 0.0]
    return [v * scale for v in n]


def mass_density(species):
    return sum(n * (m * 1e-3 / N_A) for n, m in zip(species, MOLAR_G))


def cases_drag(bodies):
    rows = []
    vdirs = [(1.0, 0.0, 0.0), (1.0, 0.05, -0.02), (0.0, 1.0, 0.0), (0.3, -0.9, 0.32), (0.0, 0.0, 1.0)]
    ci = 0
    for bname in ("graceFO_C", "cube"):
        plates = bodies[bname]["plates"]
        for vd in vdirs:
            for kind, scale, t_k in (("O", 1.0, 950.0), ("N", 0.4, 780.0)):
                for model, alpha, tw, cd, spec in (
                    ("sentman", 0.85, 300.0, None, None),
                    ("sentman", 1.0, 250.0, None, None),
                    ("maxwell", 0.85, 300.0, None, 0.2),
                    ("constantCd", None, None, 2.2, None),
                ):
                    ci += 1
                    if ci % 3 == 0:
                        shadow = [0.5 + 0.5 * math.cos(k + ci) ** 2 for k in range(len(plates))]
                    else:
                        shadow = [1.0] * len(plates)
                    sp = gas(kind, scale)
                    vrel = tuple(7600.0 * c for c in unit(vd))
                    kw = dict(alpha=alpha, tw_k=tw, cd=cd, specular=spec, shadow=shadow)
                    if model == "constantCd":
                        kw["cd_rho"] = mass_density(sp)
                    f = drag_panel_force(plates, sp, t_k, vrel, model, **kw)
                    for k in range(len(plates)):
                        rows.append([ci, bname, model,
                                     "" if alpha is None else alpha, "" if tw is None else tw,
                                     "" if cd is None else cd, "" if spec is None else spec,
                                     t_k, vrel[0], vrel[1], vrel[2]] + list(sp) +
                                    [shadow[k], k + 1, f[k][0], f[k][1], f[k][2]])
    w("mirror_drag.csv",
      ["caseId", "body", "model", "alpha", "tw_k", "cd", "specular", "t_k", "vx", "vy", "vz"] +
      ["d%d" % i for i in range(1, 9)] + ["shadow", "plate", "fx", "fy", "fz"], rows)


# ------------------------------------------------------------ radiationPanelForce
def cases_radiation(bodies):
    rows = []
    ci = 0
    sundirs = [(1.0, 0.0, 0.0), (0.6, 0.0, 0.8), (-0.3, 0.5, 0.81), (0.0, 0.0, -1.0)]
    for bname in ("graceFO_C", "cube"):
        plates = bodies[bname]["plates"]
        for sd in sundirs:
            for shadowed in (False, True):
                ci += 1
                srcs = [
                    {"dir": unit(sd), "p": 4.5606e-6, "band": "vis"},
                    {"dir": unit((-0.2, 0.1, -1.0)), "p": 1.83e-7, "band": "vis"},
                    {"dir": unit((0.1, -0.15, -1.0)), "p": 7.4e-7, "band": "ir"},
                ]
                if shadowed:
                    shadow = [[0.25 + 0.25 * ((k + m) % 3) for m in range(3)] for k in range(len(plates))]
                else:
                    shadow = [[1.0] * 3 for _ in range(len(plates))]
                forces, absorbed = radiation_panel_force(plates, srcs, shadow)
                for k in range(len(plates)):
                    flat = []
                    for m in range(3):
                        flat += [srcs[m]["dir"][0], srcs[m]["dir"][1], srcs[m]["dir"][2],
                                 srcs[m]["p"], srcs[m]["band"], shadow[k][m]]
                    rows.append([ci, bname] + flat + [k + 1, forces[k][0], forces[k][1], forces[k][2], absorbed[k]])
    hdr = ["caseId", "body"]
    for m in (1, 2, 3):
        hdr += ["s%ddx" % m, "s%ddy" % m, "s%ddz" % m, "s%dp" % m, "s%dband" % m, "s%dshadow" % m]
    hdr += ["plate", "fx", "fy", "fz", "absorbed"]
    w("mirror_radiation.csv", hdr, rows)


# ----------------------------------------------------------------- thermalRecoil
def cases_thermal(bodies):
    rows = []
    ci = 0
    for bname in ("graceFO_C", "cube"):
        plates = bodies[bname]["plates"]
        n = len(plates)
        for kind in ("uniform", "gradient", "coldface"):
            ci += 1
            if kind == "uniform":
                temps = [290.0] * n
            elif kind == "gradient":
                temps = [250.0 + 12.0 * k for k in range(n)]
            else:
                temps = [0.0 if k == 0 else 300.0 - 5.0 * k for k in range(n)]
            f = thermal_recoil(plates, temps)
            for k in range(n):
                rows.append([ci, bname, kind, k + 1, temps[k], f[k][0], f[k][1], f[k][2]])
    w("mirror_thermal.csv", ["caseId", "body", "pattern", "plate", "temp_k", "fx", "fy", "fz"], rows)


# --------------------------------------------------------------------- assemble
def cases_assemble(bodies):
    rows = []
    ci = 0
    for bname in ("graceFO_C", "cube"):
        b = bodies[bname]
        plates = b["plates"]
        n = len(plates)
        for rcm in (b["cm"], (0.0, 0.0, 0.0), (0.2, -0.1, 0.05)):
            ci += 1
            fk = [(1e-6 * math.sin(k + ci), 1e-6 * math.cos(2.0 * k - ci), 1e-6 * math.sin(0.5 * k + 2.0 * ci))
                  for k in range(n)]
            F, tau = assemble(fk, plates, rcm)
            for k in range(n):
                rows.append([ci, bname, rcm[0], rcm[1], rcm[2], k + 1,
                             fk[k][0], fk[k][1], fk[k][2], F[0], F[1], F[2], tau[0], tau[1], tau[2]])
    w("mirror_assemble.csv",
      ["caseId", "body", "rcmx", "rcmy", "rcmz", "plate", "inx", "iny", "inz",
       "Fx", "Fy", "Fz", "Tx", "Ty", "Tz"], rows)


# ------------------------------------------------------------------ bodyTorques
def cases_torques(bodies):
    rows = []
    ci = 0
    positions = [(6.87e6, 0.0, 0.0), (0.0, 6.87e6, 0.0), (0.0, 0.0, 6.87e6),
                 (3.96e6, 3.96e6, 3.96e6), (-4.86e6, 4.86e6, 0.0)]
    for bname in ("graceFO_C", "cube"):
        I = bodies[bname]["inertia"]
        for r in positions:
            for gg, mag in ((True, False), (False, True), (True, True), (False, False)):
                ci += 1
                B = (1.0e-5 * math.cos(ci), -2.0e-5 * math.sin(ci), 3.0e-5)
                dip = (0.3, 0.2, -0.1)
                tau = body_torques(I, r, GM_EARTH, B, dip, gg, mag)
                rows.append([ci, bname, r[0], r[1], r[2], GM_EARTH,
                             B[0], B[1], B[2], dip[0], dip[1], dip[2],
                             int(gg), int(mag), tau[0], tau[1], tau[2]])
    w("mirror_torques.csv",
      ["caseId", "body", "rx", "ry", "rz", "gm", "bx", "by", "bz", "mx", "my", "mz",
       "gravityGradient", "magnetic", "tx", "ty", "tz"], rows)


def main():
    bodies = load_bodies()
    print("bodies: " + ", ".join("%s (%d plates)" % (k, len(v["plates"])) for k, v in bodies.items()))
    cases_gsi()
    cases_drag(bodies)
    cases_radiation(bodies)
    cases_thermal(bodies)
    cases_assemble(bodies)
    cases_torques(bodies)


if __name__ == "__main__":
    main()
