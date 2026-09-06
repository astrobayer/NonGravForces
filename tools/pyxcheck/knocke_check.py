#!/usr/bin/env python3
"""Independent oracle for NonGravForces.env.earthRadiance / radiationSources (Knocke model), pure Python.

Knocke, Ries, Tapley (1988), Earth radiation pressure effects on satellites, Table 1:
    a(phi, t) = a0 + a1 cos(w (t - t0)) P1(sin phi) + a2 P2(sin phi)
    e(phi, t) = e0 + e1 cos(w (t - t0)) P1(sin phi) + e2 P2(sin phi)
    a0 = 0.34, a1 = 0.10, a2 = 0.29; e0 = 0.68, e1 = -0.07, e2 = -0.18; w = 2 pi / 365.25 d; t0 = 22 December.
Each surface element is a Lambertian source. Irradiance at the satellite from element dA
at distance rho, element normal n, satellite direction u (from element to satellite):
    albedo    dE = a S (n.sun)+ (n.u)+ dA / (pi rho^2)      (S = solar irradiance at Earth)
    longwave  dE = e (S/4) (n.u)+ dA / (pi rho^2)
This script integrates over a fine equal-angle grid (0.25 deg) of the whole sphere with the
visibility gate (n.u > 0) applied per element - a code path independent of the MATLAB tiles.

Output: knocke_oracle.csv, one row per case:
    case, uniform (1: a = a0, e = e0 everywhere), days_since_t0, rx, ry, rz, sunx, suny, sunz (unit), S_W_m2, E_albedo_W_m2, E_ir_W_m2
and a second file knocke_coefficients.csv restating Table 1 so the MATLAB test reads the
coefficients from an independent file rather than from the unit under test.

Written 2026-09-07 by Claude Fable 5.1 (Anthropic) from the paper's equations, not from the MATLAB.
"""
import csv, math, os

R_E = 6378137.0
A0, A1, A2 = 0.34, 0.10, 0.29
E0, E1, E2 = 0.68, -0.07, -0.18
OMEGA = 2.0 * math.pi / 365.25


UNIFORM = False   # True: a = A0, e = E0 everywhere (the point-source limit); set per case


def albedo(lat, days):
    if UNIFORM:
        return A0
    x = math.sin(lat)
    return A0 + A1 * math.cos(OMEGA * days) * x + A2 * 0.5 * (3 * x * x - 1)


def emissivity(lat, days):
    if UNIFORM:
        return E0
    x = math.sin(lat)
    return E0 + E1 * math.cos(OMEGA * days) * x + E2 * 0.5 * (3 * x * x - 1)


def irradiance(r, sun, S, days, step_deg=0.25):
    """Total albedo and longwave irradiance (W/m^2) at satellite position r (m), Sun unit vector sun."""
    nlat = int(round(180 / step_deg)); nlon = int(round(360 / step_deg))
    dlon = 2 * math.pi / nlon
    Ea = 0.0; Ei = 0.0
    for i in range(nlat):
        lat0 = -math.pi / 2 + i * math.pi / nlat; lat1 = lat0 + math.pi / nlat
        latc = 0.5 * (lat0 + lat1)
        area = R_E * R_E * dlon * (math.sin(lat1) - math.sin(lat0))
        a = albedo(latc, days); e = emissivity(latc, days)
        cl = math.cos(latc); sl = math.sin(latc)
        for j in range(nlon):
            lonc = -math.pi + (j + 0.5) * dlon
            n = (cl * math.cos(lonc), cl * math.sin(lonc), sl)
            p = (R_E * n[0], R_E * n[1], R_E * n[2])
            rho = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
            d2 = rho[0] ** 2 + rho[1] ** 2 + rho[2] ** 2
            d = math.sqrt(d2)
            cos_u = (n[0] * rho[0] + n[1] * rho[1] + n[2] * rho[2]) / d
            if cos_u <= 0.0:
                continue
            g = cos_u * area / (math.pi * d2)
            Ei += e * (S / 4.0) * g
            cos_s = n[0] * sun[0] + n[1] * sun[1] + n[2] * sun[2]
            if cos_s > 0.0:
                Ea += a * S * cos_s * g
    return Ea, Ei


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    S = 1361.0
    cases = [
        # (days since t0, r [m], sun unit)
        (0.0,   (R_E + 490e3, 0.0, 0.0),                       (1.0, 0.0, 0.0)),          # sub-solar, solstice phase 0
        (100.0, (0.0, 0.0, R_E + 490e3),                       (1.0, 0.0, 0.0)),          # over the north pole, terminator
        (200.0, (-(R_E + 490e3), 0.0, 0.0),                    (1.0, 0.0, 0.0)),          # anti-solar: albedo 0
        (45.0,  ((R_E + 490e3) * math.cos(0.7), 0.0, (R_E + 490e3) * math.sin(0.7)), (math.cos(0.4), math.sin(0.4), 0.0)),
        (300.0, (0.0, R_E + 1000e3, 0.0),                      (0.0, 1.0, 0.0)),          # higher orbit, sub-solar
        (10.0,  (R_E * 12.0, 0.0, 0.0),                        (1.0, 0.0, 0.0)),          # far field, zonal model
    ]
    uniform_cases = [
        (0.0,   (R_E * 12.0, 0.0, 0.0),                        (1.0, 0.0, 0.0)),          # uniform sphere, far field, full phase
        (0.0,   (R_E * 12.0 * math.cos(1.0), 0.0, R_E * 12.0 * math.sin(1.0)), (1.0, 0.0, 0.0)),  # uniform, far field, phase 1 rad
        (0.0,   (R_E + 490e3, 0.0, 0.0),                       (1.0, 0.0, 0.0)),          # uniform, LEO, full phase
        (0.0,   (R_E * 100.0, 0.0, 0.0),                       (1.0, 0.0, 0.0)),          # uniform, 100 R: far-field limit to 1e-2
        (0.0,   (R_E * 100.0 * math.cos(1.0), 0.0, R_E * 100.0 * math.sin(1.0)), (1.0, 0.0, 0.0)),  # uniform, 100 R, phase 1 rad
    ]
    with open(os.path.join(here, 'knocke_oracle.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'uniform', 'days_since_t0', 'rx', 'ry', 'rz', 'sunx', 'suny', 'sunz', 'S_W_m2', 'E_albedo_W_m2', 'E_ir_W_m2'])
        global UNIFORM
        allcases = [(False, c) for c in cases] + [(True, c) for c in uniform_cases]
        for k, (uni, (days, r, sun)) in enumerate(allcases):
            UNIFORM = uni
            Ea, Ei = irradiance(r, sun, S, days)
            w.writerow([k, int(uni), '%.17g' % days] + ['%.17g' % v for v in r] + ['%.17g' % v for v in sun] +
                       ['%.17g' % S, '%.17g' % Ea, '%.17g' % Ei])
            print('case', k, 'uniform', int(uni), 'E_albedo %.6f W/m2  E_ir %.6f W/m2' % (Ea, Ei))
        UNIFORM = False
    with open(os.path.join(here, 'knocke_coefficients.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name', 'value'])
        for name, val in [('a0', A0), ('a1', A1), ('a2', A2), ('e0', E0), ('e1', E1), ('e2', E2), ('period_days', 365.25), ('t0', '1981-12-22')]:
            w.writerow([name, val])
    print('wrote knocke_oracle.csv and knocke_coefficients.csv')


if __name__ == '__main__':
    main()
