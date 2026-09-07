#!/usr/bin/env python3
"""Independent restatement and evaluation of the SOLAARS-CF shadow function (Robertson 2015).

Source: R. V. Robertson (2015), Highly Physical Solar Radiation Pressure Modeling During
Penumbra Transitions, PhD thesis, Virginia Tech, Section 3.7 (eqs 3.80-3.86, Table 3.2).
Coefficients transcribed 2026-09-07 from the thesis PDF (VTechWorks handle 10919/73537)
with three independent text extractions (layout, table and raw modes of pdftotext) that
agree on every value; the sigmoid centres they give at rR = 2000 km (6.371, 6.413,
6.398 x 1e6 m) match Fig. 3.28 (transition between 6360 and 6430 km).

    nu = [1 + a1 + a2 + a1 tanh(a3 (rE - a4)) + a2 tanh(a5 (rE - a6)) + tanh(a7 (rE - a8))]
         / (2 + 2 a1 + 2 a2)                                                          (3.86)
    rR = -r . R_hat (distance behind the terminator plane, R_hat Earth -> Sun),
    rE = |r - (r . R_hat) R_hat| with the ecliptic-parallel part scaled by sO = R_e / R_p   (3.80-3.85)
    a1, a3, a5, a7, a8 = b1 exp(b2 rR) + b3 exp(b4 rR);  a2, a6 = b1 rR + b2;  a4 = b1 rR^b2 + b3
    all distances in units of 1e6 m; the fit covers rR from 2000 km to 45000 km.

Writes solaars_cf_coefficients.csv (the table) and solaars_oracle.csv (nu on a grid of rR, rE),
both read by tests/tNgfEnv.m as the independent path for the formula arithmetic.

Written 2026-09-07 by Claude Fable 5.1 (Anthropic) from the thesis, not from the MATLAB.
"""
import csv, math, os

TABLE = {   # name: (form, b1, b2, b3, b4)
    'a1': ('exp2', 0.1715, -0.1423, 0.01061, -0.01443),
    'a2': ('lin', 0.008162, 0.3401, None, None),
    'a3': ('exp2', 260.9, -0.4661, 27.81, -0.009437),
    'a4': ('pow', -0.006119, 1.176, 6.385, None),
    'a5': ('exp2', 87.56, -0.09188, 19.30, -0.01089),
    'a6': ('lin', 0.002047, 6.409, None, None),
    'a7': ('exp2', 61.98, -0.1629, 27.87, -0.02217),
    'a8': ('exp2', 6.413, -0.0002593, -0.01479, -0.1318),
}


def coeff(name, rR):
    form, b1, b2, b3, b4 = TABLE[name]
    if form == 'exp2':
        return b1 * math.exp(b2 * rR) + b3 * math.exp(b4 * rR)
    if form == 'lin':
        return b1 * rR + b2
    return b1 * rR ** b2 + b3


def nu(rR, rE):
    a = {k: coeff(k, rR) for k in TABLE}
    num = 1 + a['a1'] + a['a2'] + a['a1'] * math.tanh(a['a3'] * (rE - a['a4'])) \
        + a['a2'] * math.tanh(a['a5'] * (rE - a['a6'])) + math.tanh(a['a7'] * (rE - a['a8']))
    return num / (2 + 2 * a['a1'] + 2 * a['a2'])


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'solaars_cf_coefficients.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['name', 'form', 'b1', 'b2', 'b3', 'b4'])
        for k, (form, b1, b2, b3, b4) in TABLE.items():
            w.writerow([k, form, b1, b2, '' if b3 is None else b3, '' if b4 is None else b4])
    with open(os.path.join(here, 'solaars_oracle.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rR_Mm', 'rE_Mm', 'nu'])
        for rR in [2.0, 2.56, 3.0, 5.0, 10.0, 20.0, 45.0]:
            for rE in [6.00, 6.30, 6.35, 6.37, 6.38, 6.39, 6.40, 6.41, 6.42, 6.43, 6.45, 6.50, 6.60, 7.00]:
                w.writerow(['%.17g' % rR, '%.17g' % rE, '%.17g' % nu(rR, rE)])
    print('nu at rR=2.0: rE 6.37 -> %.4f, 6.40 -> %.4f, 6.43 -> %.4f' % (nu(2.0, 6.37), nu(2.0, 6.40), nu(2.0, 6.43)))
    print('wrote solaars_cf_coefficients.csv and solaars_oracle.csv')


if __name__ == '__main__':
    main()
