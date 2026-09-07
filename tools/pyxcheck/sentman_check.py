"""Independent kinetic-theory oracle for NonGravForces.kernel.gsiCoefficients (pure Python).

The flat-plate coefficients are obtained by NUMERICAL integration of the momentum
carried to a plate by a drifting Maxwellian, never through the closed erf/exp
form the MATLAB unit uses. Conventions (identical to the unit's help text):

  plate normal n points out of the surface into the gas; the gas streams relative
  to the plate with bulk speed U along g = -gamma n + l t, gamma = cos(theta),
  l = sin(theta), theta the angle between n and the satellite's relative velocity
  v = -g (theta = 0: head-on, theta > 90 deg: leeward); t is the tangential part of g;
  s = U / cbar, cbar = sqrt(2 k T / m) the most probable thermal speed;
  forces per unit plate area are normalised by (1/2) rho U^2:
    Cp   normal coefficient (pushes the plate along -n), incident + diffuse re-emission
    Ctau tangential coefficient (along t, the tangential part of the gas stream)
    Cp_spec  normal coefficient of a perfectly SPECULAR plate (incident doubled, no re-emission)
    CD, CL   drag (along the flow) and lift (perpendicular, in the n-t plane) implied by Cp, Ctau
  re-emission (Sentman 1961 / Koppenwallner): diffuse half-Maxwellian at T_re with
    T_re = T_i - alpha (T_i - T_w),  T_i = m U^2 / (4 k)   (energy-equivalent incident temperature),
  whose mean normal effusion speed is sqrt(pi k T_re / (2 m)).

Everything reduces to one-dimensional Gaussian integrals over the normal velocity
component cz < 0 (tangential components integrate analytically to the bulk
value), evaluated by composite Simpson on a 12-sigma window.

Writes sentman_oracle.csv with columns
  s, theta_deg, alpha, twRatio (T_w/T_inf), Cp, Ctau, Cp_spec, CD, CL
Run:  python sentman_check.py
"""
import csv
import math
import os

SIGMA = 1.0 / math.sqrt(2.0)          # std of one velocity component when cbar = 1


def simpson(f, a, b, n=40000):
    if b <= a:
        return 0.0
    h = (b - a) / n
    acc = f(a) + f(b)
    for i in range(1, n):
        acc += f(a + i * h) * (4 if i % 2 else 2)
    return acc * h / 3.0


def incident(s, gamma):
    """Number flux N/(n cbar), normal momentum flux p_i/(n m cbar^2), for cbar = 1."""
    mu = -s * gamma                      # mean normal velocity of the gas (negative = towards the plate)
    lo = mu - 12.0 * SIGMA
    if lo >= 0.0:
        return 0.0, 0.0
    hi = 0.0
    g = lambda c: math.exp(-((c - mu) ** 2) / (2 * SIGMA ** 2)) / (SIGMA * math.sqrt(2 * math.pi))
    nflux = simpson(lambda c: -c * g(c), lo, hi)
    pflux = simpson(lambda c: c * c * g(c), lo, hi)
    return nflux, pflux


def coefficients(s, theta_deg, alpha, tw_ratio):
    gamma = math.cos(math.radians(theta_deg))
    ell = math.sin(math.radians(theta_deg))
    nflux, pflux = incident(s, gamma)
    half_rho_u2 = 0.5 * s * s             # (1/2) rho U^2 with rho = n m = 1, cbar = 1
    # tangential momentum: incident molecules carry the bulk tangential velocity U l on average
    tau_i = nflux * s * ell
    # re-emission: T_i/T_inf = s^2/2 (cbar^2 = 2kT/m = 1 -> m/(4k) U^2 / T_inf = s^2/2)
    ti = 0.5 * s * s
    tre = ti - alpha * (ti - tw_ratio)
    vbar_re = math.sqrt(math.pi * tre / 4.0)   # sqrt(pi k T_re/(2 m)) with 2kT_inf/m = 1
    p_re = nflux * vbar_re
    cp = (pflux + p_re) / half_rho_u2
    ctau = tau_i / half_rho_u2
    cp_spec = 2.0 * pflux / half_rho_u2
    cd = cp * gamma + ctau * ell
    cl = cp * ell - ctau * gamma
    return cp, ctau, cp_spec, cd, cl


def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sentman_oracle.csv")
    rows = []
    for s in (0.5, 2.0, 5.0, 8.0, 12.0):
        for theta in (0.0, 30.0, 60.0, 85.0, 90.0, 120.0, 180.0):
            for alpha in (0.0, 0.85, 1.0):
                for tw in (0.3,):        # T_w = 300 K at T_inf = 1000 K
                    rows.append((s, theta, alpha, tw) + coefficients(s, theta, alpha, tw))
    with open(out, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["s", "theta_deg", "alpha", "twRatio", "Cp", "Ctau", "Cp_spec", "CD", "CL"])
        for r in rows:
            wr.writerow([repr(float(v)) for v in r])
    # two sanity anchors printed for the record: head-on hyperthermal cold plate -> Cp = 2, static gas -> p = n k T
    print(out, len(rows), "rows")
    print("s=12 head-on alpha=1 Tw/T=0.3: Cp=%.6f (2 + sqrt(pi)*sqrt(0.3)/12*... re-emission small)" % coefficients(12, 0, 1, 0.3)[0])
    cp0 = coefficients(0.001, 0, 1, 1.0)[0] * 0.5 * 0.001 ** 2   # p/(n m cbar^2) -> 1/2 = n k T /(n m cbar^2) as s -> 0
    print("static gas T_w = T: p/(n m cbar^2) = %.6f (expected 0.5)" % cp0)


if __name__ == "__main__":
    main()
