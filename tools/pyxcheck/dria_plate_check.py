"""Closed-form drag coefficient of a flat plate under diffuse reflection with incomplete
accommodation, from Walker, Mehta & Koller (2014) -- an oracle this project did not write.

WHY THIS EXISTS. Every other oracle for the drag kernel is derived here: the standard-library
mirror under ngfkernel/ and Sentman's own incident closed forms. Both are this project's reading
of the physics, so neither can catch a shared misreading. Walker et al. print a closed-form drag
coefficient for a flat plate, and at COMPLETE ACCOMMODATION it is the same quantity this toolbox's
normal coefficient computes, in a formulation that shares none of its algebra.

  [W14b] Walker A., Mehta P., Koller J., "Different Implementations of Diffuse Reflection with
         Incomplete Accommodation for Drag Coefficient Modeling", J. Spacecraft Rockets 51(5)
         1522-1532, 2014, DOI 10.2514/1.A32668.

  eq. (5)   T_ki = m U^2 / (3 k_B)                     incident kinetic temperature
  eq. (10)  T_kr = T_ki (1 - a) + a T_s                reflected kinetic temperature
  eq. (12)  CD   = (2 + 1/s^2) erf(s)
                 + 2/(sqrt(pi) s) exp(-s^2)
                 + sqrt(pi)/s sqrt(T_kr / T_inf)       flat plate, DRIA

WHERE IT IS AN IDENTITY AND WHERE IT IS NOT. This toolbox forms the re-emitted momentum from
Sentman's speed ratio in the Doornbos (2011) formulation,

  v_re/v_in = sqrt( (1/2) [ 1 + a (4 k T_w / (m U^2) - 1) ] ),

while [W14b] eq. (12) carries sqrt(T_kr/T_inf)/s. Written over the same variables the two are

  ours      sqrt( (1/2)(1 - a) + 2 a k T_w / (m U^2) )
  [W14b]    sqrt( (2/3)(1 - a) + 2 a k T_s / (m U^2) )

so they are the SAME NUMBER at a = 1 and differ only in how the unaccommodated fraction is
carried: 1/2 from the energy-flux formulation, 2/3 from a three-degree-of-freedom kinetic
temperature. The a = 1 column is therefore an identity to rounding and is asserted as one; the
rest of the table is a MEASUREMENT of that modelling difference, recorded and not asserted.

OUTPUT COLUMNS
  caseId, s, alpha, Tw_K, Tinf_K      the state
  cd_walker                           eq. (12)
  massImplied_amu                     m from s, U, T_inf -- reported so the reader can see it is
                                      atomic oxygen at the paper's conditions

Conditions follow [W14b] sect. IV: U = 7500 m/s, T_s = 300 K, T_inf = 1100 K, s ~ 7.0.

Run:  python dria_plate_check.py [--out PATH]
"""
import argparse
import math
import os

BOLTZMANN = 1.380649e-23          # J/K   (CODATA 2018, exact)
AMU = 1.66053906660e-27           # kg

U_MS = 7500.0                     # [W14b] sect. IV
TW_K = 300.0
TINF_K = 1100.0
SPEED_RATIOS = (5.0, 7.0, 9.0)
ALPHAS = (1.0, 0.95, 0.90, 0.85, 0.50, 0.0)

HEADER = ["caseId", "s", "alpha", "Tw_K", "Tinf_K", "cd_walker", "massImplied_amu"]


def implied_mass(s, u=U_MS, t_inf=TINF_K):
    """s = U sqrt(m / (2 k T_inf)) inverted for m [kg]."""
    return 2.0 * BOLTZMANN * t_inf * s * s / (u * u)


def cd_plate_walker(s, alpha, t_w=TW_K, t_inf=TINF_K, u=U_MS):
    """[W14b] eq. (12) with eqs. (5) and (10) substituted."""
    m = implied_mass(s, u, t_inf)
    t_ki = m * u * u / (3.0 * BOLTZMANN)
    t_kr = t_ki * (1.0 - alpha) + alpha * t_w
    return ((2.0 + 1.0 / (s * s)) * math.erf(s)
            + 2.0 / (math.sqrt(math.pi) * s) * math.exp(-s * s)
            + math.sqrt(math.pi) / s * math.sqrt(t_kr / t_inf))


def fmt(x):
    return repr(float(x))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(here, "dria_plate_oracle.csv"))
    args = ap.parse_args()
    rows = [",".join(HEADER)]
    case = 0
    for s in SPEED_RATIOS:
        for alpha in ALPHAS:
            case += 1
            rows.append(",".join([
                str(case), fmt(s), fmt(alpha), fmt(TW_K), fmt(TINF_K),
                fmt(cd_plate_walker(s, alpha)), fmt(implied_mass(s) / AMU)]))
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    print("wrote %d cases to %s" % (case, args.out))


if __name__ == "__main__":
    main()
