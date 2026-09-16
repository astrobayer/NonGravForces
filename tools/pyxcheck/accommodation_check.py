"""Independent oracle for two Langmuir-isotherm energy-accommodation models (pure Python).

Written from the source texts only; the MATLAB implementation was not opened.

Standard library only (math, csv, argparse, sys). Writes accommodation_oracle.csv.

---------------------------------------------------------------------------------------
MODEL 1 -- Pilinski et al. (2010), as reported in
  Pilinski, M. D. (2011), "Dynamic Gas-Surface Interaction Modeling for Satellite
  Aerodynamic Computations", PhD thesis, University of Colorado, Chapter 7.
  (Printed page numbers are the page headers; PDF page = printed page + 18.)

  Eq. (7.1), p. 129 (section 7.2):      alpha = K_L * P / (1 + K_L * P)

  P, p. 129: "a scaled static pressure component (static pressure divided by
  Boltzmann's constant)", and "The value of P was set equal to nO*Ta", i.e.
                                         P = n_O * T_a          [m^-3 K]
  with n_O the atomic-oxygen number density and T_a the atmospheric temperature,
  both from NRLMSISE-00 (p. 128). Boltzmann's constant is NOT in P.

  K_L, p. 130 (section 7.3), repeated p. 131:
                                         K_L = 4.98e-17  m^3 K^-1
  so K_L * P is dimensionless (m^3 K^-1 * m^-3 K). The thesis states K_L is tied
  to NRLMSISE-00 inputs and must be refitted for another atmosphere model (p. 131).

  Eq. (7.1) is eq. (6.9)-(6.10) (p. 113) with theta' = alpha, i.e. adsorbate
  accommodation 1 and substrate accommodation alpha_s = 0; p. 131 states the model
  "implicitly assumes a substrate accommodation coefficient of 0.0".

---------------------------------------------------------------------------------------
MODEL 2 -- Walker et al. (2014c) coverage with the Mehta et al. (2022) substrate term,
  as printed in Paul, Sheridan, Licata & Mehta (2022), arXiv:2210.08364v2.

  Eq. (10), p. 8:   C_D = f_sc C_D,ads + (1 - f_sc) C_D,surf      (not evaluated here)
  Eq. (11), p. 8:   f_sc = K_DRIA P_O / (1 + K_DRIA P_O)
                    K_DRIA = 1.44e6 ("Langmuir adsorbate constant for the DRIA model")
                    P_O    = "the partial pressure of atomic oxygen"
  Eq. (12), p. 8:   species weights chi_k m_k over k in {H, He, N, N2, O, O2}
  p. 13:            alpha_ads = 1
  Eq. (13), p. 13:  alpha_surf = 3 mu / (1 + mu)^2,
                    mu = (sum_{k=1..6} chi_k m_k) / m_surf
                    chi_k mole fraction, m_k species mass, m_surf the mass of a particle
                    of the surface lattice (= 263.3223 amu); mole fractions of
                    H, He, N, N2, O, O2 from NRLMSISE-00.

  UNIT ASSUMPTION (not printed in the paper): K_DRIA is taken in Pa^-1 and
                    P_O = n_O * k_B * T   [Pa],
  consistent with SI. (Under this reading K_DRIA * k_B = 1.988e-17 m^3 K^-1, the same
  order as Pilinski's 4.98e-17; a Torr^-1 reading would give f_sc = 0.13 for
  n_O = 1e15 m^-3, T = 1000 K, against 0.95 under the Pa^-1 reading.)

---------------------------------------------------------------------------------------
IMPLEMENTATION CHOICES (the sources leave these open; stated so the consumer can match)

  * Mole fractions are normalised over the six species of eq. (13) only:
        chi_k = n_k / (n_H + n_He + n_N + n_N2 + n_O + n_O2),
    so sum chi_k m_k is the mean mass of that six-species mixture (matches Pilinski
    eq. (6.8), p. 112). Ar and anomalous O are carried as inputs and IGNORED by both
    models, including the normalisation.
  * n_O in both P and P_O is the ordinary atomic-oxygen density; anomalous O is not added.
  * Species masses: m_k = M_k * 1e-3 / N_A [kg] from the molar masses M_k [g/mol];
    m_surf = surfaceMass_amu * u [kg]; mu is the ratio of the two in kg;
    meanMass_amu = (sum chi_k m_k) / u. N_A*u = 0.99999999965e-3 kg/mol, so
    meanMass_amu differs from the plain g/mol average by 3.5e-10 relative.
  * Constants: N_A = 6.02214076e23 mol^-1, k_B = 1.380649e-23 J/K,
    u = 1.66053906660e-27 kg.

OUTPUT COLUMNS
  caseId, He, O, N2, O2, Ar, H, N, Oanom      number densities [m^-3]
  T_K                                          temperature [K]
  surfaceMass_amu                              m_surf [amu]
  alpha_pilinski2010                           Langmuir with K = 7.50e-17 (2010 paper eq. 8)
  alpha_pilinski2011                           Langmuir with K = 4.98e-17 (2011 thesis eq. 7.1)
  coverage_walker2014                          f_sc, eq. (11) of Paul et al.
  alphaSubstrate_walker2014                    alpha_surf, eq. (13) of Paul et al.
  PO_Pa                                        n_O k_B T [Pa]
  meanMass_amu                                 six-species mean mass [amu]

Run:  python accommodation_check.py [--out PATH]
"""
import argparse
import csv
import math
import sys

# Physical constants (SI; values specified for this oracle)
AVOGADRO = 6.02214076e23          # mol^-1
BOLTZMANN = 1.380649e-23          # J K^-1
AMU = 1.66053906660e-27           # kg

# Input column order and molar masses [g/mol]
SPECIES = ("He", "O", "N2", "O2", "Ar", "H", "N", "Oanom")
MOLAR_MASS_G = {
    "He": 4.002602,
    "O": 15.9994,
    "N2": 28.0134,
    "O2": 31.9988,
    "Ar": 39.948,
    "H": 1.00794,
    "N": 14.0067,
    "Oanom": 15.9994,
}

# Species entering Paul et al. (2022) eqs. (12)-(13)
PAUL_SPECIES = ("H", "He", "N", "N2", "O", "O2")

# Model constants
# ONE MODEL, TWO PUBLISHED FITS. The 2010 paper states K in its Results section; the thesis
# refits it. Both are transcribed from their own source, which is what keeps this an oracle.
K_L_PILINSKI_2010 = 7.50e-17      # m^3 K^-1        (2010 paper, Results)
K_L_PILINSKI_2011 = 4.98e-17      # m^3 K^-1        (thesis p. 130)
K_DRIA = 1.44e6                   # Pa^-1 (assumed) (Paul et al. p. 8, eq. 11)

SURFACE_MASSES_AMU = (263.3223, 65.0)
N_CASES = 12

HEADER = ["caseId"] + list(SPECIES) + [
    "T_K", "surfaceMass_amu", "alpha_pilinski2010", "alpha_pilinski2011", "coverage_walker2014",
    "alphaSubstrate_walker2014", "PO_Pa", "meanMass_amu"]


def langmuir(kp):
    """Langmuir isotherm kp / (1 + kp) for a dimensionless kp >= 0."""
    return kp / (1.0 + kp)


def alpha_pilinski2010(n_o, t_k):
    """Pilinski, Argrow & Palo (2010) eq. 8 with P = n_O T and K = 7.50e-17 m^3/K (Results)."""
    return langmuir(K_L_PILINSKI_2010 * n_o * t_k)


def alpha_pilinski2011(n_o, t_k):
    """Thesis eq. (7.1) with P = n_O * T_a (p. 129) and K_L = 4.98e-17 m^3/K (p. 130)."""
    return langmuir(K_L_PILINSKI_2011 * n_o * t_k)


def partial_pressure_o(n_o, t_k):
    """P_O = n_O k_B T [Pa] (ideal gas; SI reading of Paul et al. p. 8)."""
    return n_o * BOLTZMANN * t_k


def coverage_walker2014(n_o, t_k):
    """Paul et al. eq. (11): f_sc = K_DRIA P_O / (1 + K_DRIA P_O)."""
    return langmuir(K_DRIA * partial_pressure_o(n_o, t_k))


def mean_mass_kg(dens):
    """sum_k chi_k m_k over H, He, N, N2, O, O2, chi normalised over those six [kg]."""
    total = 0.0
    weighted = 0.0
    for name in PAUL_SPECIES:
        n = dens[name]
        total += n
        weighted += n * (MOLAR_MASS_G[name] * 1e-3 / AVOGADRO)
    if total <= 0.0:
        raise ValueError("the six-species number density sum must be positive")
    return weighted / total


def alpha_substrate_mehta2022(dens, surface_mass_amu):
    """Paul et al. eq. (13): alpha_surf = 3 mu / (1 + mu)^2, mu = mean mass / m_surf."""
    mu = mean_mass_kg(dens) / (surface_mass_amu * AMU)
    return 3.0 * mu / ((1.0 + mu) * (1.0 + mu))


def log_space(lo, hi, i, n=N_CASES):
    """i-th of n log-spaced values from lo to hi inclusive."""
    a = math.log10(lo)
    b = math.log10(hi)
    return 10.0 ** (a + (b - a) * i / (n - 1))


def build_cases():
    """Twelve deterministic thermosphere-like cases, high altitude (i = 0) to low (i = 11)."""
    cases = []
    for i in range(N_CASES):
        dens = {
            "He": log_space(1e12, 1e14, i),
            "O": log_space(1e12, 1e17, i),
            "N2": log_space(1e10, 1e16, i),
            "O2": log_space(1e8, 1e15, i),
            "Ar": log_space(1e6, 1e9, i),       # small, ignored by both models
            "H": log_space(1e10, 1e12, i),
            "N": log_space(1e9, 1e12, i),
            "Oanom": log_space(1e9, 1e6, i),    # small, ignored by both models
        }
        t_k = 1500.0 - 900.0 * i / (N_CASES - 1)   # 1500 K down to 600 K
        m_s = SURFACE_MASSES_AMU[i % 2]           # 263.3223, 65, 263.3223, ...
        cases.append((i + 1, dens, t_k, m_s))
    return cases


def check_inputs(dens, t_k, m_s):
    for name in SPECIES:
        v = dens[name]
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError("number density of %s must be finite and >= 0" % name)
    if not (math.isfinite(t_k) and t_k > 0.0):
        raise ValueError("temperature must be finite and > 0")
    if not (math.isfinite(m_s) and m_s > 0.0):
        raise ValueError("surface mass must be finite and > 0")


def fmt(x):
    return "%.17g" % x


def main(argv=None):
    default_out = __file__.replace("\\", "/").rsplit("/", 1)
    default_out = (default_out[0] + "/" if len(default_out) == 2 else "") + "accommodation_oracle.csv"
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=default_out, help="output CSV path")
    args = parser.parse_args(argv)

    rows = []
    for case_id, dens, t_k, m_s in build_cases():
        check_inputs(dens, t_k, m_s)
        n_o = dens["O"]
        row = [str(case_id)] + [fmt(dens[s]) for s in SPECIES] + [
            fmt(t_k),
            fmt(m_s),
            fmt(alpha_pilinski2010(n_o, t_k)),
            fmt(alpha_pilinski2011(n_o, t_k)),
            fmt(coverage_walker2014(n_o, t_k)),
            fmt(alpha_substrate_mehta2022(dens, m_s)),
            fmt(partial_pressure_o(n_o, t_k)),
            fmt(mean_mass_kg(dens) / AMU),
        ]
        rows.append(row)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    sys.stdout.write("wrote %d cases to %s\n" % (len(rows), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
