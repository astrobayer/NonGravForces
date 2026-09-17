"""Independent oracle for three Langmuir-isotherm energy-accommodation models (pure Python).

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
MODEL 3 -- SESAM, Pilinski, Argrow, Palo & Bowman (2013), J. Spacecraft Rockets 50(3) 556,
  DOI 10.2514/1.A32348.

  Eq. (17), p. 561:   theta_0 = K_L(E) P_O / (1 + K_L(E) P_O)
  Eq. (18), p. 561:   alpha   = (1 - theta_0) alpha_s + theta_0
  Eq. (19), p. 561:   alpha_s = 2.4 mu / (1 + mu)^2      <- 2.4, where Paul et al. print 3.0
  Eq. (20), p. 561:   mu      = mbar / m_s,  m_s = 65 amu (p. 561 and Table C1)
  Eq. (22), p. 561:   K_L(E)  = s_o(E, Eb, Ta) K_L,o + K_L,f
  Table C1, p. 569:   Eb = 5.7 eV, K_L,o = 5e6 torr^-1, K_L,f = 3e4 torr^-1, m_s = 65 amu

  P_O IS NOT THE THERMODYNAMIC PARTIAL PRESSURE. Eq. (15), p. 560, and Appendix B, p. 569:
                      P_O = 1/2 rho_O V_r^2 C_D,in,O(s)
  where rho_O = n_O m_O and C_D,in,O is the free-molecular sphere drag coefficient with the
  re-emission term set to zero (Appendix B states exactly that construction):
                      C_D,in,O = (2s^2+1)/(sqrt(pi) s^3) exp(-s^2) + (4s^4+4s^2-1)/(2s^4) erf(s)
  with s the atomic-oxygen speed ratio V_r / sqrt(2 k_B T_a / m_O). The scan in hand renders
  the second numerator as "+1"; the classical Schaaf-Chambre sphere coefficient that
  Appendix B describes has "-1", which is what is implemented. At LEO s ~ 8 the two differ
  by 2/(2 s^4) = 2.4e-04 relative.

  THE STICKING COEFFICIENT IS COMPUTED HERE BY QUADRATURE OF ITS DEFINITION, NOT BY THE
  PAPER'S CLOSED FORM, so that this file checks the closed form rather than repeating it.
  Eq. (26a), p. 562, with the adsorption efficiency D(E) = 1 as the paper assumes:
                      s_o = int_0^Eb f(E) dE / int_0^inf f(E) dE
  and eq. (25), p. 562, whose exponent -Er + 2 sqrt(Er Ex) - Ex is -(sqrt(Er) - sqrt(Ex))^2:
                      f(E) propto exp(-(sqrt(Er) - sqrt(E))^2 / (k_B T_a))
  Substituting u = sqrt(E) removes the integrable square-root singularity at E = 0 and makes
  the integrand a Gaussian times a linear factor, which Simpson integrates to ~1e-15 here.
  The upper limit of the denominator is taken at sqrt(Er) + 12 sqrt(k_B T_a), where the
  Gaussian has fallen by exp(-144).

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

# SESAM, Table C1 (p. 569). Torr is exact by the definition of the standard atmosphere.
TORR_PA = 101325.0 / 760.0
EV_J = 1.602176634e-19            # exact since the 2019 SI redefinition
SESAM_EB_J = 5.7 * EV_J
SESAM_KLO_INV_PA = 5e6 / TORR_PA
SESAM_KLF_INV_PA = 3e4 / TORR_PA
SESAM_MS_AMU = 65.0
SESAM_GOODMAN = 2.4               # eq. (19); Paul et al. print 3.0 for the same law

SURFACE_MASSES_AMU = (263.3223, 65.0)
N_CASES = 12

HEADER = ["caseId"] + list(SPECIES) + [
    "T_K", "surfaceMass_amu", "relativeSpeed_ms",
    "alpha_pilinski2010", "alpha_pilinski2011", "coverage_walker2014",
    "alphaSubstrate_walker2014", "PO_Pa", "meanMass_amu",
    "sticking_sesam", "PO_sesam_Pa", "coverage_sesam", "alphaSubstrate_sesam"]


def langmuir(kp):
    """Langmuir isotherm kp / (1 + kp) for a dimensionless kp >= 0."""
    return kp / (1.0 + kp)


def alpha_pilinski2010(n_o, t_k):
    """Pilinski, Argrow & Palo (2010) eq. 8 with P = n_O T and K = 7.50e-17 m^3/K (Results)."""
    return langmuir(K_L_PILINSKI_2010 * n_o * t_k)


def alpha_pilinski2011(n_o, t_k):
    """Thesis eq. (7.1) with P = n_O * T_a (p. 129) and K_L = 4.98e-17 m^3/K (p. 130)."""
    return langmuir(K_L_PILINSKI_2011 * n_o * t_k)


def _simpson(f, lo, hi, n):
    """Composite Simpson on [lo, hi] with n panels (n even). Standard library only."""
    if n % 2:
        n += 1
    h = (hi - lo) / n
    total = f(lo) + f(hi)
    for i in range(1, n):
        total += (4.0 if i % 2 else 2.0) * f(lo + i * h)
    return total * h / 3.0


def sticking_sesam(er_j, t_k, eb_j=SESAM_EB_J, panels=20000):
    """Eq. (26a) BY QUADRATURE, so this file checks the paper's closed form rather than
    repeating it. With u = sqrt(E) the integrand of eq. (25) becomes 2 u exp(-(sqrt(Er)-u)^2/kT),
    which is smooth everywhere (the u = 0 singularity of dE/du is removed by the Jacobian)."""
    kt = BOLTZMANN * t_k
    ur = math.sqrt(er_j)

    def integrand(u):
        d = ur - u
        return 2.0 * u * math.exp(-d * d / kt)

    top = _simpson(integrand, 0.0, math.sqrt(eb_j), panels)
    # exp(-144) below the peak is 1e-63; the tail beyond contributes nothing at double precision
    bottom = _simpson(integrand, 0.0, ur + 12.0 * math.sqrt(kt), panels)
    return top / bottom


def po_sesam(n_o, t_k, v_r):
    """Eq. (15) with Appendix B's C_D,in,O: a RAM pressure, not n_O k_B T."""
    m_o = MOLAR_MASS_G["O"] * 1e-3 / AVOGADRO
    s = v_r / math.sqrt(2.0 * BOLTZMANN * t_k / m_o)
    cd_in = ((2.0 * s * s + 1.0) / (math.sqrt(math.pi) * s ** 3) * math.exp(-s * s)
             + (4.0 * s ** 4 + 4.0 * s * s - 1.0) / (2.0 * s ** 4) * math.erf(s))
    return 0.5 * (n_o * m_o) * v_r * v_r * cd_in


def sesam(n_o, t_k, v_r, mbar_kg):
    """Eqs (15)-(22) and (26a). Returns (s_o, P_O, theta_0, alpha_s)."""
    m_o = MOLAR_MASS_G["O"] * 1e-3 / AVOGADRO
    er = 0.5 * m_o * v_r * v_r                               # eq. (15) text: Er = 1/2 m_O V_r^2
    s_o = sticking_sesam(er, t_k)                            # eq. (26a), by quadrature
    k_l = s_o * SESAM_KLO_INV_PA + SESAM_KLF_INV_PA          # eq. (22)
    p_o = po_sesam(n_o, t_k, v_r)                            # eq. (15)
    theta = langmuir(k_l * p_o)                              # eq. (17)
    mu = mbar_kg / (SESAM_MS_AMU * AMU)                      # eq. (20)
    alpha_s = SESAM_GOODMAN * mu / (1.0 + mu) ** 2           # eq. (19)
    return s_o, p_o, theta, alpha_s


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
        # SESAM's sticking coefficient is ~1 below 5 km/s and ~0 above 11 km/s, so the cases
        # sweep 3 to 11.8 km/s: a set that stayed at 7.6 would check a near-constant and would
        # pass with the energy dependence -- the whole point of the model -- deleted.
        v_r = 3000.0 + 800.0 * i
        cases.append((i + 1, dens, t_k, m_s, v_r))
    return cases


def check_inputs(dens, t_k, m_s, v_r):
    for name in SPECIES:
        v = dens[name]
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError("number density of %s must be finite and >= 0" % name)
    if not (math.isfinite(t_k) and t_k > 0.0):
        raise ValueError("temperature must be finite and > 0")
    if not (math.isfinite(m_s) and m_s > 0.0):
        raise ValueError("surface mass must be finite and > 0")
    if not (math.isfinite(v_r) and v_r > 0.0):
        raise ValueError("relative speed must be finite and > 0")


def fmt(x):
    return "%.17g" % x


def main(argv=None):
    default_out = __file__.replace("\\", "/").rsplit("/", 1)
    default_out = (default_out[0] + "/" if len(default_out) == 2 else "") + "accommodation_oracle.csv"
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=default_out, help="output CSV path")
    args = parser.parse_args(argv)

    rows = []
    for case_id, dens, t_k, m_s, v_r in build_cases():
        check_inputs(dens, t_k, m_s, v_r)
        n_o = dens["O"]
        s_o, po_s, theta_s, alpha_s_sesam = sesam(n_o, t_k, v_r, mean_mass_kg(dens))
        row = [str(case_id)] + [fmt(dens[s]) for s in SPECIES] + [
            fmt(t_k),
            fmt(m_s),
            fmt(v_r),
            fmt(alpha_pilinski2010(n_o, t_k)),
            fmt(alpha_pilinski2011(n_o, t_k)),
            fmt(coverage_walker2014(n_o, t_k)),
            fmt(alpha_substrate_mehta2022(dens, m_s)),
            fmt(partial_pressure_o(n_o, t_k)),
            fmt(mean_mass_kg(dens) / AMU),
            fmt(s_o),
            fmt(po_s),
            fmt(theta_s),
            fmt(alpha_s_sesam),
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
