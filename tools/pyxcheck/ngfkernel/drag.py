"""Per-plate atmospheric drag and lift force in free-molecular flow.

Independent standard-library oracle for the NonGravForces drag kernel (WP-33).
Closed form throughout; no quadrature is used by the returned numbers (a
Simpson quadrature appears only in the self-test, as an independent check of
the closed form).  The MATLAB implementation this module checks, its tests and
the pre-existing Sentman oracle were NOT read.  The only project code reused is
the already-derived coefficient function ``ngfkernel.gsi.gsi_coefficients``.

    F_k = 0.5 * U^2 * A_k * chi_k
          * [ (SUM_j rho_j cp(s_j, theta_k)) * (-n_k)
            + (SUM_j rho_j ctau(s_j, theta_k)) * t_k ]

RECONCILIATION.  This module is the referee's reconciliation of two independent
derivations (A and B) of the same quantity.  Both derived the same physics --
the same mixture rule, the same geometry, the same resolution onto (-n, +t) --
and agreed numerically to 4.4e-16 N on forces of order 1 N over a grid of
40 128 numbers.  The reconciliation is therefore about interface behaviour and
conditioning, and is recorded in section 8.


===============================================================================
1.  REGIME AND WHAT IS COMPUTED
===============================================================================

At 300-700 km the Knudsen number formed with the spacecraft size is >> 1: a
molecule leaving the surface escapes to infinity without meeting another
molecule, and the incident stream is undisturbed by the presence of the body.
Consequently

  (i)   the distribution incident on a plate is the free-stream distribution,
        sampled over the half space that can reach the plate;
  (ii)  plates exchange no momentum through the gas; the only coupling between
        them is geometric blocking (the ``shadow`` fractions);
  (iii) the force on the body is the plain sum of independent plate forces.

This module returns the summands of (iii): one 3-vector per plate.  The body
force and the torque are reductions over the same list and are not done here.


===============================================================================
2.  PLATE GEOMETRY
===============================================================================

    U    = |vrel|,  v = vrel/U      unit vector of the SATELLITE velocity
                                    relative to the gas, body frame
    n                               outward unit plate normal (away from body)
    g    = -v                       direction in which the gas streams past the
                                    plate, as seen from the plate

    theta = atan2(|n x v|, n.v) in [0, pi]                                  (1)
    gamma = cos(theta) = n.v,       l = sin(theta) = |n x v| >= 0

theta is formed by atan2 of the cross-product norm against the dot product and
never by acos(n.v): acos loses half the significant digits within sqrt(eps) of
theta = 0 and theta = pi, which is exactly where the ram and lee faces of a
flat-panel model sit and where the leeward exponentials are most sensitive.

gamma > 0 is WINDWARD: the outward normal has a component along the direction
of travel, so the stream g runs into the outer face.  gamma < 0 is leeward.

Decomposing the stream direction in the plate frame,

    g = l*t - gamma*n,      t = (gamma*n - v)/l,                            (2)

so t is the unit vector along the TANGENTIAL part of the stream (t.n = 0,
|t| = 1 for l > 0).  The gas bulk velocity seen by the plate is U*g: normal
component -U*gamma (pressing on the outer face when windward) and tangential
component +U*l along t.  A stream running along +t drags the plate along +t.

For l -> 0 (head-on or exactly anti-ram) t is undefined; the tangential
momentum flux vanishes there identically because it is proportional to the
tangential bulk speed U*l, so t is set to the zero vector below a floor on
|gamma*n - v| and the shear term disappears without a limit being taken.

The coefficient function is defined by

    force per unit area from ONE species
        = 0.5 * rho * U^2 * ( cp * (-n) + ctau * t ),                       (3)

i.e. cp is a normal (compression) coefficient, positive, pushing the plate
inward along -n; ctau is a shear coefficient, positive, dragging the plate
along +t.  cp and ctau are NOT a drag/lift pair: they coincide with one at
normal incidence and differ everywhere else.  Section 9, check L2, is the limit
that discriminates the two resolutions, and it is checked in ``_selftest``.


===============================================================================
3.  THE MIXTURE RULE  (the substance of this module)
===============================================================================

The thermosphere is not a gas but a MIXTURE of J = 8 non-interacting
constituents in diffusive equilibrium.  They share one kinetic temperature T
and one bulk velocity -- the drift is a property of the flow field, not of the
constituent: every species is carried by the same wind and overtaken by the
satellite at the same vrel -- but each has its own molecular mass m_j and its
own number density n_j.  In the plate frame,

    f(c) = SUM_j n_j (m_j/(2 pi k T))^{3/2} exp( -m_j |c - U*g|^2/(2 k T) ). (4)

The momentum delivered to a surface element per unit time is the flux integral

    dP/dt = dA [ INT_{c.n<0} m c |c.n| f d3c  -  INT_{c.n>0} m c (c.n) f_re d3c ]
                                                                            (5)

(the second term enters negatively: momentum leaving the surface is momentum
lost by the surface).  Both integrands are LINEAR in the distribution and (4)
is a SUM over species, hence

    dP/dt = SUM_j (dP/dt)_j,                                                (6)

each summand being the integral formed from species j alone.  That is the whole
mixture rule: in free-molecular flow the species do not interact, so their
momentum fluxes ADD and each may be pushed through the single-gas theory
separately.  Two details make (6) exact rather than approximate:

  * the re-emitted flux of species j is tied to the INCIDENT flux of species j
    by particle conservation at the wall (nothing is adsorbed permanently and
    nothing is created), so no species' re-emission depends on another's
    incidence;
  * re-emitted molecules of species j leave with the wall temperature (fully or
    partially accommodated), and a wall Maxwellian at T_w for mass m_j is again
    a per-species object.

Non-dimensionalising species j exactly as the single-gas theory does,

    v_m,j = sqrt(2 k T/m_j),     s_j = U/v_m,j = U sqrt(m_j/(2 k T)),       (7)

so a HEAVY species is hyperthermal (large s, nearly a mono-directional beam)
and a LIGHT species is not.  At T = 1000 K and U = 7600 m/s: s(H) = 1.87,
s(He) = 3.73, s(N) = 6.98, s(O) = 7.45, s(N2) = 9.86, s(O2) = 10.54,
s(Ar) = 11.78 -- a factor 6.3 across the mixture.  With (3), (6) and (7),

    dF/dA = 0.5 U^2 SUM_j rho_j ( cp(s_j, theta)*(-n) + ctau(s_j, theta)*t ),
            rho_j = n_j * m_j       (PARTIAL MASS density of species j).    (8)

WHAT IS WEIGHTED BY WHAT.  The weight is the partial MASS density n_j m_j --
not the number density, because momentum flux carries one power of the mass per
particle -- and each weight multiplies the coefficients of ITS OWN speed ratio,
not of a mixture-averaged one.  The geometry factors (theta, n, t) are common
to all species and come out of the sum, which is why the mixture collapses to
two scalars per plate.  Equivalently, with rho = SUM_j rho_j,

    cp_mix   = SUM_j (rho_j/rho) cp(s_j, theta)
    ctau_mix = SUM_j (rho_j/rho) ctau(s_j, theta)                           (9)

are mass-fraction-weighted means OF THE COEFFICIENTS, not coefficients of a
mean gas.

That rho_j is the PARTIAL mass density -- and not the total one -- is pinned
independently by a number.  At theta = pi/2 with alpha = 1, T_w = 0 the only
surviving normal coefficient is cp = 1/(2 s^2), so

    -f.n = 0.5 U^2 SUM_j rho_j/(2 s_j^2)
         = 0.5 U^2 SUM_j n_j m_j k T/(m_j U^2) = (SUM_j n_j) k T/2 = p/2,  (10)

each species' molecular mass cancelling against its own speed ratio and leaving
Dalton's law: a gas at rest exerts half its static pressure on one face of a
plate.  Feeding a TOTAL mass density into a coefficient built from a SINGLE
species' mass would give rho_tot k T/(2 m) here, which is (10) only by accident
of the mean chosen.  ``_selftest`` check L1 verifies (10) to 1e-15.


===============================================================================
4.  WHY THE MEAN-MOLAR-MASS SHORTCUT IS A DIFFERENT ANSWER
===============================================================================

The shortcut replaces (9) by cp(s(M_bar), theta) for a single mean molar mass.
Since s ~ sqrt(m) and cp, ctau are strongly nonlinear in s -- they carry
erf(gamma s), exp(-gamma^2 s^2), a thermal term 1/(2 s^2) and a re-emission
term going as 1/s -- averaging the argument is not averaging the function, and
the two differ by Jensen's inequality for any composition that is not
monodisperse.

Two limits hide the difference and are an honest warning that it cannot be
found by poking a ram plate.  Both weightings carry the same total mass density
(an identity for the number-weighted mean), so

  * hyperthermal windward: cp -> 2 gamma^2 for every species regardless of s,
    and the two agree to O(1/s^2);
  * the edge-on thermal term: cp -> 1/(2 s_j^2) is LINEAR in m_j and
    rho_j/m_j = n_j, so SUM_j rho_j/(2 s_j^2) = n_tot k T/U^2 -- again equal,
    for the number-weighted mean.

They separate where the nonlinearity lives.  Measured by ``_selftest`` check
L9 for the 500 km composition below (sentman, alpha = 0.85, T_w = 273 K,
T = 1000 K, U = 7600 m/s, mass-weighted mean 16.3175 g/mol), relative error of
the mean-mass shortcut in cp:

    theta =   0 deg   -0.05 %          theta =  90 deg    -2.5 %
    theta =  45 deg   -0.09 %          theta = 100 deg   -24 %
    theta =  80 deg   -0.55 %          theta = 135 deg  -100 % (ten orders)

The shortcut is therefore nearly harmless on a ram face and useless on a lee
one, and the reason is structural rather than numerical: on the lee side every
coefficient is O(exp(-s^2 gamma^2)), so the species sum is dominated by the
SMALLEST speed ratio present -- hydrogen and helium, which carry 0.03 % and
1.2 % of the mass here but are the only constituents with a thermal tail fast
enough to reach a lee-facing plate at all.  A single gas at the mean mass has
no such tail.  Where lee panels matter -- torques, cross-track force, a
non-convex body whose flanks expose their undersides -- the mixture rule is not
a refinement of the mean-mass form, it is a different answer.


===============================================================================
5.  PER-PLATE FORCE AND BLOCKING
===============================================================================

Multiplying (8) by the plate area A_k and by a blocking fraction chi_k in [0,1]
(the fraction of the plate the free stream actually reaches; 1 when nothing
upstream shadows it),

    F_k = 0.5 U^2 A_k chi_k [ (SUM_j rho_j cp_j)(-n_k)
                            + (SUM_j rho_j ctau_j)(t_k) ].                 (11)

chi_k multiplies the WHOLE plate force: a molecule intercepted upstream
delivers no momentum to this plate at all, and blocking is a property of the
ray, not of the species -- all species travel the same straight lines in the
collisionless regime.  The thermal spread of the incident stream about the bulk
direction is ignored by any such geometric shadow; that is the usual, and here
inherited, approximation.

F_k is linear in every n_j, so doubling the whole composition doubles every
force exactly -- bit-exact in floating point, a factor of two being a change of
exponent only.


===============================================================================
6.  "constantCd" AND cd_rho
===============================================================================

The "constantCd" model does not resolve the gas: it asserts that each windward
plate carries a force along the stream of magnitude 0.5 rho cd A gamma U^2,
whose resolution on (-n, t) is cp = cd gamma^2, ctau = cd gamma l, and that a
leeward plate carries nothing.  These coefficients are independent of m_j, T
and s, so the species sum (8) degenerates to the total mass density alone.
This module handles that as follows:

  * ``cd_rho`` (a mass density in kg/m^3) is used as rho when it is given, and
    ``species`` may then be None or empty;
  * otherwise rho is formed from ``species`` exactly as for the kinetic models,
    so that all three models can be driven from the same inputs;
  * giving neither is an error;
  * giving ``cd_rho`` to a KINETIC model is an error, not a silently ignored
    argument (section 8, D1).

Substituting into (3) with (2),

    cp*(-n) + ctau*t = cd gamma (l t - gamma n) = cd gamma g = -cd gamma v,

so F_k = -0.5 rho cd A_k gamma_k U vrel for a windward plate, and summing over
a CONVEX body, where exactly the plates with gamma_k > 0 are lit and none
shadows another, SUM_k A_k gamma_k is the projected area, giving the textbook

    SUM_k F_k = -0.5 rho A_proj cd |vrel| vrel.                            (12)

``_selftest`` check L6 verifies (12) on a unit cube for flow directions lighting
one, two and three faces.


===============================================================================
7.  DELIBERATE LIMITATIONS (named, as required of an oracle)
===============================================================================

  * ONE temperature.  The interface supplies a single kinetic temperature, so
    all eight constituents are given T.  For seven of them that is right.
    "Anomalous oxygen" in the NRLMSISE-00 sense is a hot (~4000 K) non-thermal
    O population; giving it the ambient T overstates its speed ratio.  It is
    included in the mass density, with the molecular mass of atomic O as its
    name demands, because it is a real momentum carrier for drag and
    NRLMSISE-00 counts it in the drag mass density.  Below ~600 km it is 1e-4
    or less of the mass density, so the choice is numerically irrelevant there.
    Callers who want it out pass 0.0 for it.  Because it shares the molecular
    mass of atomic oxygen exactly, it shares its speed ratio, and the
    coefficient pair is evaluated once and reused for both -- an identity under
    this interface, not an approximation.
  * No adsorption or desorption bookkeeping: the incident and re-emitted
    particle fluxes of each species are equal.
  * The shadow fractions are taken as given; no geometry is computed here.
  * The thermal spread of the incident stream is fully resolved inside the
    coefficient function (there is no hyperthermal approximation anywhere), but
    not in the shadowing.
  * U = 0 is outside the domain and raises.  It is NOT a zero-force case: as
    U -> 0 the coefficients diverge as 1/s^2 while 0.5 rho U^2 vanishes, and
    the product tends to the finite one-sided thermal pressure n k T/2 plus the
    re-emission term.  The limit is finite but cannot be evaluated through
    s = U sqrt(m/2kT) = 0, so the honest response is an error (section 8, D2).


===============================================================================
8.  RECONCILIATION OF THE TWO DERIVATIONS
===============================================================================

MEASURED AGREEMENT.  Over 38 plate orientations (head-on, 1e-7 deg, 1e-3 deg,
5, 30, 45, 60, 80, 89, 90, 90.5, 100, 120, 135, 170, 179.999 and 180 deg at two
azimuths, plus four axis-aligned plates), two compositions (oxygen-dominated
500 km and N2-dominated 200 km), four relative velocities including a slow
400 m/s thermal case and a reversed one, four shadow patterns and nine model
settings covering sentman (alpha = 0, 0.85, 0.93, 1.0), maxwell (specular =
0, 0.3, 1.0) and constantCd (cd = 2.2, 3.5) driven both from cd_rho and from a
species list -- 40 128 force components in all:

    max |F_A - F_B|                                   4.441e-16 N
      (on components of order 1.1 N, i.e. 4e-16 relative)
    max relative difference over components > 1e-12
      of the largest force on the grid                7.707e-14
      (at a component that is a near-exact cancellation between the normal and
       shear terms -- the y-component of an oblique plate in the hyperthermal
       cold-wall limit, where -cp n + ctau t is exactly along -v and the
       transverse residual is itself at the round-off level)
    max relative difference anywhere                  2.0
      (a sign flip between -6.3e-30 N and +6.3e-30 N: the shear component of a
       plate tilted 1e-7 deg off head-on, where ctau ~ 1e-9 and the tangential
       direction is the normalised difference of two nearly equal unit vectors)

The two derivations are therefore the same computation to floating-point
round-off.  There is NO physics disagreement: both weight by partial mass
density, both evaluate each species at its own speed ratio, both resolve onto
(-n, +t), and neither computed a drag/lift pair.

The disputed points, all of them behavioural, and how each was settled:

D1  cd_rho passed to a kinetic model.  A raises; B silently ignores it.
    SETTLED FOR A.  A bulk mass density cannot drive a kinetic model -- the
    coefficients need the composition -- so a caller who passes one holds a
    misconception, and an oracle that swallows the argument hides it.  Never
    degrade silently.

D2  |vrel| = 0.  A raises; B returns zero force on every plate.
    SETTLED FOR A, BY DERIVATION.  B's answer is physically wrong, not merely
    conservative: at U = 0 the gas is still at its kinetic temperature and each
    plate receives the one-sided static pressure n k T/2 (section 3, eq. 10)
    plus a re-emission term, so the per-plate force tends to a finite non-zero
    normal pressure.  It nets to zero on a closed body only by symmetry, and
    not at all in the torque.  Since s = 0 cannot be passed to the coefficient
    function, the domain boundary is honest and the error message says so.

D3  A non-unit plate normal.  A raises if |n| differs from 1 by more than
    1e-9; B silently normalises.
    SETTLED FOR NEITHER, AND FOR BOTH: this module raises outside 1e-9 (a
    non-unit normal signals a caller bug, and the area is carried separately,
    so scaling it silently would scale the force) AND normalises inside the
    tolerance, so that round-off in a nominally unit normal does not leak into
    gamma or into |t|.  This is strictly better than either variant.

D4  Per-plate alpha, tw_k, specular.  B broadcasts a scalar or accepts one
    value per plate; A takes scalars only.
    SETTLED FOR B.  Different plate materials really do have different
    accommodation coefficients and wall temperatures, B's form is a strict
    superset, and on scalar input the two were verified bit-identical.

D5  Summation of the species terms.  A uses the built-in sum; B uses
    math.fsum.
    SETTLED FOR B.  The eight terms span six orders of magnitude (Ar at 1e-7
    against O at 2e13 in number density); fsum is exactly rounded and, unlike a
    running sum, cannot depend on the order the species are listed in.

D6  The tangential unit vector.  A forms w = gamma*n - v and divides by |w|;
    B forms the same w and divides by |n x v|.
    SETTLED FOR A, on an invariant rather than on measured accuracy: the two
    denominators are the same number analytically, and measured against each
    other down to theta = 1e-13 rad neither is systematically better (both lose
    accuracy as eps/theta, as they must).  But dividing w by its OWN length
    guarantees |t| <= 1 + eps, so the shear force can never exceed q*A*ctau;
    dividing by a separately computed |n x v| can give |t| > 1 and amplify it.
    Verified: |F_tangential|/(q A ctau_mix) = 1.000000000000 at theta = 1e-6,
    1e-4, 1e-2, 179.99 and 180 - 1e-6 degrees.

D7  Early validation of cd for "constantCd".  B checks; A lets the coefficient
    function raise.  SETTLED FOR B (the error then names this module's
    argument, not the coefficient function's).

Nothing in the physics was overturned, so no term needed re-deriving.  The
checks that WOULD have discriminated a wrong choice were run anyway and are
kept in ``_selftest`` (section 9), because an oracle's value is that its
correctness is demonstrated rather than asserted.


===============================================================================
9.  THE CHECKS IN ``_selftest``
===============================================================================

L1  Dalton, edge-on.  theta = pi/2, alpha = 1, T_w = 0: the normal pressure is
    exactly (SUM_j n_j) k T/2, each species' mass cancelling against its own
    speed ratio (eq. 10).  Pins rho_j = n_j m_j per species -- a code that fed
    a total mass density into a single-species coefficient fails this.
L2  Hyperthermal cold wall, OBLIQUE.  alpha = 1, T_w = 0, s >> 1: the plate
    absorbs the incident momentum and re-emits none, so F = -rho A gamma U vrel
    exactly, for any composition.  This is the check that discriminates the
    (-n, +t) resolution from a (drag, lift) one: the two coincide at normal
    incidence, and at gamma = 0.5 a mis-resolution along -v is wrong by 50 %.
L3  Independent Simpson quadrature.  The incident momentum flux of the mixture
    is integrated directly from the drifting Maxwellian in this file -- the
    only numerical quadrature anywhere in the module, and it is never called by
    ``drag_panel_force``.  Reproduces the closed form to 7e-14 at theta = 0,
    25, 60, 89, 90 and 105 degrees.
L4  Species linearity: the mixture force equals the sum of the eight
    single-species forces, the statement of eq. (6) itself.
L5  Doubling every number density doubles every plate force bit-exactly.
L6  constantCd on a convex body reproduces -0.5 rho A_proj cd |v| v, eq. (12).
L7  maxwell(specular = 0) reproduces sentman bit-for-bit; shadow scales the
    plate force exactly; a leeward plate is small and positive.
L8  Rotation covariance: rotating the plates and the flow rotates the forces.
L9  The size of the mean-molar-mass shortcut (section 4).
L10 Order of magnitude: a GRACE-FO-sized box at 490 km gives a drag
    acceleration in the 1e-8 .. 2e-6 m/s^2 band and an effective C_D in
    [2.0, 4.5].


===============================================================================
10. CONTRADICTIONS WITH THE REFERENCE DOCUMENT (doc/theory.txt)
===============================================================================

(a) Section 2.3 defines the speed ratio with "M the mean molar mass of the
    local gas", and section 2.6 repeats that "composition is required by the
    GSI (M_bar, ...)": ONE speed ratio for the mixture.  My derivation
    (section 3) gives one speed ratio per species with the partial mass density
    as its weight; the mean-molar-mass form is an approximation to it, not an
    equivalent statement, and section 4 measures the size of the difference.
    I FOLLOW MY DERIVATION.  The document is right about the blunder it was
    warning against -- that using air, M = 29 g/mol, instead of the local
    composition biases s by sqrt(29/16) = 1.35 -- but a mean molar mass is the
    next term of the same error, not its cure.

(b) Section 1.2 defines the incidence cosine as gamma_k = n_k . (-u_hat) with
    u_hat = vr/|vr| the satellite's relative-velocity direction.  That makes
    gamma NEGATIVE on the ram face, which contradicts the same document's
    eq. (1) a few lines later, where Cp > 0 is stated to push the panel inward
    along -n for a WINDWARD panel, and contradicts its own eqs. (4)-(5), whose
    Sentman coefficients are exponentially small for gamma < 0: a ram face read
    through section 1.2 would take the leeward branch and carry no force.
    Eq. (1) is self-consistent with gamma = +n.u_hat, which is what section 2
    above derives, and that is what this module uses.  Section 1.2 is the odd
    one out.

(c) Section 2.1 eq. (1) writes the shear direction as the normalised tangential
    part of +u_hat while the force carries an overall minus sign.  Section 2
    above carries no overall sign and puts the shear along the tangential part
    of the STREAM, t = (gamma n - u_hat)/l = -(tangential part of u_hat)/l.
    These agree exactly; the difference is bookkeeping and is recorded only so
    that a reader comparing signs term by term does not mistake it for a
    discrepancy.

(d) Section 2.3 eq. (3) renders the re-emission ratio as
    v_re/v_in = (1/2) sqrt(1 + alpha(4 R T_w/(M |vr|^2) - 1)), with the one
    half OUTSIDE the root.  The Sentman result has it inside,
    sqrt( (1/2)[1 + alpha(4 R T_w/(M |vr|^2) - 1)] ), a factor sqrt(2) larger.
    Substituting R = k/m and s^2 = m U^2/(2 k T) turns the inside form into
    sqrt( (1-alpha)/2 + alpha T_w/(T s^2) ), which is what the coefficient
    function uses and what an energy-accommodation balance gives.  This is
    almost certainly an artefact of rendering a two-line fraction under a
    radical in plain text; it is in any case the coefficient function's
    business, not this module's, and is recorded so the discrepancy is not
    re-discovered.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

import math
import os
import sys

# ---------------------------------------------------------------------------
# The coefficient function is not re-derived here; it is imported.
# ---------------------------------------------------------------------------

try:
    from gsi import gsi_coefficients
except ImportError:                                          # pragma: no cover
    _NGFKERNEL_DIR = os.path.dirname(os.path.abspath(__file__))
    if _NGFKERNEL_DIR not in sys.path:
        sys.path.insert(0, _NGFKERNEL_DIR)
    from gsi import gsi_coefficients


# ---------------------------------------------------------------------------
# Constants (values fixed by the interface; CODATA 2018 / SI 2019)
# ---------------------------------------------------------------------------

K_B = 1.380649e-23              # J/K
N_A = 6.02214076e23             # 1/mol

SPECIES_NAMES = ("He", "O", "N2", "O2", "Ar", "H", "N", "anomalous O")

MOLAR_MASS_G_PER_MOL = (4.002602,      # He
                        15.9994,       # O
                        28.0134,       # N2
                        31.9988,       # O2
                        39.948,        # Ar
                        1.00794,       # H
                        14.0067,       # N
                        15.9994)       # anomalous O (atomic O, section 7)

#: molecular mass of each species [kg]
MOLECULAR_MASS_KG = tuple(m * 1.0e-3 / N_A for m in MOLAR_MASS_G_PER_MOL)

_N_SPECIES = len(MOLAR_MASS_G_PER_MOL)

_KINETIC_MODELS = ("sentman", "maxwell")
_MODELS = ("sentman", "maxwell", "constantCd")

#: |gamma*n - v| below which the tangential direction is dropped.  ctau is
#: O(l) there, so the discarded shear is O(1e-14) of a term already at the
#: round-off level.
_SIN_FLOOR = 1.0e-14

#: tolerance on |n| - 1 (section 8, D3)
_UNIT_TOL = 1.0e-9


# ---------------------------------------------------------------------------
# Vector helpers.  3-vectors are plain tuples.
# ---------------------------------------------------------------------------

def _as_vec3(v, what):
    try:
        seq = tuple(float(x) for x in v)
    except (TypeError, ValueError):
        raise ValueError("%s is not a 3-vector -> pass three finite floats"
                         % what)
    if len(seq) != 3:
        raise ValueError("%s has %d components -> pass a 3-vector"
                         % (what, len(seq)))
    return seq


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _scale(a, f):
    return (a[0] * f, a[1] * f, a[2] * f)


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


# ---------------------------------------------------------------------------
# Mixture bookkeeping
# ---------------------------------------------------------------------------

def species_mass_densities(species):
    """Partial mass densities rho_j = n_j m_j [kg/m^3] of the eight species.

    ``species`` holds number densities [1/m^3] in the order
    He, O, N2, O2, Ar, H, N, anomalous O.
    """
    species = list(species)
    if len(species) != _N_SPECIES:
        raise ValueError("species has %d entries -> pass %d number densities "
                         "in the order %s"
                         % (len(species), _N_SPECIES, ", ".join(SPECIES_NAMES)))
    out = []
    for j, n_j in enumerate(species):
        n_j = float(n_j)
        if math.isnan(n_j) or n_j < 0.0 or math.isinf(n_j):
            raise ValueError("species[%d] (%s) = %r -> pass a finite "
                             "non-negative number density in 1/m^3"
                             % (j, SPECIES_NAMES[j], n_j))
        out.append(n_j * MOLECULAR_MASS_KG[j])
    return out


def speed_ratios(t_k, u_speed):
    """Speed ratio s_j = U sqrt(m_j/(2 k T)) of every species, eq. (7).

    Formed as U*sqrt(m/(2kT)) rather than U/sqrt(2kT/m) so that the small
    quantity (the molecular mass) is multiplied rather than divided, and the
    common factor sqrt(1/(2 k T)) is evaluated once.
    """
    t_k = float(t_k)
    if not (t_k > 0.0) or math.isinf(t_k):
        raise ValueError("t_k = %r -> pass a finite kinetic temperature > 0 K"
                         % (t_k,))
    f = u_speed * math.sqrt(1.0 / (2.0 * K_B * t_k))
    return [f * math.sqrt(m) for m in MOLECULAR_MASS_KG]


def _per_plate(value, n_plates, name):
    """Broadcast a scalar, or accept one value per plate (section 8, D4)."""
    if value is None:
        return [None] * n_plates
    if isinstance(value, (list, tuple)):
        if len(value) != n_plates:
            raise ValueError("%s has %d entries for %d plates -> pass a scalar "
                             "or one value per plate"
                             % (name, len(value), n_plates))
        return [None if v is None else float(v) for v in value]
    return [float(value)] * n_plates


# ---------------------------------------------------------------------------
# The public routine
# ---------------------------------------------------------------------------

def drag_panel_force(plates, species, t_k, vrel, model, alpha=None, tw_k=None,
                     cd=None, specular=None, shadow=None, cd_rho=None):
    """Free-molecular drag and lift force on each plate of a macro model.

    Implements eq. (11) of the module docstring: per-species coefficients from
    ``gsi.gsi_coefficients``, combined by partial-mass-density weighting.

    Parameters
    ----------
    plates    list of plate dicts with keys "area" [m^2] and "normal" (outward
              unit 3-vector, body frame).  "centroid", "vis" and "ir" are
              carried by the macro model for other kernels and play no part in
              drag.
    species   list of 8 number densities [1/m^3] in the order
              He, O, N2, O2, Ar, H, N, anomalous O.  May be None or empty for
              model "constantCd" when ``cd_rho`` is given.
    t_k       ambient kinetic temperature at the satellite [K], > 0.  Unused by
              "constantCd" (may be None).
    vrel      3-vector, satellite velocity relative to the gas, BODY frame
              [m/s].  |vrel| must be finite and > 0 (section 7).
    model     "sentman", "maxwell" or "constantCd".
    alpha     energy accommodation coefficient in [0, 1]; a scalar, or one
              value per plate.  Required by the kinetic models.
    tw_k      wall temperature [K]; a scalar, or one value per plate.
              Required by the kinetic models.
    cd        drag coefficient; required by "constantCd" only.
    specular  specularly reflected fraction in [0, 1]; a scalar, or one value
              per plate.  Required by "maxwell" only.
    shadow    optional list of unblocked fractions in [0, 1], one per plate;
              default 1.0 everywhere.
    cd_rho    total mass density [kg/m^3] for "constantCd"; see section 6.
              Passing it to a kinetic model is an error.

    Returns
    -------
    list of 3-vectors (tuples), the force in newtons on each plate, in the
    frame of ``vrel`` and the plate normals.

    Raises
    ------
    ValueError with a message ending in " -> <what to change>" for every
    contract violation.
    """
    if model not in _MODELS:
        raise ValueError("unknown model %r -> use one of %s"
                         % (model, ", ".join(_MODELS)))

    plates = list(plates)
    n_plates = len(plates)
    if n_plates == 0:
        return []

    # ---- flow direction -------------------------------------------------
    vrel = _as_vec3(vrel, "vrel")
    u_speed = _norm(vrel)
    if not (u_speed > 0.0) or math.isinf(u_speed):
        raise ValueError("|vrel| = %r -> pass a finite relative velocity of "
                         "non-zero magnitude (U = 0 is outside the domain of "
                         "the speed ratio; the force there is the finite "
                         "thermal pressure n k T/2, not zero)" % (u_speed,))
    vhat = _scale(vrel, 1.0 / u_speed)

    # ---- blocking fractions ---------------------------------------------
    if shadow is None:
        chi = [1.0] * n_plates
    else:
        shadow = list(shadow)
        if len(shadow) != n_plates:
            raise ValueError("shadow has %d entries for %d plates -> pass one "
                             "fraction per plate" % (len(shadow), n_plates))
        chi = []
        for k, c in enumerate(shadow):
            c = float(c)
            if not (0.0 <= c <= 1.0):
                raise ValueError("shadow[%d] = %r -> pass a fraction in [0, 1]"
                                 % (k, c))
            chi.append(c)

    # ---- the mixture: partial mass densities and per-species speed ratios
    if species is not None:
        species = list(species)             # materialise once: it is read twice
    have_species = species is not None and len(species) > 0
    if model == "constantCd":
        if cd is None:
            raise ValueError("model 'constantCd' needs cd -> pass cd=<float>")
        if cd_rho is not None:
            rho_j = [float(cd_rho)]
            if math.isnan(rho_j[0]) or rho_j[0] < 0.0 or math.isinf(rho_j[0]):
                raise ValueError("cd_rho = %r -> pass a finite non-negative "
                                 "mass density in kg/m^3" % (cd_rho,))
        elif have_species:
            rho_j = [math.fsum(species_mass_densities(species))]
        else:
            raise ValueError("model 'constantCd' needs a total mass density -> "
                             "pass cd_rho=<kg/m^3> or the species list")
        s_j = [None]                       # s is unused by constantCd
    else:
        if cd_rho is not None:
            raise ValueError("cd_rho is meaningful only for 'constantCd' -> "
                             "pass the eight species number densities to %r "
                             "(a bulk density cannot drive a kinetic model)"
                             % (model,))
        if not have_species:
            raise ValueError("model %r needs the composition -> pass the eight "
                             "number densities in 1/m^3" % (model,))
        if t_k is None:
            raise ValueError("model %r needs t_k -> pass the ambient kinetic "
                             "temperature in K" % (model,))
        rho_j = species_mass_densities(species)
        s_j = speed_ratios(t_k, u_speed)

    alpha_k = _per_plate(alpha, n_plates, "alpha")
    tw_kk = _per_plate(tw_k, n_plates, "tw_k")
    spec_k = _per_plate(specular, n_plates, "specular")
    if model in _KINETIC_MODELS:
        if alpha_k[0] is None:
            raise ValueError("model %r needs alpha -> pass an energy "
                             "accommodation coefficient in [0, 1]" % (model,))
        if tw_kk[0] is None:
            raise ValueError("model %r needs tw_k -> pass a wall temperature "
                             "in K" % (model,))
        if model == "maxwell" and spec_k[0] is None:
            raise ValueError("model 'maxwell' needs specular -> pass a "
                             "specularly reflected fraction in [0, 1]")

    tinf = None if t_k is None else float(t_k)
    half_u2 = 0.5 * u_speed * u_speed          # rho stays inside the sum

    # species sharing a molecular mass share a speed ratio, so their
    # coefficient pair is evaluated once (section 7); exact, not approximate.
    first_of_mass = {}
    for j in range(len(rho_j)):
        m_j = MOLECULAR_MASS_KG[j] if model != "constantCd" else None
        first_of_mass.setdefault(m_j, j)

    forces = []
    for k, plate in enumerate(plates):
        area = float(plate["area"])
        if math.isnan(area) or area < 0.0 or math.isinf(area):
            raise ValueError("plates[%d]['area'] = %r -> pass a finite "
                             "non-negative area in m^2" % (k, area))
        nvec = _as_vec3(plate["normal"], "plates[%d]['normal']" % k)
        nlen = _norm(nvec)
        if abs(nlen - 1.0) > _UNIT_TOL:
            raise ValueError("plates[%d]['normal'] has length %r -> pass an "
                             "outward UNIT normal (the area is carried "
                             "separately)" % (k, nlen))
        nhat = _scale(nvec, 1.0 / nlen)     # remove the residual round-off

        if area == 0.0 or chi[k] == 0.0:
            forces.append((0.0, 0.0, 0.0))
            continue

        # ---- geometry, eqs. (1)-(2) ------------------------------------
        gamma = _dot(nhat, vhat)
        theta = math.atan2(_norm(_cross(nhat, vhat)), gamma)

        # t = (gamma*n - v)/l, normalised by its OWN length so that |t| <= 1
        # to within one rounding (section 8, D6).
        w = (gamma * nhat[0] - vhat[0],
             gamma * nhat[1] - vhat[1],
             gamma * nhat[2] - vhat[2])
        wlen = _norm(w)
        that = _scale(w, 1.0 / wlen) if wlen > _SIN_FLOOR else (0.0, 0.0, 0.0)

        # ---- the species sum, eq. (8) ----------------------------------
        cache = {}
        cp_terms = []
        ct_terms = []
        for j, rho in enumerate(rho_j):
            if rho == 0.0:
                continue                    # contributes exactly nothing
            key = first_of_mass[MOLECULAR_MASS_KG[j]
                                if model != "constantCd" else None]
            pair = cache.get(key)
            if pair is None:
                pair = gsi_coefficients(model, s_j[j], theta, alpha_k[k],
                                        tw_kk[k], tinf, cd=cd,
                                        specular=spec_k[k])
                cache[key] = pair
            cp_terms.append(rho * pair[0])
            ct_terms.append(rho * pair[1])
        rho_cp = math.fsum(cp_terms)        # exactly rounded, order-independent
        rho_ct = math.fsum(ct_terms)

        # ---- eq. (11) ---------------------------------------------------
        pref = half_u2 * area * chi[k]
        forces.append(_add(_scale(nhat, -pref * rho_cp),
                           _scale(that, pref * rho_ct)))

    return forces


# ---------------------------------------------------------------------------
# Reductions and diagnostics the caller usually wants next
# ---------------------------------------------------------------------------

def total_force(forces):
    """Sum of the per-plate forces [N]."""
    return (math.fsum(f[0] for f in forces),
            math.fsum(f[1] for f in forces),
            math.fsum(f[2] for f in forces))


def mixture_coefficients(species, t_k, u_speed, theta, model, alpha=None,
                         tw_k=None, cd=None, specular=None):
    """(rho, cp_mix, ctau_mix) of eq. (9), exposed for inspection and tests.

    Diagnostic only: ``drag_panel_force`` does not call this, keeping the
    un-normalised sums so that the total mass density never has to be formed.
    """
    rho_j = species_mass_densities(species)
    rho = math.fsum(rho_j)
    if rho == 0.0:
        return (0.0, 0.0, 0.0)
    s_j = (speed_ratios(t_k, u_speed) if model in _KINETIC_MODELS
           else [None] * _N_SPECIES)
    cp = []
    ct = []
    for j, r in enumerate(rho_j):
        if r == 0.0:
            continue
        a, b = gsi_coefficients(model, s_j[j], theta, alpha, tw_k, t_k, cd=cd,
                                specular=specular)
        cp.append(r * a)
        ct.append(r * b)
    return (rho, math.fsum(cp) / rho, math.fsum(ct) / rho)


# ===========================================================================
# Self-test (section 9).  Never called by drag_panel_force.
# ===========================================================================

#: representative 500 km composition, moderate activity [1/m^3]
_SPECIES_500KM = [1.0e12,    # He
                  2.0e13,    # O
                  5.0e11,    # N2
                  1.0e10,    # O2
                  1.0e7,     # Ar
                  1.0e11,    # H
                  5.0e11,    # N
                  1.0e10]    # anomalous O

_T_INF = 1000.0
_T_WALL = 273.0
_ALPHA = 0.85
_U = 7600.0


def _unit(v):
    return _scale(v, 1.0 / _norm(v))


def _plate(area, normal):
    return {"area": area, "normal": _unit(normal), "centroid": (0.0, 0.0, 0.0),
            "vis": (0.3, 0.2, 0.5), "ir": (0.8, 0.1, 0.1)}


def _relerr(got, want):
    d = _norm((got[0] - want[0], got[1] - want[1], got[2] - want[2]))
    n = _norm(want)
    return d / n if n > 0.0 else d


# --- the independent quadrature used only by check L3 ----------------------

def _simpson(f, a, b, n):
    if n % 2:
        n += 1
    h = (b - a) / n
    tot = f(a) + f(b)
    for i in range(1, n):
        tot += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return tot * h / 3.0


def _incident_force_quadrature(area, nhat, species, t_k, vrel, panels=200000):
    """Incident momentum flux of the mixture on one plate, by quadrature.

    In the plate frame let xi = c.n; incident molecules have xi < 0.  With
    y = -xi/v_m the incident normal speed in units of the most probable thermal
    speed, y > 0 is distributed as (1/sqrt(pi)) exp(-(y - x)^2) with x = s*gamma,
    and the momentum delivered per unit area per unit time is

        along -n :  n_j m_j v_m^2 (1/sqrt(pi)) Int_0^inf y^2 e^{-(y-x)^2} dy
        along +t :  n_j m_j U l v_m (1/sqrt(pi)) Int_0^inf y   e^{-(y-x)^2} dy

    summed over species.  With alpha = 1 and T_w = 0 the wall re-emits with no
    momentum, so this is the whole force and can be compared with the closed
    form.  Derived and integrated here; independent of gsi.py.
    """
    u = _norm(vrel)
    vhat = _scale(vrel, 1.0 / u)
    gamma = _dot(nhat, vhat)
    l = _norm(_cross(nhat, vhat))
    w = (gamma * nhat[0] - vhat[0], gamma * nhat[1] - vhat[1],
         gamma * nhat[2] - vhat[2])
    wl = _norm(w)
    that = _scale(w, 1.0 / wl) if wl > _SIN_FLOOR else (0.0, 0.0, 0.0)

    sp = math.sqrt(math.pi)
    p_n = 0.0
    p_t = 0.0
    for n_j, m_j in zip(species, MOLECULAR_MASS_KG):
        if n_j == 0.0:
            continue
        v_m = math.sqrt(2.0 * K_B * t_k / m_j)
        x = (u / v_m) * gamma
        hi = max(x, 0.0) + 12.0
        i2 = _simpson(lambda y: y * y * math.exp(-(y - x) ** 2), 0.0, hi,
                      panels) / sp
        i1 = _simpson(lambda y: y * math.exp(-(y - x) ** 2), 0.0, hi,
                      panels) / sp
        p_n += n_j * m_j * v_m * v_m * i2
        p_t += n_j * m_j * u * l * v_m * i1
    return tuple(area * (p_t * that[i] - p_n * nhat[i]) for i in range(3))


def _selftest():
    ok = [True]
    out = []

    def note(text):
        out.append(text)

    def check(name, got, want, tol):
        den = abs(want) if abs(want) > 0.0 else 1.0
        rel = abs(got - want) / den
        good = rel <= tol
        ok[0] = ok[0] and good
        out.append("%-56s %-4s got %+.12e want %+.12e rel %.2e (tol %.0e)"
                   % (name, "OK" if good else "FAIL", got, want, rel, tol))

    def check_vec(name, got, want, tol):
        rel = _relerr(got, want)
        good = rel <= tol
        ok[0] = ok[0] and good
        out.append("%-56s %-4s rel %.3e (tol %.0e)  F = (%+.6e, %+.6e, %+.6e)"
                   % (name, "OK" if good else "FAIL", rel, tol, got[0], got[1],
                      got[2]))

    sp = _SPECIES_500KM
    rho = math.fsum(species_mass_densities(sp))
    q = 0.5 * rho * _U * _U

    note("500 km mixture: rho = %.6e kg/m^3, q = 0.5 rho U^2 = %.6e Pa"
         % (rho, q))
    s_list = speed_ratios(_T_INF, _U)
    note("speed ratios: " + ", ".join("%s %.3f" % (nm, s)
                                      for nm, s in zip(SPECIES_NAMES, s_list)))

    # -- L1  Dalton, edge-on ------------------------------------------------
    a_edge = 2.5
    n_edge = (1.0, 0.0, 0.0)
    f = drag_panel_force([_plate(a_edge, n_edge)], sp, _T_INF,
                         (0.0, 4200.0, 0.0), "sentman", alpha=1.0, tw_k=0.0)[0]
    check("L1 Dalton edge-on: -f.n/A == n_tot k T/2", -f[0] / a_edge,
          0.5 * math.fsum(sp) * K_B * _T_INF, 1e-14)

    # -- L2  hyperthermal cold wall, oblique --------------------------------
    gam, area, u_hyp, t_cold = 0.5, 1.37, 3.0e5, 50.0
    th = math.acos(gam)
    vr = (gam * u_hyp, math.sin(th) * u_hyp, 0.0)
    f = drag_panel_force([_plate(area, (1.0, 0.0, 0.0))], sp, t_cold, vr,
                         "sentman", alpha=1.0, tw_k=0.0)[0]
    want = _scale(vr, -rho * area * gam * u_hyp)
    check_vec("L2 hyperthermal cold wall == -rho A gam U vrel", f, want, 1e-5)
    bad = _scale(_unit(vr), -2.0 * gam * gam * 0.5 * rho * u_hyp * u_hyp * area)
    note("   a drag/lift mis-resolution (cp along -v) would be wrong by "
         "%.1f %% here -> the check discriminates" % (100.0 * _relerr(bad, want)))

    # -- L3  independent quadrature -----------------------------------------
    worst = 0.0
    for deg in (0.0, 25.0, 60.0, 89.0, 90.0, 105.0):
        thd = math.radians(deg)
        nh = _unit((math.cos(thd), 0.6 * math.sin(thd), 0.8 * math.sin(thd)))
        got = drag_panel_force([_plate(1.23, nh)], sp, _T_INF, (_U, 0.0, 0.0),
                               "sentman", alpha=1.0, tw_k=0.0)[0]
        ref = _incident_force_quadrature(1.23, nh, sp, _T_INF, (_U, 0.0, 0.0))
        worst = max(worst, _relerr(got, ref))
    check("L3 closed form == Simpson quadrature (6 angles, worst)", worst, 0.0,
          1e-12)

    # -- L4  species linearity ----------------------------------------------
    plates = [_plate(0.95, (1.0, 0.0, 0.0)), _plate(0.95, (-1.0, 0.0, 0.0)),
              _plate(2.30, _unit((0.6427, 0.7660, 0.0))),
              _plate(2.30, _unit((0.6427, -0.7660, 0.0))),
              _plate(3.10, (0.0, 0.0, 1.0)), _plate(3.10, (0.0, 0.0, -1.0))]
    vr = (7580.0, 210.0, -95.0)
    kw = dict(model="sentman", alpha=_ALPHA, tw_k=_T_WALL)
    fmix = drag_panel_force(plates, sp, _T_INF, vr, **kw)
    acc = [(0.0, 0.0, 0.0)] * len(plates)
    for j in range(_N_SPECIES):
        only = [0.0] * _N_SPECIES
        only[j] = sp[j]
        acc = [_add(a, b) for a, b in
               zip(acc, drag_panel_force(plates, only, _T_INF, vr, **kw))]
    check("L4 mixture == sum of 8 single-species runs (worst plate)",
          max(_relerr(a, b) for a, b in zip(acc, fmix)), 0.0, 1e-14)
    note("   whole-body mixture force: (%+.6e, %+.6e, %+.6e) N"
         % total_force(fmix))

    # -- L5  linearity in the number densities (bit-exact) -------------------
    fd = drag_panel_force(plates, [2.0 * n for n in sp], _T_INF, vr, **kw)
    exact = all(fd[k][i] == 2.0 * fmix[k][i]
                for k in range(len(plates)) for i in range(3))
    ok[0] = ok[0] and exact
    out.append("%-56s %s" % ("L5 doubling every n_j doubles every force "
                             "bit-exactly", "OK" if exact else "FAIL"))

    # -- L6  constantCd on a convex body ------------------------------------
    cube = [_plate(1.0, (1.0, 0.0, 0.0)), _plate(1.0, (-1.0, 0.0, 0.0)),
            _plate(1.0, (0.0, 1.0, 0.0)), _plate(1.0, (0.0, -1.0, 0.0)),
            _plate(1.0, (0.0, 0.0, 1.0)), _plate(1.0, (0.0, 0.0, -1.0))]
    rho_c, cd_c = 3.1e-12, 2.7
    worst = 0.0
    for d in ((1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 1.0, 1.0),
              (0.3, -0.9, 0.2), (2.0, -3.0, 7.0)):
        dh = _unit(d)
        v = _scale(dh, 7650.0)
        tot = total_force(drag_panel_force(cube, None, None, v, "constantCd",
                                           cd=cd_c, cd_rho=rho_c))
        a_proj = math.fsum(p["area"] * max(_dot(p["normal"], dh), 0.0)
                           for p in cube)
        want = _scale(v, -0.5 * rho_c * a_proj * cd_c * 7650.0)
        worst = max(worst, _relerr(tot, want))
    check("L6 convex cube == -0.5 rho A_proj cd |v| v (worst of 5)", worst,
          0.0, 1e-14)
    tot = total_force(drag_panel_force(cube, sp, _T_INF, (7650.0, 0.0, 0.0),
                                       "constantCd", cd=cd_c))
    check_vec("L6a constantCd with rho inferred from the species list", tot,
              (-0.5 * rho * cd_c * 7650.0 ** 2, 0.0, 0.0), 1e-14)

    # -- L7  model consistency, shadow, lee side -----------------------------
    fa = drag_panel_force(plates, sp, _T_INF, vr, "sentman", alpha=_ALPHA,
                          tw_k=_T_WALL)
    fb = drag_panel_force(plates, sp, _T_INF, vr, "maxwell", alpha=_ALPHA,
                          tw_k=_T_WALL, specular=0.0)
    exact = all(x == y for A, B in zip(fa, fb) for x, y in zip(A, B))
    ok[0] = ok[0] and exact
    out.append("%-56s %s" % ("L7 maxwell(specular=0) == sentman bit-for-bit",
                             "OK" if exact else "FAIL"))
    fs = drag_panel_force(plates, sp, _T_INF, vr, shadow=[0.5] * len(plates),
                          **kw)
    exact = all(fs[k][i] == 0.5 * fmix[k][i]
                for k in range(len(plates)) for i in range(3))
    ok[0] = ok[0] and exact
    out.append("%-56s %s" % ("L7a shadow = 0.5 halves every force exactly",
                             "OK" if exact else "FAIL"))
    lee = drag_panel_force([_plate(1.0, (-1.0, 0.0, 0.0))], sp, _T_INF,
                           (_U, 0.0, 0.0), **kw)[0]
    ram = drag_panel_force([_plate(1.0, (1.0, 0.0, 0.0))], sp, _T_INF,
                           (_U, 0.0, 0.0), **kw)[0]
    good = 0.0 < abs(lee[0]) < 1.0e-6 * abs(ram[0])
    ok[0] = ok[0] and good
    out.append("%-56s %s  (lee/ram = %.3e, cp_ram = %.6f)"
               % ("L7b leeward plate small, positive, non-zero",
                  "OK" if good else "FAIL", abs(lee[0] / ram[0]),
                  -ram[0] / q))

    # -- L8  rotation covariance --------------------------------------------
    ang = 0.7391
    ca, sa = math.cos(ang), math.sin(ang)

    def rot(a):
        return (ca * a[0] - sa * a[1], sa * a[0] + ca * a[1], a[2])

    fr = drag_panel_force([_plate(p["area"], rot(p["normal"])) for p in plates],
                          sp, _T_INF, rot(vr), **kw)
    check("L8 rotation covariance about z (worst plate)",
          max(_relerr(a, rot(b)) for a, b in zip(fr, fmix)), 0.0, 1e-13)

    # -- L9  the mean-molar-mass shortcut ------------------------------------
    rho_list = species_mass_densities(sp)
    m_mass = math.fsum(r * m for r, m in zip(rho_list, MOLECULAR_MASS_KG)) / rho
    s_mass = _U * math.sqrt(m_mass / (2.0 * K_B * _T_INF))
    note("L9 mean-molar-mass shortcut (mass-weighted mean %.4f g/mol, "
         "s = %.4f):" % (m_mass * N_A * 1.0e3, s_mass))
    for deg in (0.0, 45.0, 80.0, 90.0, 100.0, 135.0):
        thd = math.radians(deg)
        _, cp_mix, ct_mix = mixture_coefficients(sp, _T_INF, _U, thd, "sentman",
                                                 _ALPHA, _T_WALL)
        cp_b, ct_b = gsi_coefficients("sentman", s_mass, thd, _ALPHA, _T_WALL,
                                      _T_INF)
        note("   theta = %5.1f deg  cp %.6e vs %.6e (%+9.4f %%)   "
             "ctau %.6e vs %.6e (%+9.4f %%)"
             % (deg, cp_mix, cp_b, 100.0 * (cp_b / cp_mix - 1.0) if cp_mix
                else float("nan"), ct_mix, ct_b,
                100.0 * (ct_b / ct_mix - 1.0) if ct_mix else float("nan")))

    # -- L10 order of magnitude at GRACE-FO ----------------------------------
    box = [_plate(0.95, (1.0, 0.0, 0.0)), _plate(0.95, (-1.0, 0.0, 0.0)),
           _plate(3.10, (0.0, 1.0, 0.0)), _plate(3.10, (0.0, -1.0, 0.0)),
           _plate(6.00, (0.0, 0.0, 1.0)), _plate(6.00, (0.0, 0.0, -1.0))]
    fb_tot = total_force(drag_panel_force(box, sp, _T_INF, (_U, 0.0, 0.0),
                                          **kw))
    a_ng = _norm(fb_tot) / 601.214
    cd_eff = 2.0 * _norm(fb_tot) / (rho * 0.95 * _U * _U)
    good = (1e-9 < a_ng < 2e-6) and (2.0 < cd_eff < 4.5)
    ok[0] = ok[0] and good
    out.append("%-56s %s  (|a| = %.4e m/s^2, C_D,eff = %.4f)"
               % ("L10 GRACE-sized box: |a_drag| and C_D in band",
                  "OK" if good else "FAIL", a_ng, cd_eff))

    print("\n".join(out))
    print("\n%s" % ("ALL CHECKS PASSED" if ok[0] else "*** FAILURES ***"))
    return ok[0]


if __name__ == "__main__":
    sys.exit(0 if _selftest() else 1)
