"""Free-molecular gas-surface interaction coefficients of a flat plate.

Independent numpy-free oracle for the NonGravForces kernel (WP-33).  Standard
library only (``math``).  The physics below was derived from kinetic theory for
this module; the MATLAB implementation it checks was not read.

    force per unit plate area = 0.5 * rho * U^2 * ( cp * (-n) + ctau * t )

RECONCILIATION.  This module is the referee's reconciliation of three
independent derivations (A, B, C) of the same coefficients.  All three agreed
on the physics -- see section 8 -- so the reconciliation is about the closed
form used and its conditioning, and the winning form (section 6) is one that
none of the three used.


===============================================================================
1.  GEOMETRY AND CONVENTIONS
===============================================================================

n       outward unit normal of the plate (away from the body, into space)
v       unit vector of the SATELLITE velocity relative to the gas; U = |v_rel|
g       the direction the gas streams in as seen from the plate, g = -v
theta   angle between n and v, in [0, pi]
gamma   cos(theta) = n.v ;  l = sin(theta) >= 0
t       unit tangential part of the stream, t = (g + gamma*n)/l

From g = -v it follows that

    g = l*t - gamma*n,                                                    (1)

so the gas bulk velocity U*g seen by the plate has normal component -U*gamma
(directed INTO the outward face when gamma > 0, i.e. windward) and tangential
component +U*l along t.

Species mass m, ambient kinetic temperature T_inf, mass density rho.  Most
probable thermal speed and speed ratio

    v_m = sqrt(2 k T_inf / m),      s = U / v_m.                          (2)

Let u_n = c.n be the normal component of a molecular velocity c in the plate
frame.  The gas is a Maxwellian centred on U*g, so u_n and the tangential
components are INDEPENDENT Gaussians with

    mean(u_n) = -U*gamma = -v_m * x,   x := s*gamma,   mean(u_t) = U*l,   (3)

x being the normal speed ratio: x > 0 windward, x < 0 leeward.  Molecules
strike the outward face when u_n < 0.


===============================================================================
2.  INCIDENT FLUXES (exact; no hyperthermal approximation)
===============================================================================

Put y = -u_n/v_m, so that y > 0 for a molecule moving toward the wall and y has
the Gaussian density pi^(-1/2) exp(-(y - x)^2).  The elementary moments over
y in (0, inf) are, with

    E = exp(-x^2),      Z = erfc(-x)  ( = 1 + erf(x) ),

    (1/sqrt(pi)) Int_0^inf     exp(-(y-x)^2) dy = Z/2
    (1/sqrt(pi)) Int_0^inf y   exp(-(y-x)^2) dy = E/(2 sqrt(pi)) + x*Z/2
    (1/sqrt(pi)) Int_0^inf y^2 exp(-(y-x)^2) dy = x*E/(2 sqrt(pi))
                                                  + (x^2 + 1/2)*Z/2       (4)

(the second from Int_{-inf}^{c} w e^{-w^2} dw = -e^{-c^2}/2, the third from
Int_{-inf}^{c} w^2 e^{-w^2} dw = -(c/2) e^{-c^2} + (sqrt(pi)/4)(1 + erf c),
both with c = -x after the shift w = y - x.)

Define the two dimensionless building blocks

    M(x) = E/sqrt(pi) + x*Z                      (mass-flux moment)       (5)
    W(x) = x*E/sqrt(pi) + (x^2 + 1/2)*Z          (normal-momentum moment) (6)

both STRICTLY POSITIVE for every finite x.  Then

  * incident mass flux per unit area      Phi_m = rho * v_m * M(x)/2
  * incident normal momentum flux         p_i   = rho * v_m^2 * W(x)/2,
    directed along -n (arriving molecules carry momentum m*u_n < 0),
  * incident tangential momentum flux     tau_i = Phi_m * (U*l),
    directed along +t.

The shear needs no new integral: u_t is independent of u_n in a Maxwellian, so
<u_t u_n>_{u_n<0} = <u_t> <u_n>_{u_n<0} = (U*l) * <u_n>_{u_n<0}.

Dividing by 0.5*rho*U^2 = 0.5*rho*s^2*v_m^2 and writing

    Phi = M(x)/s        (dimensionless mass flux 2*Phi_m/(rho*U)),        (7)

    cp_inc   = W(x) / s^2                                                 (8)
    ctau_inc = l * Phi = l * M(x) / s                                     (9)

Equivalently, with P = E/s and G = 1/(2 s^2),

    cp_inc   = gamma*P/sqrt(pi) + (gamma^2 + G) * Z                       (8')
    ctau_inc = l * ( P/sqrt(pi) + gamma * Z )                             (9')

The coefficient of Z in cp_inc is (gamma^2 + G).  The 1/(2 s^2) stands ALONE,
it is NOT multiplied by gamma^2, and ctau_inc contains no 1/(2 s^2) term at
all.  Section 7 explains why this is the decisive discriminator against the
drag/lift pair, and section 5(iii) gives the limit that proves it.


===============================================================================
3.  DIFFUSE RE-EMISSION (Sentman)
===============================================================================

All reflected molecules leave with a cosine-law half-Maxwellian at T_re and a
mass flux equal to the incident one (nothing accumulates).  For a
half-Maxwellian of most probable speed V_re = sqrt(2 k T_re/m) emitted into
u_n > 0, the momentum carried per unit mass flux is

    <u_n>_flux = Int_0^inf c^2 e^{-c^2/V_re^2} dc / Int_0^inf c e^{-c^2/V_re^2} dc
               = (V_re^3 sqrt(pi)/4) / (V_re^2/2) = (sqrt(pi)/2) * V_re,  (10)

so the reaction on the plate is a FURTHER inward pressure
p_re = (sqrt(pi)/2) * V_re * Phi_m.  The cosine law is azimuthally symmetric
about n, so the re-emitted stream carries ZERO net tangential momentum:

    *** the diffuse re-emission contributes nothing to ctau. ***
    *** ctau is therefore independent of alpha, T_w and T_inf. ***

In the 0.5 rho U^2 normalisation

    cp_re = (sqrt(pi)/2) * (V_re/U) * Phi,      ctau_re = 0.             (11)

Energy accommodation, as specified: the incident stream is given the
energy-equivalent temperature T_i = m U^2/(4 k), and
T_re = T_i - alpha*(T_i - T_w).  Since T_inf * s^2 = m U^2/(2 k), that is
T_i/T_inf = s^2/2 and

    T_re/T_inf = (1 - alpha) * s^2/2 + alpha * T_w/T_inf
    V_re/U     = sqrt( (1 - alpha)/2 + alpha * T_w/(T_inf * s^2) ).      (12)

The molecular mass cancels: only (s, alpha, T_w, T_inf) enter the interface.
The second form of (12) is a sum of two NON-NEGATIVE terms, so it is exact at
alpha = 0 and alpha = 1 with no cancellation; the algebraically equivalent
sqrt( (1/2)(1 + alpha(4 k T_w/(m U^2) - 1)) ) is not, and is not used.

    SENTMAN:   cp = cp_inc + cp_re,     ctau = ctau_inc.                 (13)


===============================================================================
4.  MAXWELL AND constantCd
===============================================================================

MAXWELL.  A fraction sigma reflects specularly: u_n is reversed and u_t kept,
so such a molecule delivers 2*m*|u_n| along -n and nothing along t, and it does
not re-emit thermally.  The rest behaves as Sentman.  Exactly (no hyperthermal
approximation):

    cp   = 2*sigma*cp_inc + (1-sigma)*(cp_inc + cp_re)
         = (1 + sigma)*cp_inc + (1 - sigma)*cp_re                        (14)
    ctau = (1 - sigma)*ctau_inc                                          (15)

sigma = 0 reproduces Sentman exactly; sigma = 1 gives cp = 2*cp_inc, ctau = 0.

constantCd.  Every windward plate carries a force along the gas stream g of
magnitude 0.5*rho*cd*A*gamma*U^2, so that the sum over a body reproduces
F = -0.5*rho*cd*A_projected*|v_rel|*v_rel with A_proj = sum A_k gamma_k.
Resolving cd*gamma*g on (-n, t) with (1):

    cp = cd*gamma^2,   ctau = cd*gamma*l   (gamma > 0);   0, 0 otherwise. (16)

s, alpha, T_w and T_inf are unused.


===============================================================================
5.  LIMITS THAT PIN THE RESULT (all checked in _selftest)
===============================================================================

(i)   Head-on hyperthermal, alpha = 1, cold wall: cp -> 2.  The plate absorbs
      the whole incident momentum flux rho U^2, and rho U^2/(0.5 rho U^2) = 2.
(ii)  Same, purely specular: cp -> 4, the classical mirror-plate value.
(iii) EDGE-ON (gamma = 0, so x = 0, E = 1, Z = 1): W(0) = 1/2 gives
      cp_inc = 1/(2 s^2) exactly, i.e. a dimensional pressure
      0.5*rho*U^2/(2 s^2) = rho*v_m^2/4 = rho*k*T_inf/(2 m) = p/2, the
      one-sided thermal momentum flux of a gas at rest.  M(0) = 1/sqrt(pi)
      gives ctau_inc = l/(s sqrt(pi)) exactly.  ANY form in which the
      1/(2 s^2) term is multiplied by gamma^2 fails this check, and any form
      in which ctau carries a lone 1/(2 s^2) term fails it too.  This is the
      limit that discriminates (cp, ctau) from (C_D, C_L) -- see section 7.
(iv)  Gas at rest in equilibrium with the wall (s -> 0, alpha = 1,
      T_w = T_inf): cp_inc -> p/2 and cp_re -> p/2, so the plate feels exactly
      the static pressure p = rho k T_inf/m on the exposed face.
(v)   Hyperthermal, alpha = 1, cold plate, arbitrary theta: cp -> 2 gamma^2 and
      ctau -> 2 gamma l, i.e. exactly constantCd with cd = 2, as pure momentum
      absorption on the projected area must give.
(vi)  Leeward: everything is O(exp(-x^2)) and both coefficients stay positive.


===============================================================================
6.  CONDITIONING -- THE FORM ACTUALLY USED  (the reconciliation's substance)
===============================================================================

For a leeward plate at high speed ratio (x = -u, u >> 1) the direct forms (5)
and (6) are catastrophically cancelling differences:

    M(-u) = E/sqrt(pi) - u*erfc(u)             ~ E/(2 sqrt(pi) u^2)
    W(-u) = (u^2 + 1/2)*erfc(u) - u*E/sqrt(pi) ~ E/(2 sqrt(pi) u^3)

while the individual terms are O(E/u) and O(u E): the relative cancellation is
~2 u^2 for M and ~2 u^4 for W.  Measured against an 80-digit reference the
direct form loses up to 3.5e-11 relative at u = 26.  An asymptotic series in
1/u fixes that down to its optimal-truncation floor of ~1e-13 near u = 6.

This module uses a form with NO cancellation at all and no truncation floor.
From the continued fraction of the complementary error function (valid for
u > 0),

    sqrt(pi) * exp(u^2) * erfc(u) = 1/(u + T1),
    T_k = (k/2)/(u + T_{k+1}),   T1 = (1/2)/(u + T2),  T2 = 1/(u + T3), ... (17)

all T_k > 0.  Writing R(u) = sqrt(pi) u exp(u^2) erfc(u) = u/(u + T1),

    M(-u) * sqrt(pi)/E = 1 - R(u)          = T1/(u + T1)                 (18)

    W(-u) * sqrt(pi)/E = R(u)/(2u) - u*(1 - R(u))
                       = [1/2 - u*T1]/(u + T1)
                       = T2 / ( 2*(u + T2)*(u + T1) )                    (19)

the last step using u*T1 = u/(2*(u + T2)), whence
1/2 - u*T1 = [ (u + T2) - u ] / (2*(u + T2)) = T2/(2*(u + T2)).

Every quantity in (18) and (19) is positive: the two leading-order
cancellations of M and W have been performed ALGEBRAICALLY, once and for all,
by descending two levels of the continued fraction.  The overall exp(-u^2) is
carried as an explicit factor and never cancels.  Verified against the 80-digit
reference: <= 4e-16 relative for every u in [2, 27], versus 1e-13 (asymptotic
series) and 3.5e-11 (direct).

The continued fraction converges geometrically for u >~ 2 and slowly below it,
so:
  * x > -1   : the direct forms (5) and (6) are used; the cancellation there is
               at most a factor 2 u^4 = 2, i.e. nothing, and the measured
               error is <= 3e-16.
  * x <= -1  : the continued fraction (18)-(19), with the number of levels
               taken as 24 + int(200/u), capped at 320 (224 levels at the
               switch).  The switch at u = 1 was placed by measurement: over
               u = 1(1/8)5 the continued fraction is <= 4e-16 throughout while
               the direct form ranges from 6e-16 to 2e-13, and below u = 1 the
               continued fraction's convergence degrades (5e-12 at u = 0.5)
               while the direct form has nothing left to cancel.
  * u^2 > 700: exp(-u^2) is no longer a normal double, so the product
               exp(-u^2)*f would be built from a crushed subnormal mantissa
               (5.5e-3 relative at u = 27).  The exponents are combined
               instead, exp(-u^2 + log f), which holds ~1e-13 relative and
               underflows to exactly 0.0 only when the true value is below the
               smallest subnormal (u > 27.3, coefficients below 1e-320).

One more error source dominates once |x| is large and has nothing to do with
the cancellation: rounding u*u before exponentiating.  The rounded u^2 is wrong
by up to ulp(u^2)/2, an ABSOLUTE error of ~5e-14 at u = 22, and exp turns that
into the same RELATIVE error on the answer.  _exp_neg_square removes it with a
Veltkamp/Dekker split of the product (exact head + tail) and one first-order
correction, taking the leeward moments from ~1.5e-14 to ~1e-16.  What remains
is the amplification of the caller's own rounding of x = s*cos(theta), which is
2|x| times one ulp (4.4e-14 at x = -22) and is a property of the argument, not
of this module.

Two further conditioning choices:
  * Z is always formed as erfc(-x), NEVER as 1 + erf(x): at x = -8 the latter
    is 0.0 in double while erfc(8) = 1.1e-29.
  * (V_re/U)^2 is formed as in (12) as a sum of non-negative terms.

Everything is closed form.  The quadrature in _selftest is a verification of
the algebra only; gsi_coefficients never calls it.


===============================================================================
7.  CONTRADICTION WITH THE REFERENCE DOCUMENT (theory.txt section 2.3)
===============================================================================

theory.txt equations (4) and (5), labelled "Cp" and "C_tau" and fed into its
equation (1) as the normal and tangential coefficients, are in fact the DRAG
and LIFT coefficients.  Its (4)-(5) read, with P = E/s, Q = 1 + 1/(2 s^2),
G = 1/(2 s^2), Z = 1 + erf(gamma s),

    "Cp"    = P/sqrt(pi) + gamma*Q*Z + (gamma/2)*(V_re/U)*(sqrt(pi)*gamma*Z + P)
    "C_tau" =              l*G*Z     + (l/2)    *(V_re/U)*(sqrt(pi)*gamma*Z + P)

Rotating MY result onto the stream direction g and its perpendicular
e_L = l*n - gamma*t gives exactly those two expressions:

    C_D = (cp*(-n) + ctau*t) . g = cp*gamma + ctau*l
        = gamma^2*P/sqrt(pi) + gamma*(gamma^2+G)*Z + l^2*P/sqrt(pi) + l^2*gamma*Z
        = P/sqrt(pi) + gamma*(1+G)*Z = P/sqrt(pi) + gamma*Q*Z    [+ gamma*cp_re]

    C_L = (cp*(-n) + ctau*t) . e_L = cp*l - ctau*gamma
        = l*gamma*P/sqrt(pi) + l*(gamma^2+G)*Z
          - gamma*l*P/sqrt(pi) - gamma^2*l*Z
        = l*G*Z                                                  [+ l*cp_re]

so theory.txt's "Cp" is C_D and its "C_tau" is C_L.  The inverse rotation, if
one starts from the published pair, is cp = gamma*C_D + l*C_L,
ctau = l*C_D - gamma*C_L.  The tells are the Q = 1 + 1/(2 s^2) multiplying Z
(the true cp carries gamma^2 + 1/(2 s^2)) and the lone G*Z in "C_tau" (the true
shear has no 1/(2 s^2) term whatever).  The two pairs COINCIDE at normal
incidence (l = 0) and differ at every inclined panel; a leeward panel acquires
a spurious O(1/s^2) shear.  Following the instruction to trust the derivation
over the document: this module implements the true (cp, ctau), and the
edge-on limit 5(iii) is the test that decides it.  All three independent
derivations reached the same conclusion.


===============================================================================
8.  WHAT THE THREE DERIVATIONS DISAGREED ABOUT
===============================================================================

Nothing in the physics.  Over 36 199 grid points (s in [0.01, 100], theta over
[0, pi] including 0, pi/2 and pi, alpha in {0, 0.5, 0.85, 1}, five T_w/T_inf
ratios, sigma in {0, 0.2, 0.5, 0.9, 1}, cd in {2, 2.2, 4}) the three agreed to
3.6e-12 absolute; the largest RELATIVE spread, 2.8e-8, was entirely
floating-point conditioning in the deep leeward tail on values of order
1e-280.  None of the three had computed the drag/lift pair by mistake; all
three independently flagged the theory.txt mislabelling of section 7.
Settled by derivation and by measurement against an 80-digit reference, not by
majority:
  * the leeward form.  A and C switch to an asymptotic series only at x = -6
    and use the direct erfc difference above it, so they carry its cancellation
    through the whole band -6 < x < -1; B never switches at all.  Worst
    relative error of the normal-momentum moment W over x = -0.25 ... -27,
    evaluated at the same double x so that the caller's own rounding is
    excluded:

        A 6.4e-13     B 2.4e-10     C 6.4e-13     this module 1.0e-15

    Neither route is used here.  The continued fraction of section 6 is
    cancellation-free at every x and has no truncation floor, and the
    Veltkamp-corrected exp(-u^2) removes the error that then dominates.  The
    gap is not confined to the unphysical tail: at s = 15, theta = 1.9635
    (x = -5.7, an ordinary leeward panel at GRACE-FO speed ratios) all three
    variants give 9.0e-12 and this module 2.9e-15.
  * C returns exactly 0.0 for u^2 > 700, i.e. from u = 26.46 on, where the
    coefficients are still representable (3.8e-311 at u = 26.5).  A stops at
    u = 27, B at u^2 = 745.  This module carries the value to the last
    representable subnormal -- measured within 0.5 ulp at u = 26, 26.5 and 27 --
    and returns 0.0 only once the exact answer is below the smallest subnormal.
  * T_w = 0 K: A and C admit it, B rejects it.  Admitted: T_re = (1-alpha) T_i
    + alpha T_w is well defined and non-negative at T_w = 0, and the
    hyperthermal cold-plate limit 5(i)/5(v) is a standard test point.
  * an infinite s: A admitted it (returning nan), B and C reject it.  Rejected.
  * theta round-off: C admits a 1e-12 slack outside [0, pi], A and B do not.
    Admitted and clamped, since a caller forming theta = acos(n.v) can land a
    few ulp outside the interval.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

import math

__all__ = ["gsi_coefficients"]

_SQRT_PI = math.sqrt(math.pi)
_MODELS = ("sentman", "maxwell", "constantCd")

_CF_SWITCH = 1.0        # use the continued fraction for x <= -_CF_SWITCH
_NORMAL_LIMIT = 700.0   # exp(-u^2) is still a normal double below this
_THETA_SLACK = 1.0e-12  # round-off tolerance on the [0, pi] interval


def _exp_neg_square(u):
    """exp(-u^2) to full double precision, for 0 < u < 2^400.

    Forming u*u first loses the low half of the product: at u = 22 the rounded
    u^2 is wrong by ~5e-14 ABSOLUTE, and exp turns that straight into a 5e-14
    RELATIVE error on the result -- the dominant error of every leeward
    coefficient once |x| is large.  Veltkamp/Dekker splitting gives the product
    exactly as an unevaluated sum p + err with |err| <= ulp(p)/2, and then

        exp(-(p + err)) = exp(-p) * exp(-err) = exp(-p) * (1 - err + ...)

    with |err| < 1e-13, so a single first-order correction is enough.  Pure
    double arithmetic; no dependence on math.fma or any newer runtime.
    """
    p = u * u
    c = 134217729.0                      # 2^27 + 1
    t = c * u
    u_hi = t - (t - u)                   # exact 26-bit head
    u_lo = u - u_hi                      # exact tail
    err = ((u_hi * u_hi - p) + 2.0 * u_hi * u_lo) + u_lo * u_lo
    return math.exp(-p) * (1.0 - err)


def _cf_tails(u):
    """T1 and T2 of the erfc continued fraction (17), by backward recurrence.

    T_k = (k/2)/(u + T_{k+1}), truncated at a level chosen so that the
    truncation error is below the double round-off for every u >= _CF_SWITCH.
    Both returned values are strictly positive.
    """
    levels = 24 + int(200.0 / u)
    if levels > 320:
        levels = 320
    t = 0.0
    for k in range(levels, 1, -1):
        t = (0.5 * k) / (u + t)
    t2 = t
    t1 = 0.5 / (u + t2)
    return t1, t2


def _flux_moments(x):
    """Return (W, M) of equations (6) and (5) at the normal speed ratio x.

        M(x) = exp(-x^2)/sqrt(pi) + x*erfc(-x)              (mass flux)
        W(x) = x*exp(-x^2)/sqrt(pi) + (x^2 + 1/2)*erfc(-x)  (normal momentum)

    Both are strictly positive for every finite x.  Windward and mildly
    leeward values use the direct form; deep leeward values use the
    cancellation-free continued fraction (18)-(19).
    """
    if x > -_CF_SWITCH:
        e = math.exp(-x * x)
        z = math.erfc(-x)                     # = 1 + erf(x), stable branch
        m = e / _SQRT_PI + x * z
        w = x * e / _SQRT_PI + (x * x + 0.5) * z
        return w, m

    u = -x
    uu = u * u
    t1, t2 = _cf_tails(u)
    f_m = t1 / ((u + t1) * _SQRT_PI)                            # eq. (18)
    f_w = t2 / (2.0 * (u + t2) * (u + t1) * _SQRT_PI)           # eq. (19)
    if uu <= _NORMAL_LIMIT:
        e = _exp_neg_square(u)
        return e * f_w, e * f_m
    # Below the smallest NORMAL double the product exp(-u^2)*f would be formed
    # from a subnormal factor and lose most of its mantissa (5e-3 relative at
    # u = 27).  Combining the exponents keeps ~1e-13 relative instead, and
    # underflows to exactly 0.0 only when the answer itself is unrepresentable.
    return math.exp(-uu + math.log(f_w)), math.exp(-uu + math.log(f_m))


def gsi_coefficients(model, s, theta, alpha, tw_k, tinf_k, cd=None,
                     specular=None):
    """Free-molecular pressure and shear coefficients of a flat plate.

        force per unit plate area = 0.5*rho*U^2*( cp*(-n) + ctau*t )

    Parameters
    ----------
    model     "sentman", "maxwell" or "constantCd".
    s         molecular speed ratio U/sqrt(2 k T_inf/m) of the species, > 0
              and finite.  Unused by "constantCd" (may be None).
    theta     angle between the outward normal n and the satellite velocity
              direction v, in [0, pi] radians; theta = 0 is head-on,
              theta > pi/2 is leeward.
    alpha     energy accommodation coefficient in [0, 1].  Unused by
              "constantCd".
    tw_k      wall temperature [K], >= 0.  Unused by "constantCd".
    tinf_k    ambient gas kinetic temperature [K], > 0.  Unused by
              "constantCd".
    cd        drag coefficient; required by "constantCd" only.
    specular  specularly reflected fraction in [0, 1]; required by "maxwell"
              only.

    Returns
    -------
    (cp, ctau), both >= 0 for every admissible argument.

    Notes
    -----
    ctau is independent of alpha, tw_k and tinf_k: a cosine-law re-emission
    carries no net tangential momentum (section 3 of the module docstring).
    """
    if model not in _MODELS:
        raise ValueError("unknown model %r -> use one of %s" % (model, _MODELS))

    if theta is None:
        raise ValueError("theta is required -> pass an angle in [0, pi] radians")
    theta = float(theta)
    if not (-_THETA_SLACK <= theta <= math.pi + _THETA_SLACK):
        raise ValueError("theta = %r out of range -> pass an angle in [0, pi] "
                         "radians" % (theta,))
    if theta < 0.0:
        theta = 0.0
    elif theta > math.pi:
        theta = math.pi

    gamma = math.cos(theta)
    l = math.sin(theta)
    if l < 0.0:                       # only reachable through round-off
        l = 0.0

    # ---- constantCd: pure geometry, eq. (16) ----------------------------
    if model == "constantCd":
        if cd is None:
            raise ValueError("model 'constantCd' needs cd -> pass cd=<float>")
        cd = float(cd)
        if gamma <= 0.0:              # leeward plate contributes nothing
            return (0.0, 0.0)
        return (cd * gamma * gamma, cd * gamma * l)

    # ---- validation shared by the kinetic models ------------------------
    if s is None:
        raise ValueError("s is required -> pass a finite speed ratio > 0")
    s = float(s)
    if not (s > 0.0) or math.isinf(s):
        raise ValueError("s = %r invalid -> pass a finite speed ratio > 0"
                         % (s,))
    if alpha is None:
        raise ValueError("alpha is required -> pass a value in [0, 1]")
    alpha = float(alpha)
    if not (0.0 <= alpha <= 1.0):
        raise ValueError("alpha = %r out of range -> pass a value in [0, 1]"
                         % (alpha,))
    if tw_k is None or tinf_k is None:
        raise ValueError("tw_k and tinf_k are required -> pass temperatures "
                         "in K")
    tw_k = float(tw_k)
    tinf_k = float(tinf_k)
    if not (tw_k >= 0.0):
        raise ValueError("tw_k = %r invalid -> pass a wall temperature >= 0 K"
                         % (tw_k,))
    if not (tinf_k > 0.0):
        raise ValueError("tinf_k = %r invalid -> pass an ambient temperature "
                         "> 0 K" % (tinf_k,))

    sigma = 0.0
    if model == "maxwell":
        if specular is None:
            raise ValueError("model 'maxwell' needs specular -> pass "
                             "specular=<fraction in [0, 1]>")
        sigma = float(specular)
        if not (0.0 <= sigma <= 1.0):
            raise ValueError("specular = %r out of range -> pass a value in "
                             "[0, 1]" % (sigma,))

    # ---- kinetic coefficients, eqs. (8), (9), (11) ----------------------
    x = s * gamma                                     # normal speed ratio
    w_moment, m_moment = _flux_moments(x)

    cp_inc = w_moment / (s * s)                       # eq. (8)
    phi = m_moment / s                                # eq. (7)
    ctau_inc = l * phi                                # eq. (9)

    # re-emission speed ratio, eq. (12): a sum of non-negative terms
    v_re_over_u = math.sqrt(0.5 * (1.0 - alpha)
                            + alpha * tw_k / (tinf_k * s * s))
    cp_re = 0.5 * _SQRT_PI * v_re_over_u * phi        # eq. (11)

    if model == "sentman":                            # eq. (13)
        return (cp_inc + cp_re, ctau_inc)

    cp = (1.0 + sigma) * cp_inc + (1.0 - sigma) * cp_re      # eq. (14)
    ctau = (1.0 - sigma) * ctau_inc                          # eq. (15)
    return (cp, ctau)


# ---------------------------------------------------------------------------
# Verification of the algebra.  Never called by gsi_coefficients.
# ---------------------------------------------------------------------------

def _moments_quadrature(x, panels=40000):
    """(W, M) by cancellation-free composite Simpson quadrature of (4).

    With y the inward normal speed in units of v_m,

        M(x) = (2/sqrt(pi)) Int_0^inf y   exp(-(y-x)^2) dy
        W(x) = (2/sqrt(pi)) Int_0^inf y^2 exp(-(y-x)^2) dy
             = (2 exp(-x^2)/sqrt(pi)) Int_0^inf y^k exp(2xy - y^2) dy

    The second form carries the tiny exponential as an explicit factor, so the
    integrand stays O(1) for x << 0 and nothing cancels.  For x << 0 the
    integrand decays on the scale 1/(2|x|), so y is stretched by
    b = 1 + 2 max(0, -x) to keep the sampled peak resolved.
    """
    if x >= 0.0:
        b, z_hi = 1.0, x + 9.0
    else:
        b, z_hi = 1.0 + 2.0 * (-x), 60.0
    h = z_hi / panels
    s1 = s2 = 0.0
    for i in range(1, panels):
        y = i * h / b
        f = math.exp(2.0 * x * y - y * y)
        wt = 4.0 if (i % 2) else 2.0
        s1 += wt * y * f
        s2 += wt * y * y * f
    y = z_hi / b                        # y = 0 contributes nothing to either
    f = math.exp(2.0 * x * y - y * y)
    s1 += y * f
    s2 += y * y * f
    scale = 2.0 * math.exp(-x * x) * (h / b) / (3.0 * _SQRT_PI)
    return scale * s2, scale * s1


def _selftest():
    ok = [True]

    def check(name, got, want, tol):
        """Relative check; `want` must be non-zero (or `got` exactly zero)."""
        if want == 0.0:
            err = 0.0 if got == 0.0 else float("inf")
        else:
            err = abs(got - want) / abs(want)
        good = err <= tol
        if not good:
            ok[0] = False
        print("  %s %-46s got %.15g  want %.15g  rel %.2e"
              % ("ok  " if good else "FAIL", name, got, want, err))

    def check_le(name, got, bound):
        """Absolute check: `got` must not exceed a bound derived above."""
        good = got <= bound
        if not good:
            ok[0] = False
        print("  %s %-46s got %.3e  bound %.3e"
              % ("ok  " if good else "FAIL", name, got, bound))

    print("A. closed form vs quadrature of eq. (4)")
    worst = 0.0
    for i in range(-120, 121, 10):
        x = i / 10.0
        wc, mc = _flux_moments(x)
        wq, mq = _moments_quadrature(x)
        worst = max(worst, abs(wc - wq) / wq, abs(mc - mq) / mq)
    # the bound is the Simpson rule's own error, not the closed form's
    check_le("worst |closed/quadrature - 1| over x in [-12, 12]", worst, 5e-12)

    print("B. positivity and monotonicity of the moments")
    bad = 0
    prev_w = prev_m = -1.0
    for i in range(-2600, 2601):
        x = i / 100.0
        w, m = _flux_moments(x)
        if w < 0.0 or m < 0.0:
            bad += 1
        if w < prev_w - 1e-300 or m < prev_m - 1e-300:
            bad += 1
        prev_w, prev_m = w, m
    check_le("negative or non-monotone moment count", float(bad), 0.0)

    print("C. the two branches agree AT THE SAME x")
    # (comparing _flux_moments either side of the switch would measure dW/dx,
    #  which is O(2u + 3/u) per unit x and says nothing about the branches)
    for x in (-1.0, -1.5, -2.0, -3.0, -5.0):
        e = math.exp(-x * x)
        z = math.erfc(-x)
        w_direct = x * e / _SQRT_PI + (x * x + 0.5) * z
        m_direct = e / _SQRT_PI + x * z
        w_cf, m_cf = _flux_moments(x)          # takes the continued fraction
        check("W direct vs continued fraction at x=%.1f" % x,
              w_cf, w_direct, 3e-13)
        check("M direct vs continued fraction at x=%.1f" % x,
              m_cf, m_direct, 3e-14)

    print("D. limits of section 5")
    s = 5000.0
    cp, ct = gsi_coefficients("sentman", s, 0.0, 1.0, 0.0, 1000.0)
    check("(i)  head-on hyperthermal cold: cp = 2 + 1/s^2",
          cp, 2.0 + 1.0 / (s * s), 1e-15)
    check("(i)  head-on: ctau = 0", ct, 0.0, 0.0)
    cp, ct = gsi_coefficients("maxwell", s, 0.0, 1.0, 0.0, 1000.0,
                              specular=1.0)
    check("(ii) head-on specular: cp = 2(2 + 1/s^2)",
          cp, 2.0 * (2.0 + 1.0 / (s * s)), 1e-15)
    for s_edge in (0.5, 2.0, 8.0):
        cp_inc = _flux_moments(0.0)[0] / (s_edge * s_edge)
        check("(iii) edge-on cp_inc = 1/(2 s^2), s=%g" % s_edge,
              cp_inc, 1.0 / (2.0 * s_edge * s_edge), 0.0)
        _, ct = gsi_coefficients("sentman", s_edge, math.pi / 2, 0.85,
                                 300.0, 1000.0)
        check("(iii) edge-on ctau = 1/(s sqrt(pi)), s=%g" % s_edge,
              ct, 1.0 / (s_edge * _SQRT_PI), 1e-15)
    # (iv) s -> 0 with alpha = 1, T_w = T_inf.  Expanding (8), (11) to O(s):
    #   W(s) -> 1/2 + 2s/sqrt(pi),  M(s) -> 1/sqrt(pi) + s,  V_re/U = 1/s,
    #   cp*s^2 -> 1 + (2/sqrt(pi) + sqrt(pi)/2) s + O(s^2).
    slope = 2.0 / _SQRT_PI + _SQRT_PI / 2.0
    for s_small in (1e-2, 1e-3, 1e-4):
        cp, _ = gsi_coefficients("sentman", s_small, 0.0, 1.0, 1000.0, 1000.0)
        check_le("(iv) gas at rest: |cp s^2 - 1 - %.4f s| = O(s^2), s=%g"
                 % (slope, s_small),
                 abs(cp * s_small * s_small - 1.0 - slope * s_small),
                 5.0 * s_small * s_small)
    # (v) at s = 5000 the hyperthermal limit is reached EXACTLY in double
    # (E = 0, Z = 2), so cp = 2 gamma^2 + 1/s^2 and ctau = 2 gamma l exactly.
    s_hyper = 5000.0
    for th in (0.0, 0.4, 1.0, 1.5):
        cp, ct = gsi_coefficients("sentman", s_hyper, th, 1.0, 0.0, 1000.0)
        cp2, ct2 = gsi_coefficients("constantCd", None, th, None, None, None,
                                    2.0)
        check("(v) hyperthermal cp = cd2 + 1/s^2, th=%.1f" % th,
              cp, cp2 + 1.0 / (s_hyper * s_hyper), 1e-14)
        check("(v) hyperthermal ctau = cd2 ctau, th=%.1f" % th, ct, ct2, 1e-14)
    for th in (1.8, 2.2, 2.6, 3.0):
        cp, ct = gsi_coefficients("sentman", 8.0, th, 0.85, 300.0, 1000.0)
        check("(vi) leeward cp > 0, th=%.1f" % th, 1.0 if cp > 0 else 0.0,
              1.0, 0.0)
        check("(vi) leeward ctau > 0, th=%.1f" % th, 1.0 if ct > 0 else 0.0,
              1.0, 0.0)
    # the deep tail must still be positive, not flushed to zero prematurely
    for s_deep, th in ((8.0, 2.9), (20.0, 2.6), (26.0, 3.0)):
        cp, ct = gsi_coefficients("sentman", s_deep, th, 0.85, 300.0, 1000.0)
        check("(vi) deep tail cp > 0, s=%g th=%.1f" % (s_deep, th),
              1.0 if cp > 0 else 0.0, 1.0, 0.0)

    print("E. ctau does not depend on alpha, T_w or T_inf")
    ref = gsi_coefficients("sentman", 8.0, 1.0, 0.85, 300.0, 1000.0)[1]
    for al, tw, ti in ((0.0, 100.0, 400.0), (1.0, 900.0, 1200.0),
                       (0.5, 0.0, 700.0)):
        check("ctau invariant (alpha=%g, Tw=%g, Tinf=%g)" % (al, tw, ti),
              gsi_coefficients("sentman", 8.0, 1.0, al, tw, ti)[1], ref, 0.0)

    print("F. drag/lift rotation reproduces the published Sentman pair")
    for s, th, al, tw, ti in ((8.0, 0.7, 0.85, 300.0, 1000.0),
                              (3.0, 2.2, 0.40, 300.0, 900.0),
                              (0.5, 1.1, 1.00, 300.0, 950.0)):
        cp, ct = gsi_coefficients("sentman", s, th, al, tw, ti)
        g, l = math.cos(th), math.sin(th)
        p = math.exp(-(s * g) ** 2) / s
        z = math.erfc(-s * g)
        q = 1.0 + 0.5 / (s * s)
        gg = 0.5 / (s * s)
        vre = math.sqrt(0.5 * (1.0 - al) + al * tw / (ti * s * s))
        re = 0.5 * vre * (_SQRT_PI * g * z + p)
        check("C_D  s=%g th=%.1f" % (s, th), cp * g + ct * l,
              p / _SQRT_PI + g * q * z + g * re, 5e-14)
        check("C_L  s=%g th=%.1f" % (s, th), cp * l - ct * g,
              l * gg * z + l * re, 5e-11)

    print("G. maxwell consistency")
    a = gsi_coefficients("sentman", 8.0, 0.9, 0.85, 300.0, 1000.0)
    b = gsi_coefficients("maxwell", 8.0, 0.9, 0.85, 300.0, 1000.0,
                         specular=0.0)
    check("specular = 0 equals sentman (cp)", b[0], a[0], 0.0)
    check("specular = 0 equals sentman (ctau)", b[1], a[1], 0.0)
    c = gsi_coefficients("maxwell", 8.0, 0.9, 0.85, 300.0, 1000.0,
                         specular=1.0)
    check("specular = 1 gives cp = 2 cp_inc", c[0],
          2.0 * _flux_moments(8.0 * math.cos(0.9))[0] / 64.0, 1e-15)
    check("specular = 1 gives ctau = 0", c[1], 0.0, 0.0)

    print("SELF-TEST", "PASS" if ok[0] else "FAIL")
    return ok[0]


if __name__ == "__main__":
    raise SystemExit(0 if _selftest() else 1)
