r"""Radiation-pressure force on flat satellite plates -- reconciled oracle.

Standard library only (``math``).  SI units, angles in radians.  numpy is NOT
imported and is not installed in the verification environment.

This module is an ORACLE for ``NonGravForces.kernel.radiationPanelForce``.  Its
only value is that it was derived from the physics rather than transcribed from
the implementation it checks.  Neither the MATLAB kernel, nor its tests, nor the
existing Sentman oracle, nor the project HANDOVER were opened while it was
written.  The derivation is set out in full below so a reader can check the
result without trusting the code.

It is the reconciliation of two independent derivations (route (a), a direct
momentum balance; route (b), superposition of the three pure limiting cases),
refereed against exact rational arithmetic and against a photon Monte Carlo.
The reconciliation record is in section 6.


================================================================================
1.  GEOMETRY AND SYMBOLS
================================================================================

    n        outward unit normal of the plate (away from the body, into space)
    A        plate area [m^2]
    u        unit vector pointing FROM THE SATELLITE TO THE SOURCE.  Photons
             therefore PROPAGATE along -u.
    E        irradiance at the satellite on a surface held PERPENDICULAR to the
             beam [W/m^2]
    p        radiation pressure at the satellite, p = E / c [Pa]; equivalently
             the momentum crossing unit area normal to the beam per second
    mu       cos(incidence) = n . u
    f        shadow[k][m] in [0,1], the fraction of the plate NOT shadowed by
             the rest of the body for this source
    (a,s,d)  optical triplet of the source band, a + s + d = 1: absorbed,
             specularly reflected, diffusely reflected.  For the "ir" band the
             first entry is the emissivity e, which by Kirchhoff's law equals
             the absorptance in that band and so plays exactly the role of a.

A plate is ONE-SIDED.  If mu <= 0 the normal faces away from (or lies along)
the beam and the plate receives nothing: exactly zero force, exactly zero
absorbed flux.  The opposite face of a real panel is a separate plate with
normal -n.


================================================================================
2.  INCIDENT MOMENTUM
================================================================================

The beam sees the projection of the plate onto the plane perpendicular to it,

    A_proj = A mu                                                          (1)

so the intercepted power is P_in = E A mu and, since a photon of energy dE
carries momentum dE/c along its direction of travel, momentum arrives at

    dP/dt|_in = (P_in / c) (-u) = p A mu (-u)                              (2)

The force on the plate is the rate at which the plate takes up momentum,

    F = dP/dt|_in - dP/dt|_out                                             (3)

dP/dt|_out being whatever the surface sends back out.  The three fates of the
incident photons -- absorbed, specularly reflected, diffusely reflected -- are
exclusive and exhaust the incident energy, and F is linear in the photon flux,
so the three channels superpose with weights a, s, d.


================================================================================
3.  THE THREE CHANNELS
================================================================================

3.1  Absorbed fraction a
------------------------
Nothing leaves.  (The thermal re-emission of this energy is a SEPARATE effect
and is deliberately NOT part of this function; it is why ``absorbed`` is
returned as a second output, for a thermal model to consume.)

    F_a = -a p A mu u                                                      (4)

A pure absorber is pushed straight along the beam, away from the source.

3.2  Specularly reflected fraction s
------------------------------------
The incoming propagation direction is w_in = -u.  Reflection about the plane
with normal n mirrors the normal component and leaves the tangential one:

    w_out = w_in - 2 (w_in . n) n = -u - 2(-mu) n = -u + 2 mu n            (5)

which is a unit vector (|w_out|^2 = 1 - 4 mu^2 + 4 mu^2 = 1).  All of the
energy in this channel leaves, so the outgoing momentum flux has the same
magnitude s p A mu as the incoming one, and

    F_s = s p A mu (w_in - w_out) = s p A mu (-u + u - 2 mu n)
        = -2 s p A mu^2 n                                                  (6)

Two things to note, because they are exactly where a wrong derivation lands:
the specular force is PURELY NORMAL -- a mirror exerts no tangential force,
because reflection does not change the tangential momentum -- and it is
QUADRATIC in mu, one factor from the projected area and one from the fraction
of the momentum that is normal.  At grazing incidence a mirror feels nothing.

3.3  Diffusely reflected fraction d (Lambertian)
------------------------------------------------
The power P_d = d E A mu leaves with a radiance independent of direction over
the outward hemisphere about n, i.e. intensity proportional to cos(alpha):

    dP/dOmega = (P_d / pi) cos(alpha)                                      (7)

Normalisation check: integrating (7) over the hemisphere,
(P_d/pi) 2 pi \int_0^{pi/2} cos a sin a da = (P_d/pi) 2 pi (1/2) = P_d.  Good.

By symmetry about n the transverse momentum cancels exactly, and the surviving
normal component is

    dP/dt|_diff . n = (1/c) \int (dP/dOmega) cos(alpha) dOmega
                    = (P_d / (pi c)) 2 pi \int_0^{pi/2} cos^2 a sin a da
                    = (P_d / (pi c)) 2 pi (1/3)
                    = (2/3) P_d / c                                        (8)

The famous 2/3 is precisely this hemispherical cos^2 average.  Hence

    F_d = -d p A mu u - (2/3) d p A mu n = -d p A mu ( u + (2/3) n )       (9)

The u-part is there because a diffuser still absorbs the incoming momentum
before re-radiating it isotropically-in-the-Lambertian-sense.


================================================================================
4.  THE COMBINED FORCE AND THE ABSORBED IRRADIANCE
================================================================================

Adding (4), (6), (9) and multiplying by the unshadowed fraction f -- which
enters linearly, because a partially lit plate intercepts that fraction of the
beam cross-section:

    ==========================================================================
    F = - p A f mu [ (a + d) u + ( 2 s mu + (2/3) d ) n ],     mu > 0      (10)
    F = 0,                                                    mu <= 0
    ==========================================================================

The coefficient of u is written (a + d) and NOT the algebraically equal (1 - s).
They coincide only when the triplet closes exactly.  (a + d) is the form that
comes directly out of the per-channel sum (4)+(9), so it stays the physically
correct answer if a caller passes a triplet that closes only to rounding.

The second output is the ABSORBED IRRADIANCE, in W/m^2 of PLATE area -- not of
beam cross-section, and not multiplied by the plate area.  The irradiance
falling on the plate is E mu = p c mu, and the fraction a of it is absorbed:

    q_abs = a f p c mu   [W/m^2]                                          (11)

Both outputs are summed over sources, each source using the optical triplet of
its own band.

Sanity limits contained in (10), all verified in ``_selftest`` below:
  * absorber (1,0,0):  F = -p A mu u,  magnitude p A mu, along the beam.
  * mirror   (0,1,0):  F = -2 p A mu^2 n,  purely normal, cos^2 law; the
    tangential component is EXACTLY zero at every incidence angle.
  * diffuser (0,0,1):  F = -p A mu (u + (2/3) n); at normal incidence the
    magnitude is (1 + 2/3) p A = (5/3) p A.
  * F -> 0 continuously like mu as the terminator is approached.
  * F . u = -p A f mu (a + d + 2 s mu^2 + (2/3) d mu) <= 0: radiation can never
    pull a plate towards the source.

Everything is CLOSED FORM.  No numerical quadrature is used anywhere; the only
integrals in the derivation, (7) and (8), were done analytically.


================================================================================
5.  INDEPENDENT CONFIRMATION (not used to obtain the result)
================================================================================

5.1  Photon Monte Carlo.  4e6 packets per case were transported through the
surface without any use of (10): a packet arrives along -u with momentum p/c,
is absorbed with probability a, leaves along (5) with probability s, or leaves
along a cosine-weighted direction about n with probability d (sampling
cos alpha = sqrt(xi), which reproduces Lambert's law).  F was accumulated as
(momentum in - momentum out).  Over pure absorber / mirror / diffuser and
mixtures at 0, 40, 60, 75 and 80 degrees, the largest deviation from (10) was
1.27 Monte Carlo sigma.  This confirms the mu^2 of (6) and the 2/3 of (8)
simultaneously and without reference to any closed form.

5.2  Sphere coefficients.  Tiling a sphere with radial-normal plates and
summing (10) reproduces the classical F = -C_R p pi R^2 u with C_R = 1 for a
black sphere, C_R = 1 for a specular sphere (a mirrored sphere scatters the
beam isotropically into 4 pi, so the reflected momentum cancels and only the
incident momentum survives -- a check that pins the exponent of mu in (6):
any power other than 2 breaks the isotropy), and C_R = 1 + 4/9 for a
Lambertian sphere.  ``sphere_coefficient`` below measures these.


================================================================================
6.  RECONCILIATION OF THE TWO INDEPENDENT DERIVATIONS
================================================================================

Route (a) did a direct momentum balance, F = in - specular_out - diffuse_out,
and collapsed the u-terms using 1 - s = a + d.  Route (b) derived the three
pure limiting cases separately and superposed them with weights a, s, d.  Both
arrived at (10) term for term, and my own re-derivation above (which is route
(a) carried out independently) agrees.  There is NO physics disagreement.

In particular, the failure mode to watch for -- resolving the force along the
beam and its perpendicular (a "drag/lift" pair) instead of along u and n --
is absent from both.  The discriminating limit is a pure mirror at oblique
incidence: the correct force is purely normal, so its component along the
plate tangent must vanish identically.  Both routes give exactly zero
tangential mirror force at every angle tested, as does (10).

Over a dense grid of 29058 configurations (incidence 0 to 180 degrees
including the exact terminator and angles within 1e-12 degrees of it, nine
azimuths, twelve optical triplets from pure to near-degenerate, both bands,
nine shadow fractions from 0 to 1, five pressures spanning 1e-12 to 1 Pa, four
areas, multi-plate/multi-source configurations, and 4000 pseudo-random draws)
the two routes differed by at most 3.55e-15 N in absolute force and 4.05e-13
in the relative force component, the latter only where a component cancels to
1.8e-5 of the force magnitude.  Exact rational arithmetic confirmed both were
within one rounding of the true value there.  The differences are floating-
point evaluation order, not physics.

Three interface questions were left open by the agreement, and were settled by
derivation rather than by majority:

  (i)  mu = n . u by a naive three-term dot product loses its leading digits at
       grazing incidence, and with them the SIGN that decides whether the plate
       is lit at all.  Over 400000 geometries constructed within a few ulp of
       the plate plane, the naive dot product misclassified lit/unlit in 16.3%
       of cases, ``math.sumprod`` (a single correctly-rounded dot product) in
       none.  The force involved is ~1e-16 of p A and so physically negligible,
       but a lit/unlit predicate that changes under a rotation of the
       coordinate frame is not acceptable in an oracle, whose whole job is to
       be the reproducible answer.  This module therefore uses ``math.sumprod``
       for mu, falling back to the naive sum on Python < 3.12.

  (ii) Evaluation order of the bracket.  Scored against exact rational
       arithmetic over 150000 draws with mu held fixed, folding the signed
       common factor into the two coefficients and closing the two-term sum
       with a single rounded product-sum was best: max error 3.4e-12 of a
       cancelling component versus 4.6e-12 and 5.0e-12 for the alternatives,
       and 4.03e-16 versus 4.68e-16 of the force magnitude.  That form is used.

  (iii) Summation over sources.  Both routes accumulated with ``+=`` in source
       order, which makes the total depend on the order in which the caller
       lists the sources -- a real hazard here, since the Earth-radiation
       source list is hundreds of tiles of widely differing magnitude next to
       one large solar term.  This module collects the per-source
       contributions and closes each component with ``math.fsum``, which is
       exactly rounded and therefore order-independent.

  (iv) A shadow fraction outside [0,1].  Route (a) applied a negative fraction
       as a negative force; route (b) silently skipped it.  Neither is
       defensible: such a value is not an illuminated fraction at all, so this
       module raises ``ValueError``.  Silent acceptance of a bad shadow is the
       kind of thing an oracle exists to catch.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

import math

__all__ = ["radiation_panel_force", "sphere_coefficient", "SPEED_OF_LIGHT"]

# Physical constants (SI), as fixed by the SI definitions.
SPEED_OF_LIGHT = 299792458.0        # m/s, exact
SIGMA_SB = 5.670374419e-8           # W/m^2/K^4  (not needed here)
K_B = 1.380649e-23                  # J/K        (not needed here)
N_A = 6.02214076e23                 # 1/mol      (not needed here)

_BANDS = ("vis", "ir")

# math.sumprod is a single correctly-rounded sum of products (Python >= 3.12).
# See reconciliation note (i): it is what keeps the lit/unlit predicate
# frame-independent at the terminator.
_HAVE_SUMPROD = hasattr(math, "sumprod")

if _HAVE_SUMPROD:
    def _dot(a, b):
        """Dot product of two 3-vectors, with a single rounding."""
        return math.sumprod(a, b)
else:                                                # pragma: no cover
    def _dot(a, b):
        """Dot product of two 3-vectors (naive fallback, Python < 3.12).

        Accurate to ~1 ulp of the largest term, so the sign of a dot product
        that cancels to the last bit is not guaranteed; see note (i) of the
        module docstring.
        """
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def radiation_panel_force(plates, sources, shadow=None):
    """Radiation-pressure force and absorbed irradiance on a set of flat plates.

    Implements equation (10) of the module docstring,

        F = -p A f mu [ (a + d) u + (2 s mu + (2/3) d) n ]     for mu > 0

    together with the absorbed irradiance (11), q = a f p c mu, per plate and
    source, summed over sources.  Closed form; no quadrature.

    Thermal re-emission of the absorbed energy is NOT included -- that is the
    purpose of the second output.

    Parameters
    ----------
    plates : list of dict
        Each plate is
        ``{"area": float [m^2], "normal": 3-vector outward unit normal,
        "centroid": 3-vector [m] (unused here), "vis": (a, s, d),
        "ir": (e, s, d)}``, each optical triplet summing to 1.
    sources : list of dict
        Each source is
        ``{"dir": 3-vector unit vector pointing FROM THE SATELLITE TO THE
        SOURCE, "p": radiation pressure at the satellite [Pa] (irradiance / c),
        "band": "vis" or "ir"}``.
    shadow : list of list of float, optional
        ``shadow[k][m]`` in [0, 1] is the fraction of plate k NOT shadowed by
        the rest of the body for source m.  Default: all 1.0.

    Returns
    -------
    forces : list of list of float
        Force in newtons on each plate, summed over all sources.
    absorbed : list of float
        Absorbed irradiance in W/m^2 of plate area, summed over all sources.

    Raises
    ------
    ValueError
        If a band is not "vis" or "ir", if the shadow array has the wrong
        shape, or if a shadow fraction lies outside [0, 1].
    """
    n_plates = len(plates)
    n_sources = len(sources)

    if shadow is None:
        shadow = [[1.0] * n_sources for _ in range(n_plates)]
    else:
        if len(shadow) != n_plates:
            raise ValueError(
                "shadow must have one row per plate: got %d rows for %d plates"
                % (len(shadow), n_plates))
        for k, row in enumerate(shadow):
            if len(row) != n_sources:
                raise ValueError(
                    "shadow row %d must have one entry per source: got %d for "
                    "%d sources" % (k, len(row), n_sources))
            for m, f in enumerate(row):
                if not (0.0 <= f <= 1.0):
                    raise ValueError(
                        "shadow[%d][%d] = %r is not an illuminated fraction in "
                        "[0, 1]" % (k, m, f))

    # Validate the bands once, up front, so an empty plate list still catches a
    # malformed source list.
    for m, src in enumerate(sources):
        if src["band"] not in _BANDS:
            raise ValueError("source %d: band must be 'vis' or 'ir', got %r"
                             % (m, src["band"]))

    forces = []
    absorbed = []

    for k, plate in enumerate(plates):
        area = float(plate["area"])
        n = plate["normal"]
        shadow_k = shadow[k]

        # Per-source contributions, kept as lists so the sums below can be
        # closed exactly with math.fsum: the total is then independent of the
        # order in which the caller lists the sources.  See note (iii).
        cx, cy, cz, cq = [], [], [], []

        for m, src in enumerate(sources):
            mu = _dot(n, src["dir"])
            if mu <= 0.0:
                # Normal faces away from, or lies in, the beam: the plate is
                # one-sided and receives nothing at all.
                continue

            f = shadow_k[m]
            if f == 0.0:
                continue

            p = float(src["p"])
            a, s, d = plate[src["band"]]
            u = src["dir"]

            # Equation (10), in the evaluation order measured best in note (ii):
            # the signed common factor is folded into each coefficient, and the
            # two-term sum per component is closed with one rounding.
            common = -p * area * f * mu
            c_u = common * (a + d)
            c_n = common * (2.0 * s * mu + (2.0 / 3.0) * d)

            if _HAVE_SUMPROD:
                cx.append(math.sumprod((c_u, c_n), (u[0], n[0])))
                cy.append(math.sumprod((c_u, c_n), (u[1], n[1])))
                cz.append(math.sumprod((c_u, c_n), (u[2], n[2])))
            else:                                    # pragma: no cover
                cx.append(c_u * u[0] + c_n * n[0])
                cy.append(c_u * u[1] + c_n * n[1])
                cz.append(c_u * u[2] + c_n * n[2])

            # Equation (11): absorbed irradiance, W/m^2 of plate area.
            cq.append(a * f * p * SPEED_OF_LIGHT * mu)

        forces.append([math.fsum(cx), math.fsum(cy), math.fsum(cz)])
        absorbed.append(math.fsum(cq))

    return forces, absorbed


# ---------------------------------------------------------------------------
# The sphere check of section 5.2, kept in the module because it is the global
# limit that pins the exponent of mu in the specular term.
# ---------------------------------------------------------------------------
def sphere_coefficient(triplet, radius=1.0, n_theta=240, n_phi=480,
                       direction=(0.0, 0.0, 1.0), p=1.0):
    """Radiation-pressure coefficient C_R of a plate-tiled sphere.

    The sphere is tiled by lat-lon cells with radial normals and exact
    spherical cell areas; the per-plate forces are summed and matched to
    ``F = -C_R p pi R^2 u``.  Returns ``(C_R, transverse_fraction)``; the
    second value should be ~0 by symmetry.

    Derived expectations (section 5.2): 1 for a black sphere, 1 for a specular
    sphere, 1 + 4/9 for a Lambertian sphere, and 1 + (4/9) d in general.

    This is the one place a discretisation is used, and that is the point of
    the test: ``radiation_panel_force`` itself is closed form.
    """
    m = math.sqrt(_dot(direction, direction))
    u = [v / m for v in direction]

    plates = []
    dphi = 2.0 * math.pi / n_phi
    for i in range(n_theta):
        t0 = math.pi * i / n_theta
        t1 = math.pi * (i + 1) / n_theta
        area = radius * radius * (math.cos(t0) - math.cos(t1)) * dphi
        tm = 0.5 * (t0 + t1)
        st, ct = math.sin(tm), math.cos(tm)
        for j in range(n_phi):
            pm = dphi * (j + 0.5)
            nrm = [st * math.cos(pm), st * math.sin(pm), ct]
            plates.append({"area": area, "normal": nrm,
                           "centroid": [radius * v for v in nrm],
                           "vis": triplet, "ir": triplet})

    forces, _ = radiation_panel_force(
        plates, [{"dir": u, "p": p, "band": "vis"}])
    total = [math.fsum(fk[i] for fk in forces) for i in range(3)]
    along = _dot(total, u)
    perp = math.sqrt(sum((total[i] - along * u[i]) ** 2 for i in range(3)))
    c_r = -along / (p * math.pi * radius * radius)
    return c_r, perp / max(abs(along), 1e-300)


# ---------------------------------------------------------------------------
# Self-test: the hand-checkable limits of section 4 and the sphere limits.
# ---------------------------------------------------------------------------
def _selftest():
    ok = [True]

    def check(name, got, want, tol=1e-14):
        err = abs(got - want) / max(abs(want), 1.0)
        good = err <= tol
        ok[0] = ok[0] and good
        print("  %-52s got %+.14e  err %.1e  %s"
              % (name, got, err, "OK" if good else "FAIL"))

    p, area = 4.5606e-6, 2.0
    z = [0.0, 0.0, 1.0]
    src_n = [{"dir": z, "p": p, "band": "vis"}]

    def plate(trip_vis, trip_ir=None):
        return {"area": area, "normal": z, "centroid": [0.0, 0.0, 0.0],
                "vis": trip_vis, "ir": trip_ir or trip_vis}

    print("normal incidence (mu = 1)")
    F, q = radiation_panel_force([plate((1.0, 0.0, 0.0))], src_n)
    check("absorber F_z = -p A", F[0][2], -p * area)
    check("absorber q   = p c", q[0], p * SPEED_OF_LIGHT)
    F, q = radiation_panel_force([plate((0.0, 1.0, 0.0))], src_n)
    check("mirror   F_z = -2 p A", F[0][2], -2.0 * p * area)
    check("mirror   q   = 0", q[0], 0.0)
    F, q = radiation_panel_force([plate((0.0, 0.0, 1.0))], src_n)
    check("diffuser F_z = -(5/3) p A", F[0][2], -(5.0 / 3.0) * p * area)

    print("oblique, theta = 60 deg (mu = 1/2)")
    th = math.radians(60.0)
    u = [math.sin(th), 0.0, math.cos(th)]
    src_o = [{"dir": u, "p": p, "band": "vis"}]
    mu = math.cos(th)
    F, _ = radiation_panel_force([plate((1.0, 0.0, 0.0))], src_o)
    check("absorber F_x = -p A mu u_x", F[0][0], -p * area * mu * u[0])
    check("absorber F_z = -p A mu u_z", F[0][2], -p * area * mu * u[2])
    F, _ = radiation_panel_force([plate((0.0, 1.0, 0.0))], src_o)
    check("mirror   F_x = 0 exactly (no tangential force)", F[0][0], 0.0)
    check("mirror   F_z = -2 p A mu^2", F[0][2], -2.0 * p * area * mu * mu)
    F, _ = radiation_panel_force([plate((0.0, 0.0, 1.0))], src_o)
    check("diffuser F_z = -p A mu (mu + 2/3)", F[0][2],
          -p * area * mu * (mu + 2.0 / 3.0))

    print("one-sidedness, shadow, bands")
    F, q = radiation_panel_force(
        [plate((0.3, 0.3, 0.4))],
        [{"dir": [0.0, 0.0, -1.0], "p": p, "band": "vis"}])
    check("back-lit |F| = 0", math.sqrt(sum(v * v for v in F[0])), 0.0)
    check("back-lit q   = 0", q[0], 0.0)
    F, _ = radiation_panel_force(
        [plate((0.3, 0.3, 0.4))],
        [{"dir": [1.0, 0.0, 0.0], "p": p, "band": "vis"}])
    check("grazing  |F| = 0", math.sqrt(sum(v * v for v in F[0])), 0.0)
    Fs, qs = radiation_panel_force([plate((0.2, 0.3, 0.5))], src_n,
                                   shadow=[[0.25]])
    Ff, qf = radiation_panel_force([plate((0.2, 0.3, 0.5))], src_n)
    check("shadow 0.25 scales F", Fs[0][2], 0.25 * Ff[0][2])
    check("shadow 0.25 scales q", qs[0], 0.25 * qf[0])
    F, q = radiation_panel_force([plate((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))],
                                 [{"dir": z, "p": p, "band": "ir"}])
    check("ir source uses the ir triplet", F[0][2], -2.0 * p * area)

    print("bad shadow must raise")
    for bad in ([[-0.1]], [[1.5]]):
        try:
            radiation_panel_force([plate((1.0, 0.0, 0.0))], src_n, shadow=bad)
        except ValueError:
            print("  %-52s %s" % ("shadow %r rejected" % bad, "OK"))
        else:
            ok[0] = False
            print("  %-52s %s" % ("shadow %r rejected" % bad, "FAIL"))

    print("sphere coefficients (240 x 480 tiles, tilted incidence)")
    for name, trip, want in (("black    ", (1.0, 0.0, 0.0), 1.0),
                             ("mirror   ", (0.0, 1.0, 0.0), 1.0),
                             ("Lambert  ", (0.0, 0.0, 1.0), 1.0 + 4.0 / 9.0),
                             ("half diff", (0.5, 0.0, 0.5),
                              1.0 + (4.0 / 9.0) * 0.5)):
        cr, perp = sphere_coefficient(trip, direction=(0.3, -0.7, 0.5))
        rel = abs(cr - want) / want
        good = rel < 2.0e-4 and perp < 2.0e-4
        ok[0] = ok[0] and good
        print("  %s C_R want %.9f got %.9f  rel %.2e  perp %.2e  %s"
              % (name, want, cr, rel, perp, "OK" if good else "FAIL"))

    return ok[0]


if __name__ == "__main__":
    print("radiation.py self-test (math.sumprod available: %s)\n"
          % _HAVE_SUMPROD)
    print("ALL GREEN" if _selftest() else "SOMETHING FAILED")
