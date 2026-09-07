"""Thermal recoil force on grey Lambertian panels -- independent oracle.

Standard library only (``math``).  Written from the physics, not transcribed
from any implementation, so that it can serve as an oracle for the MATLAB
kernel ``NonGravForces.kernel.thermalRecoil``.

DERIVATION
==========

Setup
-----
A flat plate of area ``A`` has an OUTWARD unit normal ``n`` (pointing away
from the spacecraft body, into space).  Its outer face is at temperature
``T`` and has hemispherical emissivity ``eps``.  It radiates into the outward
hemisphere only; the body blocks the other half.

Step 1 -- total emitted flux (Stefan-Boltzmann)
-----------------------------------------------
A grey emitter radiates a total hemispherical exitance

    M = eps * sigma * T**4                            [W / m^2]                (1)

with ``sigma = 5.670374419e-8 W m^-2 K^-4``.  The emitted power is
``P = A * M``.

Step 2 -- Lambertian angular distribution
------------------------------------------
"Lambertian" means the radiance ``L`` [W m^-2 sr^-1] leaving the surface is
the same in every outgoing direction; the *projected* area seen from a
direction at polar angle ``theta`` from ``n`` shrinks as ``cos(theta)``, so
the power per unit solid angle follows the cosine law,

    dP / dOmega = A * L * cos(theta).                                          (2)

``L`` is fixed by requiring that (2) integrate to the exitance (1) over the
outward hemisphere.  With ``dOmega = sin(theta) dtheta dphi``,

    M = Int_0^{2pi} dphi Int_0^{pi/2} L cos(theta) sin(theta) dtheta
      = 2*pi * L * [ sin^2(theta) / 2 ]_0^{pi/2}
      = pi * L,

so

    L = M / pi = eps * sigma * T**4 / pi.                                      (3)

This is the familiar factor ``1/pi`` between exitance and radiance of a
Lambertian surface.

Step 3 -- momentum carried away
--------------------------------
A photon of energy ``E`` carries linear momentum ``E / c`` along its
direction of travel ``u``.  The rate at which momentum is carried away by the
emitted field is therefore the direction-weighted power integral

    dp/dt = (A / c) * Int_hemisphere L * cos(theta) * u dOmega.                (4)

Decompose ``u`` into a component along ``n`` and a component in the tangent
plane.  Integrating over the full azimuth ``phi`` in [0, 2pi) kills the
tangential part by symmetry, leaving only ``u . n = cos(theta)``:

    dp/dt = (A * L / c) * n * Int_0^{2pi} dphi Int_0^{pi/2} cos^2(theta) sin(theta) dtheta.

The angular integral is elementary (closed form -- no quadrature is used
anywhere in this module):

    Int_0^{2pi} dphi = 2*pi,
    Int_0^{pi/2} cos^2(theta) sin(theta) dtheta = [ -cos^3(theta) / 3 ]_0^{pi/2} = 1/3,

hence

    Int_hemisphere cos^2(theta) dOmega = 2*pi / 3.                             (5)

Substituting (3) and (5) into (4):

    dp/dt = (A / c) * (eps*sigma*T^4 / pi) * (2*pi/3) * n
          = (2/3) * A * eps * sigma * T**4 / c * n.                            (6)

Step 4 -- the recoil
---------------------
Newton's third law: the momentum given to the radiation field along ``+n`` is
taken from the plate, so the force ON THE PLATE is the negative of (6):

    F = - (2/3) * eps * sigma * T**4 * A / c * n.                              (7)

DIRECTION: anti-parallel to the outward normal.  A hot panel is pushed
*backwards*, away from the direction its photons leave in -- exactly as a
rocket is pushed opposite its exhaust.

HEMISPHERICAL FACTOR: 2/3, and it is the product of two separate pieces that
must not be confused with one another:

    * ``1/pi``   from normalising the Lambertian radiance to the exitance,
                 eq. (3) -- one power of cos in the emission law;
    * ``2*pi/3`` from projecting each ray's momentum onto ``n``, eq. (5) --
                 the second power of cos.

Only the product ``(1/pi) * (2*pi/3) = 2/3`` is physical.  Two neighbouring
wrong answers, for orientation:

    * ``1``   would be a perfectly collimated beam fired along ``n``
              (every photon contributing its full momentum);
    * ``1/2`` would be an isotropic-over-the-hemisphere emitter
              (radiance flat per unit solid angle rather than per unit
              projected area, i.e. mean of cos(theta) over the hemisphere).

The true Lambertian value 2/3 lies between them, as it must: a cosine-law
emitter is more forward-peaked than isotropic but less than collimated.

Sanity limits
-------------
    * ``T = 0``          -> F = 0 (nothing is radiated).
    * ``eps = 0``        -> F = 0 (a perfect reflector emits nothing).
    * ``A = 0``          -> F = 0.
    * ``T -> 2T``        -> |F| x 16 (quartic in temperature).
    * Two identical, equally hot plates with opposite normals give exactly
      cancelling forces; a closed isothermal box has zero net thermal recoil.

Scale check for GRACE-FO: one 0.9 m^2 panel at 290 K with eps = 0.8 gives
|F| = 2/3 * 0.8 * 5.670e-8 * 290^4 * 0.9 / c = 6.42e-7 N, i.e. 1.07e-9 m/s^2
on a 601.2 kg spacecraft.  The net TRP acceleration is set by the *imbalance*
between panels, not by any single one, which is why the published GRACE TRP
magnitude (~0.8 nm/s^2) is below this single-panel figure.

What is NOT in this model
-------------------------
Only the emissivity entry of the infrared triplet is used.  The specular and
diffuse entries of ``ir`` describe how the panel REFLECTS incident infrared
radiation, which is a separate (radiation-pressure) effect and no part of the
plate's own thermal emission.  Likewise the visible triplet plays no role
here, and the plate centroid matters only for torque, not for force.

The emission is taken to be grey (emissivity independent of wavelength) and
the outer face isothermal at ``T``; radiation from the inner face is assumed
to be exchanged inside the body and to contribute no net external recoil.

REFERENCES
----------
Rievers, B. and Laemmerzahl, C. (2011), "High precision thermal modeling of
    complex systems"; the 2/3 Lambertian recoil factor.
Woeske, F. (2019), PhD thesis, eq. (5.14):
    F_trp,k = -(2/3) * A_k * eps_k * sigma * T_k^4 / c * e_N,k.
Both agree with eq. (7) above.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

import math

__all__ = ["thermal_recoil", "SIGMA_SB", "C_LIGHT", "LAMBERT_RECOIL_FACTOR"]

# Physical constants (SI).
SIGMA_SB = 5.670374419e-8      # Stefan-Boltzmann constant, W m^-2 K^-4
C_LIGHT = 299792458.0          # speed of light in vacuum, m s^-1

# Hemispherical recoil factor for a Lambertian emitter, derived in the module
# docstring: (1/pi) * Int_hemisphere cos^2(theta) dOmega = (1/pi) * (2*pi/3).
LAMBERT_RECOIL_FACTOR = 2.0 / 3.0

_UNIT_NORMAL_TOL = 1.0e-8      # tolerance on |n| - 1 before renormalising


def _unit(vec, index):
    """Return ``vec`` as a unit 3-tuple of floats, validating it."""
    try:
        x, y, z = (float(component) for component in vec)
    except (TypeError, ValueError):
        raise ValueError(
            "plate %d: 'normal' must be a 3-vector of finite floats" % index
        )
    norm = math.sqrt(x * x + y * y + z * z)
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError(
            "plate %d: 'normal' must have non-zero finite length" % index
        )
    if abs(norm - 1.0) <= _UNIT_NORMAL_TOL:
        return (x, y, z)
    # Only the direction of n enters eq. (7); renormalise a slightly
    # off-unit normal rather than propagating the length into the force.
    return (x / norm, y / norm, z / norm)


def thermal_recoil(plates, temps):
    """Recoil force on each plate from its own thermal (infrared) emission.

    Each plate is treated as a grey Lambertian emitter radiating from its
    outward face only.  Per the derivation in the module docstring,

        F = -(2/3) * eps * sigma * T**4 * A / c * n                    [N]

    with ``eps`` the plate's infrared emissivity (the first entry of the
    ``ir`` triplet), ``sigma`` the Stefan-Boltzmann constant, ``A`` the plate
    area, ``c`` the speed of light and ``n`` the OUTWARD unit normal.  The
    force is anti-parallel to ``n``.

    Parameters
    ----------
    plates : list of dict
        Plate dicts with at least the keys ``"area"`` (m^2), ``"normal"``
        (outward unit 3-vector) and ``"ir"`` (triplet whose first entry is
        the emissivity).  ``"centroid"`` and ``"vis"``, if present, are
        ignored: they play no part in the force.
    temps : list of float
        Plate temperatures in kelvin, one per plate, each >= 0.

    Returns
    -------
    list of tuple
        One 3-vector of newtons per plate, in the same order and in the same
        frame as the plate normals (the body frame).

    Raises
    ------
    ValueError
        If ``len(temps) != len(plates)``, or a temperature is negative or
        non-finite, or an area is negative or non-finite, or an emissivity
        lies outside [0, 1], or a normal is degenerate.
    """
    plates = list(plates)
    temps = list(temps)
    if len(temps) != len(plates):
        raise ValueError(
            "temps has %d entries but plates has %d -> pass one temperature "
            "per plate" % (len(temps), len(plates))
        )

    forces = []
    for index, (plate, temperature) in enumerate(zip(plates, temps)):
        temperature = float(temperature)
        if not math.isfinite(temperature) or temperature < 0.0:
            raise ValueError(
                "plate %d: temperature must be finite and >= 0 K, got %r"
                % (index, temperature)
            )

        area = float(plate["area"])
        if not math.isfinite(area) or area < 0.0:
            raise ValueError(
                "plate %d: 'area' must be finite and >= 0 m^2, got %r"
                % (index, area)
            )

        emissivity = float(plate["ir"][0])
        if not math.isfinite(emissivity) or not 0.0 <= emissivity <= 1.0:
            raise ValueError(
                "plate %d: infrared emissivity ir[0] must lie in [0, 1], "
                "got %r" % (index, emissivity)
            )

        normal = _unit(plate["normal"], index)

        # Eq. (7).  Grouped so that the emitted power A*eps*sigma*T^4 [W] is
        # formed first, then converted to a momentum flux by 1/c.
        emitted_power_W = area * emissivity * SIGMA_SB * temperature ** 4
        magnitude_N = LAMBERT_RECOIL_FACTOR * emitted_power_W / C_LIGHT

        # The trailing "+ 0.0" turns an IEEE negative zero into +0.0 so that
        # a zero component prints and compares cleanly.
        forces.append(
            (
                -magnitude_N * normal[0] + 0.0,
                -magnitude_N * normal[1] + 0.0,
                -magnitude_N * normal[2] + 0.0,
            )
        )

    return forces
