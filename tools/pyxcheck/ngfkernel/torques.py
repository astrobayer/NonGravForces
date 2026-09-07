"""Body torques on a rigid satellite: gravity gradient and magnetic dipole.

Independent standard-library oracle for the NonGravForces kernel.  Written from
the equations derived below, not transcribed from the MATLAB implementation.

Python 3, standard library only (``math``).  SI units, angles in radians.
A 3-vector is a tuple/list of three floats; a 3x3 matrix is a list of three rows.


DERIVATION 1 -- GRAVITY-GRADIENT TORQUE
=======================================

Setting.  A rigid body of total mass ``M`` orbits a point mass ``GM`` at the
origin.  Let ``R`` be the geocentric position of the body's centre of mass and
let ``rho`` be the position of a mass element ``dm`` relative to that centre of
mass, so the geocentric position of the element is ``R + rho``.  Both vectors
are expressed in the body frame; the derivation is frame-agnostic because the
torque is a vector identity, so working in the body frame costs nothing.

Force on the element (Newton, point-mass field):

    dF = -GM * (R + rho) / |R + rho|**3 * dm

Torque about the centre of mass:

    T = INT rho x dF
      = -GM * INT rho x (R + rho) / |R + rho|**3 dm
      = -GM * INT (rho x R) / |R + rho|**3 dm            [rho x rho = 0]

Expand the denominator.  With ``R = |R|``,

    |R + rho|**-3 = R**-3 * (1 + 2*(R.rho)/R**2 + rho**2/R**2)**-1.5
                  = R**-3 * (1 - 3*(R.rho)/R**2 + O(rho**2/R**2))

THE ASSUMPTION: |rho| << R.  The body's dimensions are negligible against the
orbit radius, and the binomial series is truncated after the linear term in
rho/R.  For GRACE-FO (half-extent ~2.5 m, R ~ 6.8e6 m) rho/R ~ 4e-7, so the
neglected quadratic term is ~1e-6 relative -- utterly negligible.  Two further
assumptions ride along: the field is that of a point mass (no J2 or higher
harmonics, which contribute their own gradient), and the body is rigid so that
the inertia tensor is constant in the body frame.

Substituting, the two surviving orders are

    T ~= -GM/R**3 * INT (rho x R) dm  +  3*GM/R**5 * INT (R.rho) (rho x R) dm

The first integral vanishes: ``INT rho dm = 0`` is the definition of the centre
of mass.  This is why the leading term of the gravity-gradient torque is the
*gradient* term -- a uniform field exerts no torque about the centre of mass.
The second gives

    T = (3*GM/R**5) * [ INT rho (rho.R) dm ] x R
      = (3*GM/R**5) * (J R) x R

where ``J = INT rho rho^T dm`` is the second-moment (covariance) tensor.

Relate ``J`` to the inertia tensor ``I = INT (|rho|**2 * Id - rho rho^T) dm``:

    I = tr(J) * Id - J          =>      J = tr(J) * Id - I
    tr(I) = 3*tr(J) - tr(J) = 2*tr(J)   =>      tr(J) = tr(I)/2

so ``J R = (tr(I)/2) * R - I R``.  Crossing with ``R`` kills the first piece
(``R x R = 0``):

    (J R) x R = -(I R) x R = R x (I R)

Hence the closed form -- no quadrature anywhere in this module:

    T_gg = (3*GM / |R|**5) * ( R x (I R) )                              (1)

or, with the nadir-pointing unit vector ``u = R/|R|`` and mean motion
``n**2 = GM/|R|**3``,

    T_gg = 3*n**2 * ( u x (I u) )                                       (1')

Sanity checks built into the result.

*  If ``R`` is along a principal axis of ``I`` then ``I R || R`` and the torque
   vanishes: principal axes pointing at nadir are equilibria.  Correct.
*  Restoring sense.  Take ``I = diag(Ix, Iy, Iz)`` and rotate the body by a
   small angle ``phi`` about ``z``; the nadir direction seen in body
   coordinates then sits at ``u = (cos a, sin a, 0)`` with ``a = -phi``.
   Equation (1') gives ``T_z = 3*n**2*(Iy - Ix)*sin(a)*cos(a)``, i.e.
   ``T_z = -3*n**2*(Iy - Ix)*sin(phi)*cos(phi)``.  With ``Ix < Iy`` this
   opposes ``phi``: the MINOR-inertia axis is driven to the local vertical.
   That is the textbook gravity-gradient stable configuration, so the sign and
   the factor 3 are confirmed independently of any reference.
*  ``T_gg`` is orthogonal to ``R`` by construction -- a point-mass gradient
   torque has no component along nadir.
*  Adding a multiple of the identity to ``I`` changes nothing, since
   ``R x (c*R) = 0``.  Only inertia *differences* produce torque, as they must.

Note that circularity of the orbit is NOT required: (1) holds instantaneously
for any ``R``.  See CONTRADICTIONS below.


DERIVATION 2 -- MAGNETIC TORQUE
===============================

A magnetic dipole ``m`` in an external flux density ``B`` has interaction
energy ``U = -m.B``.  For a rigid dipole of fixed magnitude the generalised
torque about an axis ``e`` is ``T.e = -dU/dtheta`` for a rotation ``theta``
about ``e``; with ``dm/dtheta = e x m`` this gives

    T.e = d(m.B)/dtheta = (e x m).B = e.(m x B)      for every e

    =>  T_mag = m x B                                                   (2)

Equivalently, from the Lorentz force on a current loop:
``T = INT r x (I dl x B)``, which for a small planar loop of current ``I`` and
vector area ``A`` reduces to ``(I*A) x B = m x B``.

Sanity check.  ``m`` along ``+z``, ``B`` along ``+x``: (2) gives ``T`` along
``+y``, and a right-handed rotation about ``+y`` carries ``z -> x``, i.e. it
turns ``m`` toward ``B``, lowering ``U = -m.B``.  Correct sign; the opposite
convention ``B x m`` would be anti-restoring.

Note (2) is exact for a point dipole in a *uniform* field.  A field gradient
adds a net force ``grad(m.B)`` and, through it, a torque about the centre of
mass of order ``|r_dipole - r_com| * |grad B|``; both are ignored here, which
is the standard treatment and is what the interface (a single ``b_body``
sample, no gradient) permits.


CONTRADICTIONS WITH THE REFERENCE DOCUMENTS
===========================================

None of substance; the derivation above and both documents agree on the
formulae.  Two remarks:

1.  Woeske Eq. (2.13) presents ``T_gg = 3*GM/|rb|**5 * (rb x I rb)`` as an
    approximation valid "for a simple spherical gravitational potential and a
    circular satellite orbit".  The derivation above shows the circular-orbit
    qualifier is unnecessary: (1) is the instantaneous point-mass gradient
    torque for ANY orbit, eccentric included, because nothing in the expansion
    referred to the orbit's shape or to the mean motion.  The genuine
    approximations are ``|rho| << R`` and the point-mass field.  This is an
    overstatement in the document, not an error in the formula.

2.  Woeske Eq. (2.12) gives the general form for an arbitrary field via the
    gravity-gradient (Jacobi) matrix ``Gg = dg/dr``, ``T = INT rho x Gg rho
    dm``.  Specialising that to ``Gg = -(GM/R**3)*(Id - 3*u u^T)`` reproduces
    (1) exactly, which is an independent confirmation of the factor 3 and the
    sign.  This module implements only the point-mass case, per the interface.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

import math

__all__ = ["body_torques", "gravity_gradient_torque", "magnetic_torque"]


# ---------------------------------------------------------------------------
# small vector / matrix helpers (standard library only)
# ---------------------------------------------------------------------------

def _vec3(v, name):
    """Validate and return a 3-vector as a tuple of floats."""
    try:
        a, b, c = v
    except (TypeError, ValueError):
        raise ValueError("%s must be a sequence of three numbers" % name)
    return (float(a), float(b), float(c))


def _mat3(m, name):
    """Validate and return a 3x3 matrix as a tuple of three row tuples."""
    try:
        r0, r1, r2 = m
    except (TypeError, ValueError):
        raise ValueError("%s must be three rows of three numbers" % name)
    return (_vec3(r0, name + " row 1"),
            _vec3(r1, name + " row 2"),
            _vec3(r2, name + " row 3"))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _matvec(m, v):
    return (m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2])


def _norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


# ---------------------------------------------------------------------------
# the two torques
# ---------------------------------------------------------------------------

def gravity_gradient_torque(inertia, r_body, gm):
    """Gravity-gradient torque, Eq. (1): 3*GM/|r|**5 * (r x I r).

    inertia  3x3 inertia tensor in the body frame, kg m^2 (list of three rows)
    r_body   geocentric position expressed in the BODY frame, m
    gm       Earth gravitational parameter, m^3/s^2

    Returns the torque about the centre of mass in the body frame, N m.
    """
    inr = _mat3(inertia, "inertia")
    r = _vec3(r_body, "r_body")
    gm = float(gm)

    rmag = _norm(r)
    if rmag <= 0.0:
        raise ValueError("r_body must be a non-zero vector "
                         "-> the gravity-gradient torque is singular at r = 0")
    if gm <= 0.0:
        raise ValueError("gm must be positive -> pass the Earth GM in m^3/s^2")

    # Scale by 1/|r| three times and once more for the two extra powers, so the
    # factor never overflows for large |r| and stays well conditioned.
    scale = 3.0 * gm / (rmag ** 5)
    ir = _matvec(inr, r)
    t = _cross(r, ir)
    return (scale * t[0], scale * t[1], scale * t[2])


def magnetic_torque(dipole, b_body):
    """Magnetic torque, Eq. (2): m x B.

    dipole   satellite magnetic dipole moment in the body frame, A m^2
    b_body   magnetic flux density at the satellite in the body frame, T

    Returns the torque in the body frame, N m.
    """
    m = _vec3(dipole, "dipole")
    b = _vec3(b_body, "b_body")
    return _cross(m, b)


def body_torques(inertia, r_body, gm,
                 b_body=None, dipole=None,
                 gravity_gradient=False, magnetic=False):
    """Sum of the selected body torques in the body frame, N m.

    inertia           3x3 inertia tensor, body frame, kg m^2
    r_body            geocentric position expressed in the body frame, m
    gm                Earth gravitational parameter, m^3/s^2
    b_body            magnetic flux density in the body frame, T
    dipole            satellite magnetic dipole moment in the body frame, A m^2
    gravity_gradient  include 3*GM/|r|**5 * (r x I r)
    magnetic          include dipole x b_body

    With both flags false the result is exactly (0.0, 0.0, 0.0); the two terms
    are independent and simply add, so any subset is a valid request.
    """
    tx = ty = tz = 0.0

    if gravity_gradient:
        t = gravity_gradient_torque(inertia, r_body, gm)
        tx += t[0]
        ty += t[1]
        tz += t[2]

    if magnetic:
        if b_body is None or dipole is None:
            raise ValueError("magnetic=True needs both b_body and dipole "
                             "-> pass the field in T and the moment in A m^2")
        t = magnetic_torque(dipole, b_body)
        tx += t[0]
        ty += t[1]
        tz += t[2]

    return (tx, ty, tz)


# ---------------------------------------------------------------------------
# self-check: cases whose answers are predicted by hand in the docstring above
# ---------------------------------------------------------------------------

def _selfcheck():
    GM = 3.986004418e14
    inertia = [[1000.0, 0.0, 0.0],
               [0.0, 2000.0, 0.0],
               [0.0, 0.0, 3000.0]]

    # --- hand case A ------------------------------------------------------
    # r at 45 deg between body x and y, |r| = 7.0e6 m.
    #   u = (1,1,0)/sqrt(2);  I u = (1000, 2000, 0)/sqrt(2)
    #   u x (I u) = (0, 0, 2000/2 - 1000/2) = (0, 0, 500)
    #   T = 3*GM/|r|^3 * (0,0,500)
    #     = 3*3.986004418e14 / 3.43e20 * 500
    #     = 3.486301240e-6 * 500 = 1.743150620e-3 N m about +z
    R = 7.0e6
    r = (R / math.sqrt(2.0), R / math.sqrt(2.0), 0.0)
    t = body_torques(inertia, r, GM, gravity_gradient=True)
    expect = (0.0, 0.0, 3.0 * GM / R ** 3 * 500.0)
    assert abs(t[0]) < 1e-15 and abs(t[1]) < 1e-15
    assert abs(t[2] - expect[2]) < 1e-15 * abs(expect[2]) + 1e-18
    assert abs(t[2] - 1.743150620e-3) < 1e-11, t[2]

    # --- hand case B: magnetic -------------------------------------------
    #   m = (0.2, 0, 0) A m^2,  B = (0, 3e-5, 0) T
    #   m x B = (0, 0, 0.2*3e-5) = (0, 0, 6e-6) N m
    t = body_torques(inertia, r, GM, b_body=(0.0, 3e-5, 0.0),
                     dipole=(0.2, 0.0, 0.0), magnetic=True)
    assert t == (0.0, 0.0, 6e-6), t

    # --- both on: the sum -------------------------------------------------
    t = body_torques(inertia, r, GM, b_body=(0.0, 3e-5, 0.0),
                     dipole=(0.2, 0.0, 0.0),
                     gravity_gradient=True, magnetic=True)
    assert abs(t[2] - (1.743150620e-3 + 6e-6)) < 1e-11, t

    # --- neither on -------------------------------------------------------
    assert body_torques(inertia, r, GM) == (0.0, 0.0, 0.0)

    # --- principal axis at nadir gives zero torque ------------------------
    for axis in ((R, 0.0, 0.0), (0.0, R, 0.0), (0.0, 0.0, R)):
        t = body_torques(inertia, axis, GM, gravity_gradient=True)
        assert max(abs(c) for c in t) < 1e-18, (axis, t)

    # --- spherical inertia gives zero torque for any r --------------------
    sph = [[900.0, 0.0, 0.0], [0.0, 900.0, 0.0], [0.0, 0.0, 900.0]]
    t = body_torques(sph, (1.0e6, -2.0e6, 6.5e6), GM, gravity_gradient=True)
    assert max(abs(c) for c in t) < 1e-18, t

    # --- torque is perpendicular to r ------------------------------------
    rr = (1.0e6, -2.0e6, 6.5e6)
    t = body_torques(inertia, rr, GM, gravity_gradient=True)
    dot = sum(a * b for a, b in zip(t, rr))
    assert abs(dot) < 1e-9 * _norm(t) * _norm(rr) + 1e-20, dot

    # --- adding c*Identity to I changes nothing --------------------------
    shifted = [[inertia[i][j] + (5000.0 if i == j else 0.0)
                for j in range(3)] for i in range(3)]
    t2 = body_torques(shifted, rr, GM, gravity_gradient=True)
    for a, b in zip(t, t2):
        assert abs(a - b) < 1e-12 * (abs(a) + 1e-30), (a, b)

    # --- restoring sense about z (minor axis seeks nadir) -----------------
    # body rotated by +phi about z  =>  nadir seen in body at a = -phi
    phi = 0.05
    u = (math.cos(-phi), math.sin(-phi), 0.0)
    t = body_torques(inertia, (R * u[0], R * u[1], u[2]), GM,
                     gravity_gradient=True)
    n2 = GM / R ** 3
    assert abs(t[2] - (-3.0 * n2 * (2000.0 - 1000.0)
                       * math.sin(phi) * math.cos(phi))) < 1e-15
    assert t[2] < 0.0            # opposes +phi: x (minor axis) seeks nadir

    # --- magnetic torque is perpendicular to both m and B -----------------
    m, b = (0.31, -0.12, 0.44), (1.2e-5, 2.7e-5, -3.1e-5)
    t = magnetic_torque(m, b)
    assert abs(sum(a * c for a, c in zip(t, m))) < 1e-18
    assert abs(sum(a * c for a, c in zip(t, b))) < 1e-18

    # --- errors -----------------------------------------------------------
    for bad in (lambda: body_torques(inertia, (0.0, 0.0, 0.0), GM,
                                     gravity_gradient=True),
                lambda: body_torques(inertia, r, GM, magnetic=True)):
        try:
            bad()
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError")

    return True


if __name__ == "__main__":
    _selfcheck()
    print("torques.py self-check: all cases pass")
