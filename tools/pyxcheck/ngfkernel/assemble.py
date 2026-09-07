"""Rigid-body assembly of per-panel forces into a net force and net torque.

Independent oracle for the ``assemble`` kernel of NonGravForces.  Written from
first principles (rigid-body statics), standard library only.

DERIVATION
==========

Setup
-----
A satellite body is modelled as a finite set of flat plates, indexed by k.  The
surface interaction (drag, SRP, ERP, thermal re-radiation, ...) has already been
reduced, for each plate, to a single resultant force vector

    F_k                                                            [N]

applied at the plate's centre of pressure.  In the panel model the centre of
pressure of a *flat, uniformly illuminated / uniformly impinged* plate coincides
with its area centroid, because every surface element of the plate receives the
same momentum flux per unit area and carries the same outward normal n_k; the
force per unit area is therefore constant over the plate and its resultant acts
through the first area moment, i.e. the centroid

    c_k = (1/A_k) * integral over plate of r dA                    [m]

That is why the plate dict's "centroid" is the point of application.  (The
optical and gas-surface coefficients enter F_k, not this step; nothing in this
module needs "area", "normal", "vis" or "ir".)

Net force
---------
Force is a free vector.  Newton's second law for the whole body follows from
summing the momentum balance of each surface element; all internal constraint
forces of the rigid body cancel pairwise (Newton's third law), so only the
external surface loads survive:

    F = sum_k F_k                                                  (1)

The sum is over vectors, component by component.  F does not depend on where the
forces are applied, and in particular it does not depend on r_cm.

Net torque about the centre of mass
-----------------------------------
The moment (torque) of a single force F_k applied at the point c_k, taken about
a reference point r_cm, is by definition the cross product of the moment arm
with the force,

    tau_k = (c_k - r_cm) x F_k                                     (2)

with the moment arm drawn FROM the reference point TO the point of application.
This is the standard right-handed convention: the resulting tau is the axial
vector whose direction is given by the right-hand rule about which the force
tends to rotate the body, and whose magnitude is |r_perp| |F_k|, the force times
the perpendicular lever arm.  With that convention the rotational equation of
motion of the rigid body reads

    dL/dt = tau ,      L = I omega            (angular momentum about the c.o.m.)

so that a positive tau_z increases L_z; the opposite sign convention would
require dL/dt = -tau, which is not used anywhere in mechanics.  Equation (2) is
the definition adopted here.

Moments about a common point add as vectors (they are the antisymmetric part of
the first moment of the load distribution, which is additive), hence

    tau = sum_k (c_k - r_cm) x F_k                                 (3)

Equations (1) and (3) are the whole content of this module.

Sign convention, stated once
----------------------------
    moment arm = (point of application) - (reference point),
    tau = arm x F,  right-handed cross product
        a x b = (a2 b3 - a3 b2,  a3 b1 - a1 b3,  a1 b2 - a2 b1).

Properties that follow (useful as metamorphic checks)
-----------------------------------------------------
(a) Reference-point transfer.  For two reference points p and q,

        tau(q) = sum_k (c_k - q) x F_k
               = sum_k (c_k - p) x F_k + sum_k (p - q) x F_k
               = tau(p) + (p - q) x F                              (4)

    so shifting the reference point changes the torque by (p - q) x F only.

(b) A pure couple (F = 0) has a torque independent of the reference point; by
    (4) tau(q) = tau(p) for all p, q.

(c) Linearity.  assemble is linear in the list {F_k}: scaling all F_k by lambda
    scales both outputs by lambda, and the outputs of a concatenated plate list
    are the sums of the outputs of the parts.  This is why the "one kernel"
    design is legitimate: force and torque are two reductions over the same
    per-panel array, and no separate torque model is needed.

(d) Translation invariance.  Translating the whole geometry, centroids and
    reference point alike, by any vector d leaves both outputs unchanged, since
    (c_k + d) - (r_cm + d) = c_k - r_cm.

(e) A force whose line of action passes through r_cm (c_k - r_cm parallel to
    F_k) contributes zero torque.

Method: closed form.  No quadrature is used anywhere.

Numerical note
--------------
The sums are accumulated component-wise in the given plate order with plain
floating-point addition, i.e. exactly the order a scalar loop would use.  No
BLAS-style reassociation, no pairwise or compensated summation, so a bitwise
comparison against a column-wise MATLAB sum is meaningful.

Units: F_k in N, centroids and r_cm in m, hence tau in N m.

Created 2026-09-07 by Claude Opus 5 (Anthropic).
"""

__all__ = ["assemble", "cross", "torque_about", "shift_reference"]


def cross(a, b):
    """Right-handed cross product of two 3-vectors."""
    a1, a2, a3 = a
    b1, b2, b3 = b
    return (a2 * b3 - a3 * b2,
            a3 * b1 - a1 * b3,
            a1 * b2 - a2 * b1)


def assemble(fk, plates, rcm):
    """Net force and net torque about ``rcm`` of a set of per-plate forces.

    Parameters
    ----------
    fk : list of 3-vectors
        Force on each plate, in N, in the body frame.
    plates : list of plate dicts
        Only the "centroid" entry is used: the point of application of the
        corresponding force, in m, in the same frame.
    rcm : 3-vector
        Reference point (the centre of mass), in m, same frame.

    Returns
    -------
    (force, torque) : tuple of two 3-tuples
        force  = sum_k F_k                       [N]
        torque = sum_k (c_k - rcm) x F_k         [N m]

    Raises
    ------
    ValueError
        If ``fk`` and ``plates`` differ in length, or a vector is not 3-long.
    """
    if len(fk) != len(plates):
        raise ValueError("assemble: len(fk)=%d does not match len(plates)=%d"
                         % (len(fk), len(plates)))
    if len(rcm) != 3:
        raise ValueError("assemble: rcm must be a 3-vector")

    x0, y0, z0 = (float(rcm[0]), float(rcm[1]), float(rcm[2]))

    fx = fy = fz = 0.0
    tx = ty = tz = 0.0

    for k, (f, plate) in enumerate(zip(fk, plates)):
        if len(f) != 3:
            raise ValueError("assemble: fk[%d] must be a 3-vector" % k)
        c = plate["centroid"]
        if len(c) != 3:
            raise ValueError("assemble: plates[%d]['centroid'] must be a "
                             "3-vector" % k)

        f1, f2, f3 = float(f[0]), float(f[1]), float(f[2])

        # moment arm: from the reference point to the point of application
        r1 = float(c[0]) - x0
        r2 = float(c[1]) - y0
        r3 = float(c[2]) - z0

        fx += f1
        fy += f2
        fz += f3

        tx += r2 * f3 - r3 * f2
        ty += r3 * f1 - r1 * f3
        tz += r1 * f2 - r2 * f1

    return (fx, fy, fz), (tx, ty, tz)


def torque_about(fk, plates, point):
    """Net torque about an arbitrary reference point (convenience wrapper)."""
    return assemble(fk, plates, point)[1]


def shift_reference(torque, force, p, q):
    """Transfer a torque from reference point ``p`` to reference point ``q``.

    Implements equation (4):  tau(q) = tau(p) + (p - q) x F.
    """
    d = (p[0] - q[0], p[1] - q[1], p[2] - q[2])
    m = cross(d, force)
    return (torque[0] + m[0], torque[1] + m[1], torque[2] + m[2])


# ---------------------------------------------------------------------------
# Self-test: cases whose answer is predicted by hand in the report.
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    def plate(centroid):
        # only "centroid" matters here; the rest is filled for interface parity
        return {"area": 1.0, "normal": (0.0, 0.0, 1.0), "centroid": centroid,
                "vis": (0.4, 0.3, 0.3), "ir": (0.8, 0.1, 0.1)}

    def show(name, got, want):
        ok = all(abs(g - w) <= 1e-15 * max(1.0, abs(w))
                 for g, w in zip(got, want))
        print("%-34s %-42s %s" % (name, tuple(got), "OK" if ok else
                                  "FAIL expected %s" % (tuple(want),)))
        return ok

    allok = True

    # Case 1 -- single unit force, unit lever arm.
    # F = (0,0,-1) N applied at c = (1,0,0) m, rcm = origin.
    # tau = (1,0,0) x (0,0,-1) = (0*(-1)-0*0, 0*0-1*(-1), 1*0-0*0) = (0,1,0)
    F, T = assemble([(0.0, 0.0, -1.0)], [plate((1.0, 0.0, 0.0))],
                    (0.0, 0.0, 0.0))
    allok &= show("1 force", F, (0.0, 0.0, -1.0))
    allok &= show("1 torque", T, (0.0, 1.0, 0.0))

    # Case 2 -- pure couple: F1 = (0,1,0) at (1,0,0), F2 = (0,-1,0) at (-1,0,0).
    # F = 0.  tau1 = (1,0,0)x(0,1,0) = (0,0,1);  tau2 = (-1,0,0)x(0,-1,0) = (0,0,1)
    # tau = (0,0,2) N m, and by property (b) it is the same about any point.
    fk = [(0.0, 1.0, 0.0), (0.0, -1.0, 0.0)]
    pl = [plate((1.0, 0.0, 0.0)), plate((-1.0, 0.0, 0.0))]
    F, T = assemble(fk, pl, (0.0, 0.0, 0.0))
    allok &= show("2 force", F, (0.0, 0.0, 0.0))
    allok &= show("2 torque", T, (0.0, 0.0, 2.0))
    F2, T2 = assemble(fk, pl, (3.7, -2.1, 0.9))
    allok &= show("2 torque (shifted rcm)", T2, (0.0, 0.0, 2.0))

    # Case 3 -- line of action through rcm gives zero torque (property (e)).
    F, T = assemble([(0.0, 0.0, -1.0)], [plate((1.0, 0.0, 0.0))],
                    (1.0, 0.0, 0.0))
    allok &= show("3 torque", T, (0.0, 0.0, 0.0))

    # Case 4 -- reference-point transfer, equation (4), on a generic load.
    fk = [(1.0, -2.0, 0.5), (-0.25, 0.75, 3.0), (2.0, 0.0, -1.5)]
    pl = [plate((0.6, -0.3, 1.1)), plate((-0.9, 0.4, -0.2)),
          plate((0.05, 0.05, 0.05))]
    p = (0.0, 0.0, 0.0)
    q = (0.021, -0.013, 0.007)
    Fp, Tp = assemble(fk, pl, p)
    Fq, Tq = assemble(fk, pl, q)
    allok &= show("4 torque about q", Tq, shift_reference(Tp, Fp, p, q))
    allok &= show("4 force invariant", Fq, Fp)

    print("ALL OK" if allok else "FAILURES PRESENT")
