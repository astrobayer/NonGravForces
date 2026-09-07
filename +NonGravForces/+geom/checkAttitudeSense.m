function [ok, report] = checkAttitudeSense(q_body2gcrs, r, v, options)
%CHECKATTITUDESENSE Physics test of a BODY->GCRS quaternion series: body +Z must point to nadir, body +X along the flight line.
%
%   [ok, report] = NonGravForces.geom.checkAttitudeSense(q, r, v)
%   [ok, report] = NonGravForces.geom.checkAttitudeSense(q, r, v, 'NadirTolerance_deg', 5)
%
%   The pitfall this guards (HANDOVER C.1 and the LP-2 sense test): SCA1B is
%   ICRF->SRF and must be conjugated EXACTLY ONCE. A wrong sense is invisible
%   to every algebraic check (unit norm, orthogonality, round trips) and
%   visible only against physics: for an Earth-pointing satellite the body
%   +Z axis rotated into the GCRS with the quaternion must lie within a few
%   degrees of -r/|r|. This function measures that angle for every epoch,
%   and the same angle with the quaternion conjugated, so the report shows
%   that the alternative would have failed (rule T6: a refusal is paired
%   with the case that must succeed).
%
% INPUT ARGUMENTS
%   Name          Size  Type    Units  Description
%   ----          ----  ----    -----  -----------
%   q_body2gcrs   4xK   double  -      Quaternions to test, scalar first
%   r             3xK   double  m      GCRS position at the same epochs
%   v             3xK   double  m/s    GCRS velocity (flight-line sense)
%
% NAME-VALUE ARGUMENTS
%   Name                Size  Type    Units  Description
%   ----                ----  ----    -----  -----------
%   NadirTolerance_deg  1x1   double  deg    Largest acceptable nadir angle (default 5:
%                                            GRACE-FO pointing offsets are ~1 deg,
%                                            attitude-control deadbands smaller)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type     Units  Description
%   ----     ----  ----     -----  -----------
%   ok       1x1   logical  -      every epoch's nadir angle is within tolerance
%   report   1x1   struct   -      .nadirAngle_deg (1xK, body +Z vs -r/|r|),
%                                  .nadirAngleConjugate_deg (1xK, same with q*),
%                                  .alongTrackAngle_deg (1xK, body +X vs v/|v|),
%                                  .xSense ("alongVelocity" | "againstVelocity" |
%                                  "mixed", from the along-track angles),
%                                  .maxNadirAngle_deg, .tolerance_deg, .K,
%                                  .provenance
%
% ACCURACY
%   Exact geometry (a rotation and two angles per epoch, taken as atan2 of
%   the cross and dot products, which keeps 1e-14 deg near 0 where acos
%   would lose sqrt(eps) ~ 1e-6 deg); the verdict depends only on the
%   tolerance, which is a statement about the mission's pointing, not about
%   the arithmetic.
%
% ERRORS
%   NonGravForces:checkAttitudeSense:inputShape   q not 4xK, r or v not 3xK, or NaN/Inf
%   referenceSystem:*                             propagated unwrapped
%
% EXAMPLE
%   eop = referenceSystem.readEop(file); [ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
%   [~, q] = NonGravForces.io.readSCA1B(sca);
%   [ok, rep] = NonGravForces.geom.checkAttitudeSense(q, r, v); rep.maxNadirAngle_deg   % ~0.9 deg on 2020-01-01
%
% See also NonGravForces.io.readSCA1B, NonGravForces.geom.AttitudeProvider, referenceSystem.rotationFromQuaternion.
%
% KEYWORDS attitude sense nadir conjugation check
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    q_body2gcrs double
    r double
    v double
    options.NadirTolerance_deg (1,1) double {mustBePositive} = 5
end
unit = "checkAttitudeSense";
K = size(q_body2gcrs, 2);
if ~(ismatrix(q_body2gcrs) && size(q_body2gcrs, 1) == 4 && isequal(size(r), [3 K]) && isequal(size(v), [3 K]))
    NonGravForces.fail(unit, "inputShape", "q must be 4xK, r and v 3xK with the same K -> got q %s, r %s, v %s", mat2str(size(q_body2gcrs)), mat2str(size(r)), mat2str(size(v)));
end
if any(~isfinite([q_body2gcrs(:); r(:); v(:)]))
    NonGravForces.fail(unit, "inputShape", "q, r or v contain NaN/Inf -> screen the inputs");
end
R = referenceSystem.rotationFromQuaternion(q_body2gcrs);
Rc = referenceSystem.rotationFromQuaternion(referenceSystem.invertQuaternion(q_body2gcrs));
nadir = -r ./ vecnorm(r);
vhat = v ./ vecnorm(v);
ez = reshape(pagemtimes(R, [0; 0; 1]), 3, K);
ezc = reshape(pagemtimes(Rc, [0; 0; 1]), 3, K);
ex = reshape(pagemtimes(R, [1; 0; 0]), 3, K);
ang = @(a, b) atan2d(vecnorm(cross(a, b, 1)), sum(a .* b, 1));   % well conditioned near 0 and 180 deg (acos loses sqrt(eps))
report.nadirAngle_deg = ang(ez, nadir);
report.nadirAngleConjugate_deg = ang(ezc, nadir);
report.alongTrackAngle_deg = ang(ex, vhat);
if all(report.alongTrackAngle_deg < 90)
    report.xSense = "alongVelocity";
elseif all(report.alongTrackAngle_deg > 90)
    report.xSense = "againstVelocity";
else
    report.xSense = "mixed";
end
report.maxNadirAngle_deg = max(report.nadirAngle_deg);
report.tolerance_deg = options.NadirTolerance_deg;
report.K = K;
report.provenance = NonGravForces.provenance("geom.checkAttitudeSense", 'Options', options);
ok = report.maxNadirAngle_deg <= options.NadirTolerance_deg;
end
