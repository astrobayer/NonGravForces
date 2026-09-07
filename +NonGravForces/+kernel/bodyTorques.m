function [tau, info] = bodyTorques(macro, r_body, options)
%BODYTORQUES Gravity-gradient and magnetic torques on the body, in the body frame; both opt-in.
%
%   [tau, info] = NonGravForces.kernel.bodyTorques(macro, r_body, 'GravityGradient', true)
%   [tau, info] = NonGravForces.kernel.bodyTorques(macro, r_body, 'Magnetic', true, 'B_body', B, 'Dipole_Am2', m)
%
%   Torques that act on the body as a whole, not on its surface, and therefore
%   never enter the accelerometer (they enter the attitude, and through the
%   centre-of-mass coupling of the sensor the observable):
%
%       tau_gg = 3 GM / |r|^5  r x (I r)          (Wertz 1978 eq. 17-30; theory §7)
%       tau_m  = m x B                            (dipole m [A m^2], field B [T])
%
%   with r the geocentric position of the satellite IN THE BODY FRAME and I the
%   inertia tensor of the macro model (JPL D-56935 Table 4 for GRACE-FO). Both
%   terms are off unless requested (Options.torques.gravityGradient/.magnetic
%   are false in every preset); with neither enabled the result is zero and
%   info.enabled says so.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   macro    1x1   MacroModel     Panel model with a finite .inertia (needed by the gravity gradient)
%   r_body   3xK   double  m      Geocentric position rotated into the body frame (R' r_gcrs)
%
% NAME-VALUE ARGUMENTS
%   Name             Size      Type     Units   Description
%   ----             ----      ----     -----   -----------
%   GravityGradient  1x1       logical  -       default false
%   Magnetic         1x1       logical  -       default false
%   B_body           3xK       double   T       Magnetic field in the body frame (IGRF via env.magneticField
%                                              or MAG1B), required when Magnetic is true
%   Dipole_Am2       3x1|3xK   double   A m^2   Residual plus commanded dipole in the body frame, required
%                                              when Magnetic is true
%   GM               1x1       double   m^3/s^2 default referenceSystem.iersConstants().gmEarth
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   tau      3xK   double  N m    Sum of the enabled torques, body frame
%   info     1x1   struct  -      .gravityGradient (3xK), .magnetic (3xK), .enabled (struct of logicals),
%                                 .GM, .provenance
%
% ACCURACY
%   Closed forms; rounding 1e-16. The gravity gradient uses the point-mass
%   field: the J2 contribution to the gradient at 490 km is 3e-3 of the
%   monopole term and is neglected (documented). For GRACE-FO the gravity-
%   gradient torque is ~1e-4 N m, comparable to the aerodynamic torque.
%
% ERRORS
%   NonGravForces:bodyTorques:inputMacro     macro is not a MacroModel, or its inertia is not finite (gravity gradient requested)
%   NonGravForces:bodyTorques:inputShape     r_body not 3xK; B_body not 3xK; Dipole_Am2 not 3x1 or 3xK
%   NonGravForces:bodyTorques:inputMagnetic  Magnetic requested without B_body or Dipole_Am2
%   NonGravForces:bodyTorques:inputRange     |r| = 0 or GM <= 0
%   NonGravForces:bodyTorques:numericNaN     NaN/Inf in an input
%
% REFERENCES
%   Wertz J.R. (ed., 1978) Spacecraft Attitude Determination and Control,
%   Sect. 17.2 (gravity-gradient torque, eq. 17-30) and 17.3 (magnetic torque);
%   Wen H.Y. et al. (2019) JPL D-56935 Table 4 (inertia).
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); R = referenceSystem.rotationFromQuaternion(q);   % body -> GCRS
%   rb = squeeze(pagemtimes(pagetranspose(R), reshape(r, 3, 1, [])));
%   tau = NonGravForces.kernel.bodyTorques(mm, rb, 'GravityGradient', true);
%
% See also NonGravForces.kernel.assemble, NonGravForces.geom.graceFO, referenceSystem.iersConstants.
%
% KEYWORDS torque gravity gradient magnetic dipole body
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    r_body double
    options.GravityGradient (1,1) logical = false
    options.Magnetic (1,1) logical = false
    options.B_body double = []
    options.Dipole_Am2 double = []
    options.GM double = []
end
unit = "bodyTorques";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
end
if ~(ismatrix(r_body) && size(r_body, 1) == 3)
    NonGravForces.fail(unit, "inputShape", "r_body must be 3xK -> got %s", mat2str(size(r_body)));
end
if any(~isfinite(r_body(:)))
    NonGravForces.fail(unit, "numericNaN", "r_body contains NaN/Inf -> screen the state");
end
K = size(r_body, 2);
tauG = zeros(3, K); tauM = zeros(3, K);
if isempty(options.GM)
    GM = referenceSystem.iersConstants().gmEarth;
else
    GM = options.GM;
    if ~isscalar(GM) || ~isfinite(GM) || GM <= 0
        NonGravForces.fail(unit, "inputRange", "GM must be a positive scalar -> got %s", mat2str(GM));
    end
end
if options.GravityGradient
    I = macro.inertia;
    if any(~isfinite(I(:)))
        NonGravForces.fail(unit, "inputMacro", "the gravity gradient needs a finite inertia tensor -> build the MacroModel with 'Inertia' (graceFO carries JPL D-56935 Table 4)");
    end
    rn = vecnorm(r_body);
    if any(rn <= 0)
        NonGravForces.fail(unit, "inputRange", "|r_body| must be > 0 -> check the position");
    end
    Ir = I * r_body;                                % 3xK
    cr = [r_body(2, :) .* Ir(3, :) - r_body(3, :) .* Ir(2, :); r_body(3, :) .* Ir(1, :) - r_body(1, :) .* Ir(3, :); r_body(1, :) .* Ir(2, :) - r_body(2, :) .* Ir(1, :)];
    tauG = (3 * GM ./ rn.^5) .* cr;
end
if options.Magnetic
    B = options.B_body; m = options.Dipole_Am2;
    if isempty(B) || isempty(m)
        NonGravForces.fail(unit, "inputMagnetic", "Magnetic needs 'B_body' (3xK, T) and 'Dipole_Am2' (3x1 or 3xK) -> pass the field (env.magneticField or MAG1B) and the dipole");
    end
    if ~isequal(size(B), [3 K]) || ~(isequal(size(m), [3 1]) || isequal(size(m), [3 K]))
        NonGravForces.fail(unit, "inputShape", "B_body must be 3xK and Dipole_Am2 3x1 or 3xK -> got %s and %s", mat2str(size(B)), mat2str(size(m)));
    end
    if any(~isfinite(B(:))) || any(~isfinite(m(:)))
        NonGravForces.fail(unit, "numericNaN", "B_body or Dipole_Am2 contain NaN/Inf -> screen the inputs");
    end
    m = m + zeros(3, K);
    tauM = [m(2, :) .* B(3, :) - m(3, :) .* B(2, :); m(3, :) .* B(1, :) - m(1, :) .* B(3, :); m(1, :) .* B(2, :) - m(2, :) .* B(1, :)];
end
tau = tauG + tauM;
info.gravityGradient = tauG;
info.magnetic = tauM;
info.enabled = struct('gravityGradient', options.GravityGradient, 'magnetic', options.Magnetic);
info.GM = GM;
info.provenance = NonGravForces.provenance("kernel.bodyTorques", 'Options', struct('GravityGradient', options.GravityGradient, 'Magnetic', options.Magnetic, 'GM', GM));
end
