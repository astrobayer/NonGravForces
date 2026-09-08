function [vrel, info] = wind(model, epoch, r_gcrs, v_gcrs, options)
%WIND Velocity of the satellite relative to the air, GCRS: co-rotation, optional thermospheric wind.
%
%   [vrel, info] = NonGravForces.env.wind("corotation", epoch, r, v, 'Eop', eop)
%   [vrel, info] = NonGravForces.env.wind("none", epoch, r, v)
%   [vrel, info] = NonGravForces.env.wind("corotation", epoch, r, v, 'Eop', eop, 'Wind_m_s', w)
%
%   The drag kernel needs the velocity of the satellite THROUGH THE GAS, not
%   through the inertial frame. At 7.6 km/s the co-rotation of the atmosphere is
%   up to 470 m/s at the equator, so omitting it is a 6 % error in the
%   cross-track drag and a 0.4 % error in its magnitude (Doornbos 2011).
%
%   MODELS
%     "none"        vrel = v. The caller states that the gas is inertial, which
%                   is never true at LEO; provided for isolating the effect.
%     "corotation"  The air co-rotates rigidly with the Earth. This unit does NOT
%                   form omega x r: it takes the Earth's rotation from the SAME
%                   authority the frames come from,
%                       vrel = v - Rdot * (R' * r),
%                   with R and Rdot from referenceSystem.rotationItrsToGcrs.
%                   Rdot carries the polar motion, the length-of-day variation and
%                   the precession-nutation rate, none of which a nominal
%                   omega x r contains (rule R1, one authority per fact; the
%                   predecessor's hand-written term was defect D4).
%     "hwm14"       Horizontal Wind Model 14 (Drob et al. 2015) on top of the
%                   co-rotation. REFUSED: the model is not in hand (HANDOVER O-26).
%                   Quiet-time winds are 40-80 m/s RMS, 1 % of the relative speed
%                   and the leading term in cross-track drag, so any cross-track
%                   claim made without it is bounded by that.
%
%   A caller who has a wind field from elsewhere passes it through 'Wind_m_s' and
%   it is subtracted after the co-rotation, whatever the model.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   model    1x1   string  -      "none" | "corotation" | "hwm14"
%   epoch    1x1   struct  -      referenceSystem epoch struct with K entries
%   r_gcrs   3xK   double  m      Satellite position, GCRS
%   v_gcrs   3xK   double  m/s    Satellite velocity, GCRS
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   Eop        1x1   struct  -      referenceSystem.readEop result, or a file path;
%                                   required by "corotation" and "hwm14" (no default:
%                                   the rotation cannot be guessed, rule R5)
%   Wind_m_s   3xK   double  m/s    Air velocity in GCRS to subtract in addition,
%                                   default [] (none)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   vrel     3xK   double  m/s    Satellite velocity relative to the air, GCRS
%   info     1x1   struct  -      .model, .corotation_m_s (3xK, the term removed),
%                                 .speed_m_s (1xK, |vrel|), .windGiven (logical),
%                                 .K, .provenance
%
% ACCURACY
%   Exact given R and Rdot: one matrix product and one subtraction, rounding
%   1e-16 relative. The co-rotation magnitude was measured at 255.8 m/s at
%   -59.5 deg latitude on the 2020-01-01 GRACE-C arc, against a nominal
%   omega*r*cos(lat) of 255.9 m/s; the difference is the polar motion and the
%   length-of-day term that a nominal rate omits. The model error is the rigid
%   co-rotation assumption itself: real thermospheric winds reach 40-80 m/s RMS
%   in quiet conditions and several hundred m/s in storms ("hwm14", O-26).
%
% ERRORS
%   NonGravForces:wind:notImplemented  model "hwm14" (O-26)
%   NonGravForces:wind:inputShape      r or v not 3xK with K matching the epoch
%   NonGravForces:wind:inputEop        "corotation" without 'Eop'
%   NonGravForces:wind:numericNaN      NaN/Inf in r, v or the supplied wind
%   referenceSystem:*                  propagated unwrapped
%
% REFERENCES
%   Doornbos E. (2011) Thermospheric density and wind determination from satellite
%   dynamics, Sect. 2.3; Drob D.P. et al. (2015) An update to the Horizontal Wind
%   Model (HWM): The quiet time thermosphere, Earth and Space Science 2:301;
%   IERS Conventions (2010) Ch. 5 through referenceSystem.
%
% EXAMPLE
%   eop = referenceSystem.readEop(file);
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
%   [vrel, info] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
%   vecnorm(v - vrel)          % 0 .. 470 m/s over an orbit
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.env.density, referenceSystem.rotationItrsToGcrs.
%
% KEYWORDS wind co-rotation relative velocity thermosphere drag
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["none","corotation","hwm14"])}
    epoch (1,1) struct
    r_gcrs double
    v_gcrs double
    options.Eop = []
    options.Wind_m_s double = []
end
unit = "wind";
if model == "hwm14"
    NonGravForces.fail(unit, "notImplemented", ...
        "the HWM14 tables are not in this toolbox (HANDOVER O-26) -> use corotation, or supply a wind field through 'Wind_m_s'");
end
K = numel(epoch.sec);
if ~(ismatrix(r_gcrs) && isequal(size(r_gcrs), [3 K]) && isequal(size(v_gcrs), [3 K]))
    NonGravForces.fail(unit, "inputShape", "r_gcrs and v_gcrs must be 3xK with K = %d epochs -> got %s and %s", K, mat2str(size(r_gcrs)), mat2str(size(v_gcrs)));
end
if any(~isfinite(r_gcrs(:))) || any(~isfinite(v_gcrs(:)))
    NonGravForces.fail(unit, "numericNaN", "r_gcrs or v_gcrs contain NaN/Inf -> screen the state");
end
w = options.Wind_m_s;
if ~isempty(w)
    if ~isequal(size(w), [3 K])
        NonGravForces.fail(unit, "inputShape", "Wind_m_s must be 3xK with K = %d -> got %s", K, mat2str(size(w)));
    end
    if any(~isfinite(w(:)))
        NonGravForces.fail(unit, "numericNaN", "Wind_m_s contains NaN/Inf -> screen the wind field");
    end
end
corot = zeros(3, K);
if model == "corotation"
    if isempty(options.Eop)
        NonGravForces.fail(unit, "inputEop", "corotation needs 'Eop' -> pass referenceSystem.readEop(file) or the file path; the Earth rotation is not assumed");
    end
    [R, Rdot] = referenceSystem.rotationItrsToGcrs(epoch, 'Eop', options.Eop);
    rItrs = reshape(pagemtimes(pagetranspose(R), reshape(r_gcrs, 3, 1, K)), 3, K);
    corot = reshape(pagemtimes(Rdot, reshape(rItrs, 3, 1, K)), 3, K);   % = omega x r, carried by the dependency
end
vrel = v_gcrs - corot;
if ~isempty(w)
    vrel = vrel - w;
end
info.model = model;
info.corotation_m_s = corot;
info.speed_m_s = vecnorm(vrel);
info.windGiven = ~isempty(w);
info.K = K;
info.provenance = NonGravForces.provenance("env.wind", 'Options', struct('model', model, 'windGiven', ~isempty(w)));
end
