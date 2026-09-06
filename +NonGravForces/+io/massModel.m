function mass = massModel(satellite, options)
%MASSMODEL Spacecraft mass model: handbook constant, or MAS1B-driven variation about it.
%
%   mass = NonGravForces.io.massModel("C")
%   mass = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epMas, 'GasMass_kg', tank.total_kg)
%   m_kg = mass.at(epoch)                        % 1xK
%
%   Returns the object NonGravForces.Model consumes for the mass. Two modes:
%
%   "constant" (default): the handbook value of NonGravForces.const.graceFO
%   (601.214 kg GRACE-C, 601.209 kg GRACE-D), or 'Anchor_kg' if given, at every
%   epoch. Standard practice for accelerometer work (HANDOVER C.3).
%
%   "mas1b": m(t) = anchor + gas(t) - gas(t_1), where gas(t) is the cold-gas
%   tank mass from NonGravForces.io.readMAS1B, linearly interpolated between
%   its (hourly) samples. The anchor is taken to hold at the FIRST sample of
%   the series; the variation within the series is then applied exactly. The
%   evaluator holds the nearest sample value for up to two median sample
%   intervals beyond either end (the first MAS1B record of a day falls at
%   3599.99 s, so the first hour would otherwise be unreachable) and refuses
%   any epoch further out.
%
%   MAS1B IS NOT THE SPACECRAFT MASS (HANDOVER C.1): this function is the
%   only place tank masses become a mass, and it refuses to return anything
%   at or below 500 kg.
%
% INPUT ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   satellite   1x1   string  -      "C" or "D".
%
% NAME-VALUE ARGUMENTS
%   Name         Size  Type    Units  Description
%   ----         ----  ----    -----  -----------
%   Mode         1x1   string  -      "constant" (default) or "mas1b".
%   Epoch        1x1   struct  -      referenceSystem epoch struct of the MAS1B
%                                     samples (mas1b only; strictly increasing).
%   GasMass_kg   1xN   double  kg     Tank total at those epochs (mas1b only;
%                                     finite; readMAS1B's tank.total_kg).
%   Anchor_kg    1x1   double  kg     Mass at the anchor (constant: the value;
%                                     mas1b: the mass at the first sample).
%                                     Default: the handbook constant.
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   mass     1x1   struct  -      .kind ("constant"|"mas1b"), .satellite,
%                                 .value_kg (the anchor), .source ("JPL D-56935
%                                 handbook" | "caller" | "MAS1B about <source>"),
%                                 .span (1x2 decimal MJD in .scale, [-Inf Inf]
%                                 for constant), .scale (epoch scale of the
%                                 samples, "" for constant), .at (function
%                                 handle: epoch struct -> 1xK kg), .provenance.
%
% ACCURACY
%   constant: the handbook value; GRACE-FO consumed about 30 kg of cold gas
%   over the mission, so a fixed value is wrong by at most 5e-2 relative
%   across years and by less than 1e-5 within a day. mas1b: the anchor
%   assumption (handbook value valid at the first sample) leaves the same
%   systematic offset but removes the daily variation; linear interpolation
%   between hourly samples of a slowly varying mass is exact to the sample
%   noise (mass_tnk_err, a few grams).
%
% ERRORS
%   NonGravForces:massModel:inputMas1b    Mode "mas1b" without Epoch/GasMass_kg,
%                                         mismatched lengths, non-finite gas
%                                         mass, or non-increasing epochs
%   NonGravForces:massModel:inputAnchor   Anchor_kg <= 500 kg
%   NonGravForces:massModel:inputEpoch    (from .at) query epoch further than two
%                                         sample intervals beyond the series, or
%                                         not an epoch struct
%   NonGravForces:massModel:numericRange  (from .at) a resulting mass <= 500 kg
%   referenceSystem:*                     propagated unwrapped (scale conversion)
%
% REFERENCES
%   JPL D-56935 (masses); HANDOVER C.1, C.3.
%
% EXAMPLE
%   m = NonGravForces.io.massModel("C");
%   m.at(referenceSystem.mjdToEpoch(58849.5, "GPS"))   % 601.214
%
% See also NonGravForces.io.readMAS1B, NonGravForces.const, NonGravForces.Model.
%
% KEYWORDS mass model handbook MAS1B propellant
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    satellite (1,1) string {mustBeMember(satellite, ["C","D"])}
    options.Mode (1,1) string {mustBeMember(options.Mode, ["constant","mas1b"])} = "constant"
    options.Epoch = []
    options.GasMass_kg (1,:) double = []
    options.Anchor_kg (1,:) double = []
end
unit = "massModel";
c = NonGravForces.const;
if isempty(options.Anchor_kg)
    anchor = c.graceFO.("mass" + satellite);
    source = "JPL D-56935 handbook";
else
    anchor = options.Anchor_kg(1);
    source = "caller";
end
if ~(isscalar(anchor) && isfinite(anchor) && anchor > 500)
    NonGravForces.fail(unit, "inputAnchor", "Anchor_kg = %g must be a finite mass above 500 kg -> a GRACE-FO satellite weighs about 601 kg; tank masses (~30 kg) are not a spacecraft mass", anchor);
end
mass.kind = options.Mode;
mass.satellite = satellite;
mass.value_kg = anchor;
switch options.Mode
    case "constant"
        mass.source = source;
        mass.span = [-Inf Inf];
        mass.scale = "";
        mass.at = @(epoch) localConstant(epoch, anchor, unit);
    case "mas1b"
        ep = options.Epoch; gas = options.GasMass_kg;
        if ~(isstruct(ep) && all(isfield(ep, ["day","sec","scale"]))) || isempty(gas)
            NonGravForces.fail(unit, "inputMas1b", "Mode mas1b needs 'Epoch' (referenceSystem epoch struct) and 'GasMass_kg' -> pass readMAS1B's epoch and tank.total_kg");
        end
        tGas = double(ep.day(:)) + double(ep.sec(:)) / 86400;
        if numel(tGas) ~= numel(gas) || ~all(isfinite(gas)) || numel(tGas) < 2 || any(diff(tGas) <= 0)
            NonGravForces.fail(unit, "inputMas1b", "GasMass_kg (%d) and Epoch (%d) must have equal length >= 2, finite masses and strictly increasing epochs -> screen the MAS1B series first", ...
                numel(gas), numel(tGas));
        end
        step = median(diff(tGas));
        mass.source = "MAS1B about " + source;
        mass.span = [tGas(1) tGas(end)];
        mass.scale = string(ep.scale);
        mass.at = @(epoch) localMas1b(epoch, anchor, tGas, gas(:), 2 * step, string(ep.scale), unit);
end
optRec = options; optRec.Epoch = mass.span;   % hash the span, not the series
mass.provenance = NonGravForces.provenance("io.massModel", 'Options', optRec);
end

function m = localConstant(epoch, anchor, unit)
K = localCount(epoch, unit);
m = repmat(anchor, 1, K);
end

function m = localMas1b(epoch, anchor, tGas, gas, hold, scale, unit)
localCount(epoch, unit);
if string(epoch.scale) ~= scale
    epoch = referenceSystem.convertTimeScale(epoch, scale);
end
t = double(epoch.day(:)) + double(epoch.sec(:)) / 86400;
lo = tGas(1) - hold; hi = tGas(end) + hold;
out = t < lo | t > hi;
if any(out)
    k = find(out, 1);
    NonGravForces.fail(unit, "inputEpoch", "epoch %d (MJD %.5f %s) is outside the MAS1B series [%.5f, %.5f] plus the hold of %.3f d -> read the neighbouring MAS1B day as well", ...
        k, t(k), scale, tGas(1), tGas(end), hold);
end
tc = min(max(t, tGas(1)), tGas(end));          % hold the end values inside the margin
m = anchor + (interp1(tGas, gas, tc, 'linear') - gas(1));
m = m(:)';
if any(m <= 500)
    NonGravForces.fail(unit, "numericRange", "a mass of %.3f kg resulted -> the anchor or the gas series is wrong; MAS1B totals are ~30 kg, the spacecraft ~601 kg", min(m));
end
end

function K = localCount(epoch, unit)
if ~(isstruct(epoch) && all(isfield(epoch, ["day","sec","scale"])))
    NonGravForces.fail(unit, "inputEpoch", "at() takes a referenceSystem epoch struct (.day .sec .scale) -> build one with daySecondToEpoch or mjdToEpoch");
end
K = numel(epoch.day);
end
