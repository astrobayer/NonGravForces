function [T, info] = perRevolution(obs, ref, options)
%PERREVOLUTION Accelerometer bias and diagonal scale fitted per orbital revolution, gated per window.
%
%   T = NonGravForces.app.perRevolution(obs, ref)
%   [T, info] = NonGravForces.app.perRevolution(obs, ref, 'Period_s', 5640, 'NBoot', 50, 'Effects', effects)
%
%   The per-revolution parameters the attribution regresses on its drivers
%   (N-03c; CONCEPT 6.2 and 7; the v3 plan's WP-16 `perRev`). The arc is cut
%   into contiguous windows of 'Period_s' from its first usable epoch, and in
%   each window NonGravForces.app.fitCalibration fits a CONSTANT bias and a
%   DIAGONAL scale, with the CRN filter, gates, noise model and bootstrap it
%   always uses -- one estimator, applied per window. The cross-track scale
%   is held at its a priori unless 'Fit' says otherwise, as in fitCalibration.
%
%   A GATED WINDOW IS NaN WITH ITS REASON, NEVER A NUMBER (v3 WP-16). Two
%   gates decide:
%     "window"            fewer than 'MinWindowFraction' of the rows the
%                         window should hold are usable (a gap, the filter's
%                         NaN edges, a short last window);
%     "insufficientData"  timeSeriesPreparation.observability refused an axis
%                         (its minimum-data gate inside fitCalibration);
%     "collinear"         fitCalibration's VIF gate refused (only with a
%                         'Scale' that estimates several columns per row).
%   Every other error propagates: a gate is a statement about the data, an
%   error elsewhere is a statement about the call.
%
%   EFFECT AMPLITUDES. With 'Effects' (the per-effect a priori on the raw
%   grid, e.g. the drag, radiation and thermal arcs a Model gives through
%   Model.budget), T carries amp_<name>: the RMS of each effect per window and
%   axis, the "per-effect a priori amplitudes" CONCEPT 6.2 regresses on.
%
% INPUT ARGUMENTS
%   Name  Size  Type                        Units  Description
%   ----  ----  ----                        -----  -----------
%   obs   1x1   timeSeriesPreparation.Arc   m/s^2  Observation, as for app.fitCalibration
%   ref   1x1   timeSeriesPreparation.Arc   m/s^2  Model on the kernel's raw grid, as for app.fitCalibration
%
% NAME-VALUE ARGUMENTS
%   Name               Size  Type    Units  Description
%   ----               ----  ----    -----  -----------
%   Period_s           1x1   double  s      Window length, the orbital period; default 5640 (~490 km)
%   MinWindowFraction  1x1   double  -      Usable share of a window's rows below which it is gated;
%                                           default 0.8
%   NBoot              1x1   double  -      Bootstrap replicates per window, 0 for formal sigmas;
%                                           default 50
%   Seed               1x1   double  -      Bootstrap seed; window w uses Seed + w; default 0
%   Fit                1xN   cell    -      Further app.fitCalibration name-value pairs (Kernel,
%                                           FixCrossTrack, Noise, Robust, ...); 'BiasBasis',
%                                           'ScaleBasis', 'Span', 'NBoot' and 'Seed' are this unit's
%   Effects            1xE   struct  m/s^2  .name (string) .arc (timeSeriesPreparation.Arc, 3 axes,
%                                           any epochs); default none
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   T     RxN   table   -      One row per window: window, t0, t1, tMid (int64 GPS ns), usable
%                              (share of the window's rows), scale, sigmaScale (Rx3, NaN where
%                              held or gated), bias, sigmaBias (Rx3, m/s^2), gate (Rx3 string:
%                              "pass", "warn" -- fitted with a VIF above VifMax on one model
%                              column -- or the reason it is NaN), and amp_<name> (Rx3, m/s^2)
%                              per effect. UserData: .period_s .nBoot .sigma ("bootstrap" or
%                              "formal")
%   info  1x1   struct  -      .calibrations (1xR cell of app.Calibration, [] where gated)
%                              .messages (Rx1 string, the gate's message where gated)
%
% ACCURACY
%   Each window has app.fitCalibration's accuracy on its own rows. A window is
%   a short arc: its scale rests on one revolution of signal, so its sigma is
%   larger than a day's by roughly sqrt(revolutions), and its bootstrap block
%   length is estimated on fewer rows. The window edges are fixed in time,
%   not tied to the ascending node; a per-revolution parameter therefore
%   averages over one revolution whatever its phase.
%
% ERRORS
%   NonGravForces:perRevolution:inputArc      obs is not a timeSeriesPreparation.Arc
%   NonGravForces:perRevolution:inputFit      'Fit' sets a pair this unit owns, or is not name-value
%   NonGravForces:perRevolution:inputEffects  an effect without .name and a 3-axis .arc
%   NonGravForces:perRevolution:dataEmpty     no usable observation epoch
%   NonGravForces:fitCalibration:*            an estimator error that is not a gate
%   timeSeriesPreparation:*, lsqvce:*         as app.fitCalibration
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md 6.2 and 7; the v3 plan
%   (docs/inputs/HANDOVER_AccelCal_v3.md) WP-16.
%
% EXAMPLE
%   D = zeros(3); D(1, 1) = 0.02;
%   [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.5, 'ScaleDrift', D);
%   T = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
%   T(:, ["tMid" "scale" "gate"])
%
% See also NonGravForces.app.attributeRevolutions, NonGravForces.app.fitCalibration, NonGravForces.app.calibrate.
%
% KEYWORDS calibration per-revolution accelerometer scale bias window gate attribution
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    obs (1,1)
    ref (1,1)
    options.Period_s (1,1) double {mustBePositive, mustBeFinite} = 5640
    options.MinWindowFraction (1,1) double {mustBePositive, mustBeLessThanOrEqual(options.MinWindowFraction, 1)} = 0.8
    options.NBoot (1,1) double {mustBeInteger, mustBeNonnegative} = 50
    options.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.Fit cell = {}
    options.Effects struct = struct('name', {}, 'arc', {})
end
unit = "perRevolution";
if ~isa(obs, 'timeSeriesPreparation.Arc') || size(obs.obs, 2) ~= 3
    NonGravForces.fail(unit, "inputArc", "obs must be one timeSeriesPreparation.Arc with three acceleration columns -> got a %s", class(obs));
end
if mod(numel(options.Fit), 2) ~= 0 || any(ismember(lower(string(options.Fit(1:2:end))), ["biasbasis" "scalebasis" "span" "nboot" "seed"]))
    NonGravForces.fail(unit, "inputFit", "'Fit' must be name-value pairs and may not set BiasBasis, ScaleBasis, Span, NBoot or Seed -> those belong to the per-revolution design");
end
for e = options.Effects(:)'
    if ~(isstring(e.name) || ischar(e.name)) || ~isa(e.arc, 'timeSeriesPreparation.Arc') || size(e.arc.obs, 2) ~= 3
        NonGravForces.fail(unit, "inputEffects", "each effect needs a .name and a three-axis timeSeriesPreparation.Arc in .arc -> check 'Effects'");
    end
end

usable = obs.weights > 0 & ~obs.filled & all(isfinite(obs.obs), 2);
if ~any(usable)
    NonGravForces.fail(unit, "dataEmpty", "the observation has no usable epoch -> check the arc's weights and NaN windows");
end
P = int64(round(options.Period_s * 1e9));
tStart = obs.t(find(usable, 1));
tLast = obs.t(find(usable, 1, 'last'));
R = max(1, ceil(double(tLast - tStart + 1) / double(P)));
dt = double(median(diff(obs.t))) / 1e9;
nExpected = options.Period_s / dt;

t0 = tStart + int64(0:R - 1)' * P; t1 = t0 + P; tMid = t0 + P / 2;
use = nan(R, 1); scale = nan(R, 3); sScale = nan(R, 3); bias = nan(R, 3); sBias = nan(R, 3);
gate = strings(R, 3); messages = strings(R, 1); cals = cell(1, R);
nE = numel(options.Effects);
amp = repmat({nan(R, 3)}, 1, nE);
ax = ["x" "y" "z"];
for w = 1:R
    rows = obs.t >= t0(w) & obs.t < t1(w);
    use(w) = nnz(usable & rows) / nExpected;
    for k = 1:nE
        a = options.Effects(k).arc;
        in = a.t >= t0(w) & a.t < t1(w);
        if any(in), amp{k}(w, :) = sqrt(mean(a.obs(in, :).^2, 1, "omitnan")); end
    end
    if use(w) < options.MinWindowFraction
        gate(w, :) = "window";
        messages(w) = sprintf("%.3g of the window's rows are usable, below %g", use(w), options.MinWindowFraction);
        continue
    end
    obsW = obs.derive(obs.t(rows), obs.obs(rows, :), "NonGravForces.app.perRevolution.window", struct('window', w));
    try
        [cal, fi] = NonGravForces.app.fitCalibration(obsW, ref, options.Fit{:}, 'BiasBasis', "constant", ...
            'NBoot', options.NBoot, 'Seed', options.Seed + w);
    catch err
        switch err.identifier
            case "timeSeriesPreparation:observability:insufficientData", reason = "insufficientData";
            case "NonGravForces:fitCalibration:numericCollinear",        reason = "collinear";
            otherwise, rethrow(err);
        end
        gate(w, :) = reason; messages(w) = string(err.message);
        continue
    end
    cals{w} = cal;
    for i = 1:3
        if any(fi.names{i} == "S_" + ax(i) + ax(i))
            scale(w, i) = cal.scale(i, i);
            if options.NBoot > 0, sScale(w, i) = cal.uncertainty.sigmaScaleBootstrap(i, i);
            else, sScale(w, i) = cal.uncertainty.sigmaScaleFormal(i, i);
            end
        end
        bias(w, i) = cal.bias.coef(i, 1);
        if ~isempty(cal.bias.cov{i}), sBias(w, i) = sqrt(cal.bias.cov{i}(1, 1)); end
        gate(w, i) = cal.gates.verdict(i);
    end
end

T = table((1:R)', t0, t1, tMid, use, scale, sScale, bias, sBias, gate, ...
    'VariableNames', {'window', 't0', 't1', 'tMid', 'usable', 'scale', 'sigmaScale', 'bias', 'sigmaBias', 'gate'});
for k = 1:nE
    T.("amp_" + string(options.Effects(k).name)) = amp{k};
end
sigmaKind = "formal"; if options.NBoot > 0, sigmaKind = "bootstrap"; end
T.Properties.UserData = struct('period_s', options.Period_s, 'nBoot', options.NBoot, 'sigma', sigmaKind);
info = struct('calibrations', {cals}, 'messages', messages);
end
