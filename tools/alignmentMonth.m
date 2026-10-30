function T = alignmentMonth(options)
%ALIGNMENTMONTH Per day and axis, the lag between ACT1B and the model, raw and with the bias removed, per effect.
%
%   T = alignmentMonth                                   % 2020-01, GRACE-C and -D, from the DATAPOOL
%   T = alignmentMonth('Days', 1:3, 'Satellites', "C")    % a quick trial
%
%   WHY IT EXISTS. The 0.17.0 baseline run (release-checks, run 40911119628) raised 80
%   timeSeriesPreparation:checkAlignment warnings inside app.calibrate, one per day and axis
%   above 1 s. The lags were NOT common to the axes, as a time-tag error would make them:
%   along-track -12..+6 s, cross-track and radial growing over the month from 3-7 s to the
%   60 s edge of the search, on both satellites. Two explanations fit that pattern:
%     (1) the estimator: checkAlignment removes only the MEAN of each column before the cross-
%         correlation, and ACT1B carries its whole bias (y ~ -1.13e-5 m/s^2) with its drift, which
%         on a weak axis can pull the correlation peak towards the edge;
%     (2) the model: on y and z radiation and thermal recoil dominate, and a phase error there
%         (thermal re-emission, shadow transitions) reads as a lag.
%   This tool separates them. For each day it evaluates the Model along GNV1B exactly as
%   app.calibrate does (pchip onto the 10 Hz raw grid, CRN-filtered with the ACC kernel onto
%   the ACT1B epochs), and measures checkAlignment's lag and peak correlation per axis for:
%     raw         ACT1B against the filtered model -- what the calibration run reported;
%     linear      ACT1B minus b0 + b1*tau, fitted per axis jointly with a model scale
%                 (y = s*m + b0 + b1*tau, the baseline's "linear" bias per day), against the model;
%     quadratic   the same with b2*tau^2 as well;
%     drag        the linear-bias-corrected ACT1B minus s * (radiation + thermal), against the
%                 filtered drag alone;
%     nondrag     the linear-bias-corrected ACT1B minus s * drag, against radiation + thermal.
%   If the lags vanish in "linear", the warnings are the estimator's and the timing is right; if they
%   remain and sit in "nondrag", the model's radiation/thermal phase is what lags. tau is the day's
%   time in days from its middle epoch. No series is shifted, nothing is calibrated for keeps.
%
%   IT LIVES IN tools/ like regressionBaseline.m: a month in the standard preset is ~1 h on the CI
%   runner, a diagnostic, not a gate. release-checks.yml runs it ("what: alignment") and prints the
%   written CSV into the log.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default getenv("READING_DATA_POOL")
%   Month       1x1   string  -      "yyyy-mm"; default "2020-01"
%   Satellites  1xS   string  -      Default ["C" "D"]
%   Days        1xD   double  -      Days of the month to use; default [] = every day folder present
%   Preset      1x1   string  -      Model preset; default "standard"
%   MaxLag      1x1   double  s      checkAlignment's search bound; default 60 (timeSeriesPreparation's Config)
%   Output      1x1   string  -      CSV written; default tools/alignment_month.csv
%
% OUTPUT ARGUMENTS
%   Name  Size  Type   Units  Description
%   ----  ----  ----   -----  -----------
%   T     Nx9   table  -      satellite, day, variant, axis, lag_s, peak, verdict, scale (s of the joint
%                             fit, NaN for raw), rmsRatio (rms of the corrected obs over rms of the raw
%                             obs after its mean; 1 for raw)
%
% ERRORS
%   NonGravForces:alignmentMonth:dataMissing  no pool or month folder, or the density-scale inputs
%
% EXAMPLE
%   T = alignmentMonth('Pool', "E:\DATAPOOL", 'Days', 1:2, 'Satellites', "C");
%
% Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = string(getenv("READING_DATA_POOL"))
    options.Month (1,1) string = "2020-01"
    options.Satellites (1,:) string = ["C" "D"]
    options.Days double = []
    options.Preset (1,1) string = "standard"
    options.MaxLag (1,1) double {mustBePositive} = 60
    options.Output (1,1) string = ""
end
unit = "alignmentMonth";
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'tools', 'alignment_month.csv');
if strlength(options.Output) > 0, out = options.Output; end
monthDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(options.Month, 5));
if strlength(options.Pool) == 0 || ~isfolder(monthDir)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL month folder %s -> set READING_DATA_POOL or pass 'Pool'", monthDir);
end
eop = referenceSystem.readEop(fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
days = dir(fullfile(monthDir, options.Month + "-*"));
days = string({days([days.isdir]).name});
if ~isempty(options.Days), days = days(ismember(str2double(extractAfter(days, 8)), options.Days)); end

% The preset's density scale, as regressionBaseline builds it.
o = NonGravForces.Options(options.Preset);
scaleArgs = {};
if o.drag.densityScale ~= "none"
    swarm = string(fullfile(options.Pool, "Swarm", "DNSxPOD", extractBefore(options.Month, 5), "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_" + options.Month + ".csv"));
    f30 = fullfile(options.Pool, "SpaceWeather", "F30", "radio_flux_adjusted_observation.txt");
    if ~all(isfile(swarm)) || ~isfile(f30)
        NonGravForces.fail(unit, "dataMissing", "preset %s takes drag.densityScale %s and the pool lacks %s -> fetch them with setupNonGravForces('FetchPool')", ...
            options.Preset, o.drag.densityScale, strjoin([swarm(~isfile(swarm)) string(f30(~isfile(f30)))], ", "));
    end
    D = NonGravForces.env.densityDriverScale(sw, NonGravForces.io.spaceWeather(f30));
    switch o.drag.densityScale
        case "driver",   scaleArgs = {'DensityScale', D};
        case "swarm",    scaleArgs = {'DensityScale', NonGravForces.env.densityScale(swarm, sw, 'Model', o.drag.density)};
        case "combined", scaleArgs = {'DensityScale', NonGravForces.env.densityScale(swarm, sw, 'Model', o.drag.density, 'Fallback', D)};
    end
end
K = timeSeriesPreparation.crnKernel("ACC");
ids = ["timeSeriesPreparation:checkAlignment:offset" "timeSeriesPreparation:checkAlignment:grossOffset" ...
    "timeSeriesPreparation:checkAlignment:undetermined"];
was = warning('off', ids(1)); for id = ids(2:end), was(end + 1) = warning('off', id); end %#ok<AGROW>
cleanup = onCleanup(@() warning(was)); %#ok<NASGU>     % the states the caller had, restored

rows = {};
variants = ["raw" "linear" "quadratic" "drag" "nondrag"];
ax = ["x" "y" "z"];
for sat = options.Satellites
    ts = tic;
    [acc, orb, att, used] = localReadMonth(monthDir, days, sat, eop);
    M = NonGravForces.Model(NonGravForces.geom.graceFO(sat), o, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, scaleArgs{:});
    for k = 1:numel(acc)
        F = localFilteredModel(M, acc(k), K, orb(k), att);
        y = acc(k).obs;
        tau = double(acc(k).t - acc(k).t(round(end / 2))) / 86400e9;
        reps = cell(1, numel(variants)); s = NaN(numel(variants), 3); ratio = ones(numel(variants), 3);
        reps{1} = timeSeriesPreparation.checkAlignment(acc(k), F.total, MaxLag = options.MaxLag);
        for v = 2:numel(variants)
            yc = NaN(size(y)); refObs = NaN(size(y));
            for j = 1:3
                m = F.total.obs(:, j);
                ok = isfinite(y(:, j)) & isfinite(m);
                switch variants(v)
                    case "quadratic", X = [m, ones(size(m)), tau, tau .^ 2];
                    otherwise,        X = [m, ones(size(m)), tau];
                end
                c = localFit(X(ok, :), y(ok, j));
                yb = y(:, j) - X(:, 2:end) * c(2:end);        % the bias curve removed
                s(v, j) = c(1);
                switch variants(v)
                    case "drag"
                        yc(:, j) = yb - c(1) * (F.radiation.obs(:, j) + F.thermal.obs(:, j));
                        refObs(:, j) = F.drag.obs(:, j);
                    case "nondrag"
                        yc(:, j) = yb - c(1) * F.drag.obs(:, j);
                        refObs(:, j) = F.radiation.obs(:, j) + F.thermal.obs(:, j);
                    otherwise
                        yc(:, j) = yb;
                        refObs(:, j) = m;
                end
                y0 = y(ok, j) - mean(y(ok, j)); yv = yc(ok, j) - mean(yc(ok, j));
                ratio(v, j) = sqrt(mean(yv .^ 2) / mean(y0 .^ 2));
            end
            reps{v} = timeSeriesPreparation.checkAlignment(localArc(acc(k).t, yc, acc(k).names, "ACT1B_CORRECTED"), ...
                localArc(acc(k).t, refObs, acc(k).names, "MODEL_" + upper(variants(v))), MaxLag = options.MaxLag);
        end
        for v = 1:numel(variants)
            for j = 1:3
                rows(end + 1, :) = {sat, used(k), variants(v), ax(j), reps{v}.lag(j), reps{v}.peak(j), reps{v}.verdict(j), s(v, j), ratio(v, j)}; %#ok<AGROW>
            end
        end
        fprintf('alignmentMonth %s %s: raw lag %s s | linear %s s | nondrag %s s\n', sat, used(k), ...
            localLags(reps{1}), localLags(reps{2}), localLags(reps{5}));
    end
    fprintf('alignmentMonth %s: %d days in %.0f s\n', sat, numel(acc), toc(ts));
end
T = cell2table(rows, 'VariableNames', {'satellite', 'day', 'variant', 'axis', 'lag_s', 'peak', 'verdict', 'scale', 'rmsRatio'});
writetable(T, out);

% The summary: per satellite, variant and axis, the median |lag|, the largest |lag| and the median peak.
fprintf('\nalignmentMonth summary (median |lag| s / max |lag| s / median peak / days aligned within 1 s)\n');
for sat = options.Satellites
    for v = variants
        parts = strings(1, 3);
        for j = 1:3
            r = T(T.satellite == sat & T.variant == v & T.axis == ax(j), :);
            L = abs(r.lag_s);
            parts(j) = sprintf('%s %6.2f / %5.1f / %.3f / %2d of %2d', ax(j), median(L, 'omitnan'), max(L), median(r.peak, 'omitnan'), ...
                nnz(L <= 1), height(r));
        end
        fprintf('  %s %-9s %s\n', sat, v, strjoin(parts, "  |  "));
    end
end
end

% =========================================================================================
function F = localFilteredModel(M, acc, K, orb, att)
%LOCALFILTEREDMODEL The Model and its effects along the orbit, onto the 10 Hz raw grid by pchip and
% through the CRN kernel onto the ACT1B epochs -- the construction of app.calibrate's localModelReference.
ep = orb.epoch;
if string(ep.scale) ~= "GPS", ep = referenceSystem.convertTimeScale(ep, "GPS"); end
tOrb = referenceSystem.epochTo(ep, "gpsNanoseconds");
step = int64(round(1e9 / K.fs));
half = int64(K.Nh) * step;
t0 = acc.t(1);
lo = max(acc.t(1) - half, tOrb(1)); hi = min(acc.t(end) + half, tOrb(end));
k0 = idivide(lo - t0, step, "ceil"); k1 = idivide(hi - t0, step, "floor");
tRaw = t0 + (k0:k1).' * step;
use = find(tOrb >= tRaw(1), 1) - 1:find(tOrb <= tRaw(end), 1, 'last') + 1;
use = use(use >= 1 & use <= numel(tOrb));
epUse = referenceSystem.daySecondToEpoch([ep.day(use), ep.sec(use)], ep.scale);
[a, ai] = M.accel(epUse, orb.r(:, use), orb.v(:, use), att.at(epUse));
ts = double(tOrb(use) - t0) / 1e9; tq = double(tRaw - t0) / 1e9;
names = cellstr(acc.names);
F = struct();
F.total = timeSeriesPreparation.crnApply(localArc(tRaw, pchip(ts, a, tq).', names, "MODEL"), K, OutputEpochs = acc.t);
for e = ["drag" "radiation" "thermal"]
    F.(e) = timeSeriesPreparation.crnApply(localArc(tRaw, pchip(ts, ai.(e), tq).', names, "MODEL_" + upper(e)), K, OutputEpochs = acc.t);
end
end

function c = localFit(X, y)
%LOCALFIT Least squares on column-equilibrated X (the model column is ~1e-8, the others ~1).
n = vecnorm(X);
c = ((X ./ n) \ y) ./ n.';
end

function s = localLags(rep)
s = strjoin(compose("%.2f", rep.lag), " ");
end

function arc = localArc(t, X, names, product)
%LOCALARC An arc with the observation's column names, as app.calibrate builds its reference.
names = cellstr(names);
u = struct();
for k = 1:numel(names), u.(names{k}) = 'm/s2'; end
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", u, "fileName", product, "product", product);
arc = timeSeriesPreparation.Arc(t, X, aux);
end

function [acc, orb, att, used] = localReadMonth(monthDir, days, sat, eop)
%LOCALREADMONTH The ACT1B arcs, GNV1B orbits and one SCA1B attitude provider of the days present
% (regressionBaseline's reader; a day without records is skipped).
names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
acc = timeSeriesPreparation.Arc.empty; orb = struct('epoch', {}, 'r', {}, 'v', {}); used = strings(1, 0);
dayS = []; secS = []; qS = [];
for d = days
    f = @(p) fullfile(monthDir, d, p + "_" + d + "_" + sat + "_04.txt");
    if ~all(arrayfun(@(p) isfile(f(p)), ["ACT1B" "GNV1B" "SCA1B"])), continue, end
    try
        [epC, a] = NonGravForces.io.readACC1B(f("ACT1B"));
        [ep, r, v] = NonGravForces.io.readGNV1B(f("GNV1B"), 'Eop', eop);
        [e, q] = NonGravForces.io.readSCA1B(f("SCA1B"));
    catch err
        if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end
        continue
    end
    if isempty(a), continue, end
    acc(end + 1) = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', struct("obsNames", {names}, ...
        "variableNames", {names}, "variableFieldNames", {names}, "variableUnits", units, "fileName", "ACT1B_" + d + "_" + sat + "_04.txt", "product", "ACT1B")); %#ok<AGROW>
    orb(end + 1) = struct('epoch', ep, 'r', r, 'v', v); %#ok<AGROW>
    dayS = [dayS; e.day(:)]; secS = [secS; e.sec(:)]; qS = [qS, q]; %#ok<AGROW>
    used(end + 1) = d; %#ok<AGROW>
end
att = NonGravForces.geom.AttitudeProvider(referenceSystem.daySecondToEpoch([dayS, secS], "GPS"), qS);
end
