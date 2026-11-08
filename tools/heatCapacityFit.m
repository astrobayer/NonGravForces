function T = heatCapacityFit(options)
%HEATCAPACITYFIT Fit GRACE-FO's wall heat capacities per panel group to ACT1B through the transient thermal model.
%
%   T = heatCapacityFit                                   % 2020-01, every third day, GRACE-C and -D
%   T = heatCapacityFit('Days', 1:2, 'Satellites', "C", 'Factors', [0.5 1 2])   % a quick trial
%
%   WHY IT EXISTS. tools/alignmentMonth.m (PR #126, run 41042736099) found that the y/z lag of
%   ACT1B against the standard model is thermal inertia: with thermal.model "transient" the lags
%   fall from 8-18 s to 2-4 s, and the per-day scale moves (z 0.82 -> 1.09). The transient model's
%   heat capacities are Woeske's (2019, Table 6.1) layer-stack estimates (geom.graceFO), so they
%   now set S_zz. The owner's ruling (2026-10-30): fit them, then switch the standard preset to
%   transient with the fitted values and re-freeze the baseline once.
%
%   WHAT IT FITS. Three groups with one factor each on Woeske's value: "array" (StarboardOuter,
%   PortOuter, Zenith: solar cells + insulation + honeycomb), "kapton" (Front, Rear and the two
%   inner overhang panels: Kapton + insulation + honeycomb) and "radiator" (Nadir: teflon foil).
%   The lumped walls do not exchange heat, so a panel's temperature depends on ITS OWN capacity
%   alone: one transient Model per factor value c, with every panel scaled by c, gives every
%   panel's temperature at that c, and any combination of group factors is a SUM of per-group
%   recoils (kernel.thermalRecoil on that group's panels). Each group's recoil is pchip'd onto
%   the 10 Hz raw grid and CRN-filtered onto the ACT1B epochs (as app.calibrate does); the
%   filter is linear, so the sums are taken after it. Drag and radiation do not depend on the
%   capacities and come from the steady-state model of the same preset.
%
%   For every combination and every day, each axis is fitted as y = s*m + b0 + b1*tau (the
%   baseline's linear bias per day, one scale per day), over the day minus its first 'SpinUp_s'
%   (each day's transient model starts from the steady state at its first epoch). The score is
%   the RMS residual over all fitted samples per axis; the reported best minimises rms_y + rms_z,
%   per satellite and for both together. checkAlignment's median |lag| and the median per-day
%   scale are reported beside it, and the steadyState model is the row with NaN factors.
%
%   IT LIVES IN tools/: it is a measurement, run on the CI runner by release-checks.yml
%   ("what: heatcapacity", ~2 h for the default days) which prints the written CSV.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default "": NonGravForces.dataPool, DATAPOOL then READING_DATA_POOL
%   Month       1x1   string  -      "yyyy-mm"; default "2020-01"
%   Satellites  1xS   string  -      Default ["C" "D"]
%   Days        1xD   double  -      Days of the month; default 1:3:31 (every third day)
%   Factors     1xF   double  -      Factors on Woeske's heat capacities per group; default [0.25 0.5 1 2 4]
%   SpinUp_s    1x1   double  s      Start of each day left out of the fit and the lag; default 10800
%   Output      1x1   string  -      CSV written; default tools/heat_capacity_fit.csv
%
% OUTPUT ARGUMENTS
%   Name  Size  Type   Units  Description
%   ----  ----  ----   -----  -----------
%   T     Nx13  table  -      satellite ("C", "D", "both"), cArray, cKapton, cRadiator (NaN = steadyState),
%                             rms_x/y/z [m/s^2], scale_x/y/z (median per-day scale), lag_x/y/z (median |lag| s)
%
% ERRORS
%   NonGravForces:dataPool:dataPoolUnset  no 'Pool' and neither DATAPOOL nor READING_DATA_POOL set
%   NonGravForces:heatCapacityFit:dataMissing  no pool or month folder, or the density-scale inputs
%
% EXAMPLE
%   T = heatCapacityFit('Pool', "E:\DATAPOOL", 'Days', 1, 'Satellites', "C", 'Factors', [0.5 1 2]);
%
% Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = ""
    options.Month (1,1) string = "2020-01"
    options.Satellites (1,:) string = ["C" "D"]
    options.Days double = 1:3:31
    options.Factors (1,:) double {mustBePositive} = [0.25 0.5 1 2 4]
    options.SpinUp_s (1,1) double {mustBeNonnegative} = 10800
    options.Output (1,1) string = ""
end
unit = "heatCapacityFit";
options.Pool = NonGravForces.dataPool(options.Pool, 'Option', "Pool");   % refuses when none is named
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'tools', 'heat_capacity_fit.csv');
if strlength(options.Output) > 0, out = options.Output; end
monthDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(options.Month, 5));
if strlength(options.Pool) == 0 || ~isfolder(monthDir)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL month folder %s -> set DATAPOOL or pass 'Pool'", monthDir);
end
eop = referenceSystem.readEop(fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
days = dir(fullfile(monthDir, options.Month + "-*"));
days = string({days([days.isdir]).name});
if ~isempty(options.Days), days = days(ismember(str2double(extractAfter(days, 8)), options.Days)); end

% The standard preset's density scale, as regressionBaseline builds it.
oSS = NonGravForces.Options("standard", "thermal.model", "steadyState");
oTr = NonGravForces.Options("standard", "thermal.model", "transient");
swarm = string(fullfile(options.Pool, "Swarm", "DNSxPOD", extractBefore(options.Month, 5), "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_" + options.Month + ".csv"));
f30 = fullfile(options.Pool, "SpaceWeather", "F30", "radio_flux_adjusted_observation.txt");
if ~all(isfile(swarm)) || ~isfile(f30)
    NonGravForces.fail(unit, "dataMissing", "the standard preset takes drag.densityScale combined and the pool lacks %s -> fetch them with setupNonGravForces('FetchPool')", ...
        strjoin([swarm(~isfile(swarm)) string(f30(~isfile(f30)))], ", "));
end
D = NonGravForces.env.densityDriverScale(sw, NonGravForces.io.spaceWeather(f30));
scaleArgs = {'DensityScale', NonGravForces.env.densityScale(swarm, sw, 'Model', oSS.drag.density, 'Fallback', D)};
K = timeSeriesPreparation.crnKernel("ACC");
ids = ["timeSeriesPreparation:checkAlignment:offset" "timeSeriesPreparation:checkAlignment:grossOffset" ...
    "timeSeriesPreparation:checkAlignment:undetermined"];
was = warning('off', ids(1)); for id = ids(2:end), was(end + 1) = warning('off', id); end %#ok<AGROW>
cleanup = onCleanup(@() warning(was)); %#ok<NASGU>     % the states the caller had, restored

groups = struct('name', ["array" "kapton" "radiator"], ...
    'panels', {{["StarboardOuter" "PortOuter" "Zenith"], ["Front" "Rear" "StarboardInner" "PortInner"], "Nadir"}});
G = numel(groups.name); F = numel(options.Factors);
[iA, iK, iR] = ndgrid(1:F, 1:F, 1:F);
combos = [iA(:) iK(:) iR(:)];                         % factor indices per group
C = size(combos, 1) + 1;                              % + the steadyState row (last)
acc3 = @() struct('sse', zeros(C, 3), 'n', zeros(C, 3), 'scale', {cell(C, 3)}, 'lag', {cell(C, 3)});
S = struct();
for sat = options.Satellites
    ts = tic;
    [acc, orb, att, used] = localReadMonth(monthDir, days, sat, eop);
    macro = NonGravForces.geom.graceFO(sat);
    modelArgs = {'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, scaleArgs{:}};
    Mss = NonGravForces.Model(macro, oSS, modelArgs{:});
    Mtr = cell(1, F);
    for f = 1:F, Mtr{f} = NonGravForces.Model(localScaled(macro, options.Factors(f)), oTr, modelArgs{:}); end
    member = cellfun(@(p) ismember([macro.panels.name], p), groups.panels, 'UniformOutput', false);
    R = acc3();
    for k = 1:numel(acc)
        ss = localEval(Mss, acc(k), K, orb(k), att, macro, {});
        base = ss.drag + ss.radiation;
        Gser = cell(G, F);
        for f = 1:F
            ev = localEval(Mtr{f}, acc(k), K, orb(k), att, Mtr{f}.macro, member);
            for g = 1:G, Gser{g, f} = ev.group{g}; end
        end
        y = acc(k).obs;
        tday = double(acc(k).t - acc(k).t(1)) / 1e9;
        keep = tday >= options.SpinUp_s;
        tau = double(acc(k).t - acc(k).t(round(end / 2))) / 86400e9;
        for c = 1:C
            if c < C
                m = base + Gser{1, combos(c, 1)} + Gser{2, combos(c, 2)} + Gser{3, combos(c, 3)};
            else
                m = base + ss.thermal;
            end
            yc = NaN(size(y)); mm = NaN(size(y));
            for j = 1:3
                ok = keep & isfinite(y(:, j)) & isfinite(m(:, j));
                X = [m(:, j), ones(size(tau)), tau];
                cf = localFit(X(ok, :), y(ok, j));
                res = y(ok, j) - X(ok, :) * cf;
                R.sse(c, j) = R.sse(c, j) + sum(res .^ 2);
                R.n(c, j) = R.n(c, j) + nnz(ok);
                R.scale{c, j}(end + 1) = cf(1);
                yc(ok, j) = y(ok, j) - X(ok, 2:end) * cf(2:end);
                mm(ok, j) = m(ok, j);
            end
            rep = timeSeriesPreparation.checkAlignment(localArc(acc(k).t, yc, acc(k).names, "ACT1B_CORRECTED"), ...
                localArc(acc(k).t, mm, acc(k).names, "MODEL"));
            for j = 1:3, R.lag{c, j}(end + 1) = abs(rep.lag(j)); end
        end
        fprintf('heatCapacityFit %s %s: done (%d combinations)\n', sat, used(k), C);
    end
    S.(sat) = R;
    fprintf('heatCapacityFit %s: %d days in %.0f s\n', sat, numel(acc), toc(ts));
end

% The table: per satellite and for both together.
rows = {};
sats = string(fieldnames(S))';
labels = sats; if numel(sats) > 1, labels(end + 1) = "both"; end
for lab = labels
    if lab == "both"
        R = S.(sats(1));
        for s = sats(2:end)
            R.sse = R.sse + S.(s).sse; R.n = R.n + S.(s).n;
            R.scale = cellfun(@(a, b) [a b], R.scale, S.(s).scale, 'UniformOutput', false);
            R.lag = cellfun(@(a, b) [a b], R.lag, S.(s).lag, 'UniformOutput', false);
        end
    else
        R = S.(lab);
    end
    for c = 1:C
        if c < C, fac = options.Factors(combos(c, :)); else, fac = NaN(1, 3); end
        rms = sqrt(R.sse(c, :) ./ R.n(c, :));
        sc = cellfun(@(v) median(v, 'omitnan'), R.scale(c, :));
        lg = cellfun(@(v) median(v, 'omitnan'), R.lag(c, :));
        rows(end + 1, :) = {lab, fac(1), fac(2), fac(3), rms(1), rms(2), rms(3), sc(1), sc(2), sc(3), lg(1), lg(2), lg(3)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'satellite', 'cArray', 'cKapton', 'cRadiator', 'rms_x', 'rms_y', 'rms_z', ...
    'scale_x', 'scale_y', 'scale_z', 'lag_x', 'lag_y', 'lag_z'});
writetable(T, out);

% The summary: steadyState, Woeske's capacities (all factors 1) and the five best by rms_y + rms_z.
fprintf('\nheatCapacityFit summary (rms x/y/z nm/s^2 | median scale x/y/z | median |lag| x/y/z s)\n');
for lab = labels
    r = T(T.satellite == lab, :);
    score = r.rms_y + r.rms_z;
    ss = isnan(r.cArray);
    wo = r.cArray == 1 & r.cKapton == 1 & r.cRadiator == 1;
    score(ss) = Inf;
    [~, order] = sort(score);
    nBest = min(5, nnz(~ss));
    pick = [find(ss); find(wo); order(1:nBest)];
    tags = [repmat("steadyState", nnz(ss), 1); repmat("Woeske x1", nnz(wo), 1); "best " + (1:nBest)'];
    for i = 1:numel(pick)
        q = r(pick(i), :);
        fprintf('  %-4s %-12s c = [%5.2f %5.2f %5.2f]  rms %5.2f %5.2f %5.2f | scale %.3f %.3f %.3f | lag %5.1f %5.1f %5.1f\n', ...
            lab, tags(i), q.cArray, q.cKapton, q.cRadiator, 1e9 * q.rms_x, 1e9 * q.rms_y, 1e9 * q.rms_z, ...
            q.scale_x, q.scale_y, q.scale_z, q.lag_x, q.lag_y, q.lag_z);
    end
end
end

% =========================================================================================
function mm = localScaled(macro, c)
%LOCALSCALED The same macro model with every panel's heat capacity multiplied by c.
p = macro.panels;
for k = 1:numel(p), p(k).heatCapacity_J_m2K = c * p(k).heatCapacity_J_m2K; end
mm = NonGravForces.geom.MacroModel(p, 'Name', macro.name, 'Inertia', macro.inertia, 'Mass_kg', macro.mass_kg, ...
    'CentreOfMass', macro.centreOfMass, 'ProofMassOffset', macro.proofMassOffset, 'Reference', macro.reference);
end

function ev = localEval(M, acc, K, orb, att, macro, member)
%LOCALEVAL The Model along the orbit, onto the 10 Hz raw grid by pchip and through the CRN kernel onto
% the ACT1B epochs -- app.calibrate's construction; with 'member', the thermal recoil of each panel
% group alone (kernel.thermalRecoil on the model's own temperatures), filtered the same way.
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
[a, info] = M.accel(epUse, orb.r(:, use), orb.v(:, use), att.at(epUse));
ts = double(tOrb(use) - t0) / 1e9; tq = double(tRaw - t0) / 1e9;
names = cellstr(acc.names);
filt = @(x) timeSeriesPreparation.crnApply(localArc(tRaw, pchip(ts, x, tq).', names, "MODEL"), K, OutputEpochs = acc.t).obs;
ev = struct('total', filt(a), 'drag', filt(info.drag), 'radiation', filt(info.radiation), 'thermal', filt(info.thermal));
if ~isempty(member)
    Fk = NonGravForces.kernel.thermalRecoil(macro, info.temperature_K);     % 3 x N x K, body frame
    m = reshape(info.mass_kg, 1, []);
    ev.group = cell(1, numel(member));
    for g = 1:numel(member)
        ev.group{g} = filt(reshape(sum(Fk(:, member{g}, :), 2), 3, []) ./ m);
    end
end
end

function c = localFit(X, y)
%LOCALFIT Least squares on column-equilibrated X (the model column is ~1e-8, the others ~1).
n = vecnorm(X);
c = ((X ./ n) \ y) ./ n.';
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
