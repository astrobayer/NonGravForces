%% Tutorial 14 - A real month from the data pool
% The committed excerpts are minutes long; a calibration that means something is days long. Whole
% days and months live in the DATA POOL, the folder shared by every astrobayer toolbox that the
% environment variable DATAPOOL names (|setupNonGravForces('FetchPool', ...)| or the toolbox
% installer fill it). This tutorial reads five days of January 2020 for both GRACE-FO satellites
% from the pool, calibrates each satellite jointly over the days - one scale, a bias per day - with
% the standard Model and the thruster firings of THR1B masked, and looks at the residual day by
% day. The days are chosen to contain a known event: GRACE-D came out of a safe mode on
% 2020-01-22 and re-acquired its attitude on thrusters.
%
% *What you will learn*
%
% * how to resolve the pool with |NonGravForces.dataPool| and where the products lie in it;
% * how to read whole Level-1B days (ACT1B, GNV1B, SCA1B, THR1B) into the Arcs, orbits and
%   attitude |app.calibrate| takes, and how to skip a day without records;
% * how to build the standard preset's density scale from the pool's Swarm month and CLS F30;
% * how to calibrate a campaign with THR1B, read the per-day residuals, and see what the firings
%   do to one day;
% * which tools in |tools/| run the full month (|regressionBaseline|, |dayDiagnosis|).
%
% *Requirements:* the data pool, with 2020-01 of GRACE-FO Level-1B (days 2020-01-22..26 of C and
% D), the Swarm DNSxPOD month 2020-01 and CLS F30; the EOP and space weather are taken from the
% pool when it has them, from the committed excerpts otherwise. Run time: about 15 to 25 minutes on
% the CI runner (the standard Model at 1 Hz over whole days dominates).
%
% Requirements: data pool
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Resolve the pool
% |NonGravForces.dataPool| returns the pool root: an explicit argument, else DATAPOOL, else the
% older READING_DATA_POOL, and it REFUSES (|NonGravForces:dataPool:dataPoolUnset|) when none is
% set - nothing is ever written into the toolbox folder as a fallback. The layout is by product:
% |GRACEFO/L1B/RL04/<yyyy>/<yyyy-mm-dd>/| for Level-1B days, |EOP/finals2000A.all|,
% |SpaceWeather/SW-All.csv| and |SpaceWeather/F30/|, |Swarm/DNSxPOD/<yyyy>/|.
pool = NonGravForces.dataPool();
root = fileparts(which('setupNonGravForces'));
monthDir = fullfile(pool, "GRACEFO", "L1B", "RL04", "2020");
eopFile = fullfile(pool, "EOP", "finals2000A.all");
swFile = fullfile(pool, "SpaceWeather", "SW-All.csv");
if ~isfile(eopFile), eopFile = fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'); end
if ~isfile(swFile), swFile = fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'); end
swarm = string(fullfile(pool, "Swarm", "DNSxPOD", "2020", "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01.csv"));
f30 = fullfile(pool, "SpaceWeather", "F30", "radio_flux_adjusted_observation.txt");
assert(isfolder(monthDir), 'the pool holds no GRACE-FO Level-1B for 2020: fetch the month (setupNonGravForces ''PoolMonths'')');
assert(all(isfile(swarm)) && isfile(f30), 'the pool lacks the Swarm month 2020-01 or CLS F30: setupNonGravForces(''FetchPool'', ...)');
fprintf('pool %s\nEOP %s\nspace weather %s\n', pool, eopFile, swFile);
days = "2020-01-" + compose("%02d", 22:26);
fprintf('days: %s\n', strjoin(days, ', '));

%% Step 2 - The auxiliary data and the density scale
% The EOP rotate GNV1B (Earth-fixed) into the GCRS; the space weather drives NRLMSISE-00. The
% standard preset multiplies the density by an observed daily scale: Swarm's density over
% NRLMSISE-00's on every day Swarm flew, a scale predicted from the F30 solar flux otherwise
% ('Fallback'). 'Model' names the density backend the scale is formed for, and it must be the
% preset's.
eop = referenceSystem.readEop(eopFile);
sw = NonGravForces.io.spaceWeather(swFile);
o = NonGravForces.Options("standard");
D = NonGravForces.env.densityDriverScale(sw, NonGravForces.io.spaceWeather(f30));
S = NonGravForces.env.densityScale(swarm, sw, 'Model', o.drag.density, 'Fallback', D);
fprintf('preset standard: density %s x scale %s, thermal %s\n', o.drag.density, o.drag.densityScale, o.thermal.model);

%% Step 3 - Read the days of both satellites
% Each day folder holds the products as |<PRODUCT>_<yyyy-mm-dd>_<sat>_04.txt|. Per satellite and
% day: ACT1B into an Arc, GNV1B into an orbit struct, SCA1B into one attitude provider over all
% days, THR1B into an Arc for the mask stage. A day whose files hold no records
% (|gracefoReading:readL1b:noData|, e.g. GRACE-D in the safe mode of 2020-01-19..21) is skipped
% and listed; every other reader error stops the run. A THR1B file without records is a day
% without firings ([]). This is the loop of |tools/regressionBaseline.m| (its |localReadMonth|),
% written out.
names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
data = struct();
for sat = ["C" "D"]
    acc = timeSeriesPreparation.Arc.empty; orb = struct('epoch', {}, 'r', {}, 'v', {});
    thr = cell(1, 0); used = strings(1, 0); skipped = strings(1, 0);
    dayS = []; secS = []; qS = [];
    for d = days
        fp = @(p) fullfile(monthDir, d, p + "_" + d + "_" + sat + "_04.txt");
        if ~all(arrayfun(@(p) isfile(fp(p)), ["ACT1B" "GNV1B" "SCA1B" "THR1B"]))
            skipped(end + 1) = d + " (files missing)"; %#ok<SAGROW>
            continue
        end
        try
            [epC, a] = NonGravForces.io.readACC1B(fp("ACT1B"));
            [ep, r, v] = NonGravForces.io.readGNV1B(fp("GNV1B"), 'Eop', eop);
            [e, q] = NonGravForces.io.readSCA1B(fp("SCA1B"));
        catch err
            if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end
            skipped(end + 1) = d + " (no records)"; %#ok<SAGROW>
            continue
        end
        th = [];
        try
            [tt, ot, at] = gracefoReading.readL1b(fp("THR1B"), "TimeFormat", "gpsns");
            th = timeSeriesPreparation.Arc(tt, ot, at);
        catch err
            if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end
        end
        aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
            "variableUnits", units, "fileName", "ACT1B_" + d + "_" + sat + "_04.txt", "product", "ACT1B");
        acc(end + 1) = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', aux); %#ok<SAGROW>
        orb(end + 1) = struct('epoch', ep, 'r', r, 'v', v); %#ok<SAGROW>
        thr{end + 1} = th; %#ok<SAGROW>
        dayS = [dayS; e.day(:)]; secS = [secS; e.sec(:)]; qS = [qS, q]; %#ok<AGROW>
        used(end + 1) = d; %#ok<SAGROW>
    end
    att = NonGravForces.geom.AttitudeProvider(referenceSystem.daySecondToEpoch([dayS, secS], "GPS"), qS);
    data.(sat) = struct('acc', acc, 'orb', orb, 'thr', {thr}, 'att', att, 'used', used, 'skipped', skipped);
    fprintf('GRACE-%s: %d days read (%s), ACT1B records per day %s\n', sat, numel(used), strjoin(used, ' '), ...
        mat2str(arrayfun(@(x) numel(x.t), acc)));
    if ~isempty(skipped), fprintf('   skipped: %s\n', strjoin(skipped, ', ')); end
end

%% Step 4 - Calibrate each satellite jointly over the days
% One |app.calibrate| call per satellite: the Arc array, the Model of THAT satellite (its own
% panel model and mass), the orbits as a struct array, one attitude provider, and 'Thrusters' as a
% cell with one THR1B Arc per day. The fit is the regression baseline's: a linear bias per day,
% one common scale; the bootstrap over days is short ('NBoot' 20; with fewer than ten days the
% estimator warns that a day bootstrap rests on few days). Each call evaluates the standard
% Model at 1 Hz over every day and is what takes the time.
out = struct();
for sat = ["C" "D"]
    dd = data.(sat);
    M = NonGravForces.Model(NonGravForces.geom.graceFO(sat), o, 'Attitude', dd.att, 'Eop', eop, ...
        'SpaceWeather', sw, 'DensityScale', S);
    t0 = tic;
    res = NonGravForces.app.calibrate(dd.acc, M, 'Orbit', dd.orb, 'Attitude', dd.att, 'Thrusters', dd.thr, ...
        'Fit', {'NBoot', 20, 'BiasBasis', "linear"});
    out.(sat) = res;
    c = res.calibration{find(res.used, 1)};
    fprintf('GRACE-%s (%.0f s): S_xx %.4f +- %.4f, S_zz %.4f +- %.4f (day bootstrap), %d of %d days used\n', sat, toc(t0), ...
        c.scale(1, 1), c.uncertainty.sigmaScaleBootstrap(1, 1), c.scale(3, 3), c.uncertainty.sigmaScaleBootstrap(3, 3), ...
        nnz(res.used), numel(res.used));
    fprintf('   residual RMS over the days [%.2f %.2f %.2f] nm/s^2\n', 1e9 * c.uncertainty.residualRms);
end

%% Step 5 - The residual day by day
% |res.fit.residual| holds one residual Arc per USED day (the fit's observation minus its
% prediction, NaN where not fitted). Its RMS per day and axis is the first thing to read: a
% residual over several days can be carried by one of them. The bias coefficients are each day's
% own (value at the day's first epoch and slope per day).
fprintf('%-4s %-11s %9s %9s %9s %12s %12s\n', 'sat', 'day', 'rms x', 'rms y', 'rms z', 'b_x(t0)', 'b_z(t0)');
rmsTab = struct();
for sat = ["C" "D"]
    res = out.(sat);
    idxU = find(res.used);
    rr = nan(numel(idxU), 3);
    for j = 1:numel(idxU)
        y = res.fit.residual{j}.obs;
        for i = 1:3
            yi = y(isfinite(y(:, i)), i);
            rr(j, i) = sqrt(mean(yi .^ 2));
        end
        b = res.calibration{idxU(j)}.bias.coef;
        fprintf('%-4s %-11s %9.2f %9.2f %9.2f %12.3e %12.3e\n', sat, data.(sat).used(idxU(j)), 1e9 * rr(j, :), b(1, 1), b(3, 1));
    end
    rmsTab.(sat) = rr;
end
fprintf('(rms in nm/s^2, bias in m/s^2)\n');
figure
ax = ["x" "y" "z"];
for i = 1:3
    subplot(1, 3, i)
    plot(1:size(rmsTab.C, 1), 1e9 * rmsTab.C(:, i), 'o-', 1:size(rmsTab.D, 1), 1e9 * rmsTab.D(:, i), 's-'), grid on
    xlabel('day of the campaign'), ylabel('residual RMS (nm/s^2)'), title("axis " + ax(i))
end
legend(["GRACE-C" "GRACE-D"], 'Location', 'best')

%% Step 6 - What the firings do to one day
% GRACE-D's 2020-01-22 is the first day after a safe mode; |tools/dayDiagnosis.m| found its ACT1B
% resuming around noon and the attitude re-acquired on thrusters, with over a hundred firings in
% one hour, and the day's cross-track residual dropping to below 1 nm/s^2 once every firing is
% cut (CLAUDE.md, PR #131). Here that day alone is calibrated twice against the same reference -
% the Model's reference of the first call, handed back as an Arc (|diagnostics.rawReference|),
% so the Model is not evaluated twice - once with its THR1B masked and once without.
iD = find(data.D.used == "2020-01-22", 1);
if isempty(iD)
    fprintf('GRACE-D 2020-01-22 is not in the pool; skipped\n');
else
    dd = data.D;
    MD = NonGravForces.Model(NonGravForces.geom.graceFO("D"), o, 'Attitude', dd.att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
    thrArgs = {};
    if ~isempty(dd.thr{iD}), thrArgs = {'Thrusters', dd.thr{iD}}; end
    r1 = NonGravForces.app.calibrate(dd.acc(iD), MD, 'Orbit', dd.orb(iD), 'Attitude', dd.att, thrArgs{:}, ...
        'Fit', {'NBoot', 0, 'BiasBasis', "linear"});
    r0 = NonGravForces.app.calibrate(dd.acc(iD), r1.diagnostics.rawReference, 'Fit', {'NBoot', 0, 'BiasBasis', "linear"});
    fprintf('GRACE-D 2020-01-22, residual RMS [x y z] in nm/s^2:\n');
    fprintf('   THR1B masked:     [%.2f %.2f %.2f], %d samples with zero weight\n', 1e9 * r1.calibration.uncertainty.residualRms, ...
        nnz(r1.arc.weights == 0));
    fprintf('   THR1B not masked: [%.2f %.2f %.2f], %d samples with zero weight\n', 1e9 * r0.calibration.uncertainty.residualRms, ...
        nnz(r0.arc.weights == 0));
    th = dd.thr{iD};
    if ~isempty(th)
        on = startsWith(th.names, "on_time_");
        fprintf('   THR1B records with a firing that day: %d\n', nnz(any(th.obs(:, on) > 0, 2)));
    end
    res0 = r0.diagnostics.fit.residual.obs(:, 2); res1 = r1.diagnostics.fit.residual.obs(:, 2);
    th0 = double(dd.acc(iD).t - dd.acc(iD).t(1)) / 3600e9;
    figure
    plot(th0, 1e9 * res0, th0, 1e9 * res1), grid on
    xlabel('hours since the first ACT1B record of the day'), ylabel('cross-track residual (nm/s^2)')
    legend(["THR1B not masked" "THR1B masked"]), title('GRACE-D, 2020-01-22')
end

%% Step 7 - The whole month: the tools
% Two tools in |tools/| take this further on the pool; they are development instruments, run on
% the CI runner by |.github/workflows/release-checks.yml|, not parts of the toolbox's surface:
%
% * |regressionBaseline| calibrates every day of 2020-01 for C and D exactly as Step 4 does
%   (THR1B passed since PR #134, 100 day-bootstrap replicates) and compares every scale and bias
%   with the frozen |tools/baseline_reference.json|; a release must leave each within three of its
%   sigmas. A quick trial: |regressionBaseline('Days', 22:24, 'Satellites', "D", 'NBoot', 0, 'Mode', "trial")|.
% * |dayDiagnosis| takes one day apart hour by hour - residual, firings, attitude rate, flags -
%   and is where the GRACE-D 2020-01-22 finding of Step 6 came from:
%   |dayDiagnosis('Cases', ["D 2020-01-22" "C 2020-01-22"])|.
%
% Neither is called here: the month takes over an hour.
fprintf('frozen baseline: %s\n', fullfile(root, 'tools', 'baseline_reference.json'));
ref = jsondecode(fileread(fullfile(root, 'tools', 'baseline_reference.json')));
for sat = ["C" "D"]
    if isfield(ref.result, sat)
        rs = ref.result.(sat);
        fprintf('GRACE-%s over %d days of 2020-01: S_xx %.4f +- %.4f, S_zz %.4f +- %.4f, residual RMS [%.2f %.2f %.2f] nm/s^2\n', ...
            sat, numel(rs.usedDays), rs.scale(1), rs.scaleSigma(1), rs.scale(3), rs.scaleSigma(3), 1e9 * rs.residualRms);
    end
end
