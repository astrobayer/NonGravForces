function R = dragScaleMonth(options)
%DRAGSCALEMONTH The model's along-track drag scale against ACT1B over whole days of a month (Huckfeldt et al. 2024).
%
%   R = dragScaleMonth                                   % 2020-01, GRACE-C, from READING_DATA_POOL
%   R = dragScaleMonth('Days', 1:3, 'Step_s', 30)        % a quick trial
%
%   STEP 2 OF THE DRAG INVESTIGATION (HANDOVER O-19). tNgfModel/modelDragDiagnosis found
%   on the committed 3 x 30 min that GROOPS's along-track drag is 0.48 of this toolbox's,
%   a constant scale, and that ACT1B wants k_drag ~0.2 -- weakly determined, because a
%   30 min arc with its own bias carries little drag variation. This tool repeats the
%   ACT1B side on WHOLE DAYS from the DATAPOOL, where an orbit's worth of drag variation
%   per day pins the factor.
%
%   METHOD, per day and the way Huckfeldt et al. (2024, ASR 73, Sec. 3.5, eqs. 7-10)
%   validate their model: ACT1B with the THR1B firings cut, one bias per day and axis,
%   scale 1, and the RMS of model minus calibrated ACT1B. Beside it, along-track,
%   k_drag from ACT1B - radiation - thermal = b_day + k_drag drag (radiation and thermal
%   held at 1), per day and jointly over the month with one bias per day (the per-day
%   demeaning IS the bias estimate, so both are closed form). ACT1B's own scale is part
%   of k_drag: tools/baseline_reference.json's S_xx is the calibration-side view of the
%   same product.
%
%   THRUSTER MASK: every firing of THR1B is cut over [t - 0.5 s, t + onTime + 7 s]
%   (Huckfeldt Sec. 2.1, Fig. 3), widened on each side by the CRN kernel's half span of
%   70.3 s, through which ACT1B carries the spike. A day without THR1B in the pool is
%   used unmasked and flagged.
%
%   VARIANTS, one option changed at a time against the preset: msis21, constantCd 2.4
%   (Huckfeldt's base case) and both; each reports its along-track bias-only RMS and
%   k_drag.
%
%   IT LIVES IN tools/, like regressionBaseline.m: it reads the DATAPOOL, which a commit's
%   CI does not, and takes minutes per day. .github/workflows/drag-month.yml runs it on
%   the self-hosted runner and prints the tables into the job log.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default getenv("READING_DATA_POOL")
%   Month       1x1   string  -      "yyyy-mm"; default "2020-01"
%   Satellite   1x1   string  -      "C" (default) or "D"
%   Days        1xD   double  -      Days of the month; default [] = every day folder present
%   Step_s      1x1   double  s      Evaluation grid on the ACT1B epochs; default 10
%   Preset      1x1   string  -      The reference preset; default "standard"
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   R     1x1   struct  -      .days (table: day, nUsed, maskedFraction, thr1b, rmsX/Y/Z
%                              of the preset, kDrag of the preset, rmsX and kDrag per variant)
%                              .month (table: variant, drag mean ratio, median daily rmsX,
%                              joint kDrag) .noDataDays .seconds
%
% ERRORS
%   NonGravForces:dragScaleMonth:dataMissing  no pool or no month folder
%   gracefoReading:*, referenceSystem:*       propagated, except gracefoReading:readL1b:noData,
%                                             which lists the day in noDataDays
%
% EXAMPLE
%   R = dragScaleMonth('Pool', "E:\DATAPOOL", 'Days', 1:2, 'Step_s', 30);
%
% Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = string(getenv("READING_DATA_POOL"))
    options.Month (1,1) string = "2020-01"
    options.Satellite (1,1) string {mustBeMember(options.Satellite, ["C" "D"])} = "C"
    options.Days double = []
    options.Step_s (1,1) double {mustBePositive, mustBeInteger} = 10
    options.Preset (1,1) string = "standard"
end
unit = "dragScaleMonth";
root = fileparts(fileparts(mfilename('fullpath')));
monthDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(options.Month, 5));
if strlength(options.Pool) == 0 || ~isfolder(monthDir)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL month folder %s -> set READING_DATA_POOL or pass 'Pool'", monthDir);
end
eop = referenceSystem.readEop(fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
dd = dir(fullfile(monthDir, options.Month + "-*"));
days = string({dd([dd.isdir]).name});
if ~isempty(options.Days), days = days(ismember(str2double(extractAfter(days, 8)), options.Days)); end

o = NonGravForces.Options(options.Preset);
vName = ["preset" "msis21" "constantCd 2.4" "msis21 + Cd 2.4"];
vOpt = {o, o.set("drag.density", "msis21"), o.set("drag.gsi", "constantCd").set("drag.cd", 2.4), ...
    o.set("drag.density", "msis21").set("drag.gsi", "constantCd").set("drag.cd", 2.4)};
nV = numel(vName);
macro = NonGravForces.geom.graceFO(options.Satellite);
models = cell(1, nV);
for v = 1:nV
    models{v} = NonGravForces.Model(macro, vOpt{v}, 'Eop', eop, 'SpaceWeather', sw);
end
halfCrn = (1407 - 1) / 2 * 0.1;                % s, half span of timeSeriesPreparation.crnKernel("ACC")

t0 = tic;
rows = {}; noData = strings(1, 0);
% month accumulators, per variant: sums for the joint k_drag and for the drag-mean ratio
sDY = zeros(1, nV); sDD = zeros(1, nV); sDrag = zeros(1, nV); nAll = 0;
for d = days
    td = tic;
    f = @(p) fullfile(monthDir, d, p + "_" + d + "_" + options.Satellite + "_04.txt");
    if ~all(arrayfun(@(p) isfile(f(p)), ["ACT1B" "GNV1B" "SCA1B"])), continue, end
    try
        [epa, am] = NonGravForces.io.readACC1B(f("ACT1B"));
        [ep, r, v] = NonGravForces.io.readGNV1B(f("GNV1B"), 'Eop', eop);
        [epq, q] = NonGravForces.io.readSCA1B(f("SCA1B"));
    catch err
        if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end
        noData(end + 1) = d; %#ok<AGROW>
        continue
    end
    [tf, loc] = ismember(ep.sec, epa.sec);
    tf = tf & mod(round(ep.sec), options.Step_s) == 0;
    if nnz(tf) < 100, noData(end + 1) = d; continue, end %#ok<AGROW>
    e = referenceSystem.daySecondToEpoch([ep.day(tf) ep.sec(tf)], ep.scale);
    att = NonGravForces.geom.AttitudeProvider(epq, q);
    qk = att.at(e);
    act = am(:, loc(tf));
    t = ep.sec(tf);
    keep = true(1, nnz(tf)); thr = isfile(f("THR1B"));
    if thr
        [~, th] = NonGravForces.io.readTHR1B(f("THR1B"));
        k = any(th.fired, 1);
        tFire = reshape(th.epoch.sec(k), 1, []); onT = max(th.onTime_s(:, k), [], 1);
        cut = false(numel(t), 1);
        for j = 1:numel(tFire)                 % a loop, not a K x F matrix: a day has ~10^3 firings
            cut = cut | (t(:) >= tFire(j) - 0.5 - halfCrn & t(:) <= tFire(j) + onT(j) + 7 + halfCrn);
        end
        keep = ~cut';
    end
    rmsX = zeros(1, nV); kD = zeros(1, nV); rmsYZ = [NaN NaN];
    for vv = 1:nV
        [a, info] = models{vv}.accel(e, r(:, tf), v(:, tf), qk);
        res = a(:, keep) - act(:, keep);
        res = res - mean(res, 2);
        rmsX(vv) = sqrt(mean(res(1, :).^2));
        if vv == 1, rmsYZ = sqrt(mean(res(2:3, :).^2, 2))'; end
        y = act(1, keep) - info.radiation(1, keep) - info.thermal(1, keep);
        dr = info.drag(1, keep);
        yt = y - mean(y); dt = dr - mean(dr);
        kD(vv) = (dt * yt') / (dt * dt');
        sDY(vv) = sDY(vv) + dt * yt'; sDD(vv) = sDD(vv) + dt * dt'; sDrag(vv) = sDrag(vv) + sum(dr);
    end
    nAll = nAll + nnz(keep);
    rows(end + 1, :) = {d, nnz(keep), 1 - mean(keep), thr, rmsX(1), rmsYZ(1), rmsYZ(2), kD(1), rmsX, kD}; %#ok<AGROW>
    fprintf('dragScaleMonth %s %s: %5d epochs (%4.1f %% masked%s), rms %.2f / %.2f / %.2f nm/s^2, k_drag %.3f; variants rms_x %s, k %s  (%.0f s)\n', ...
        options.Satellite, d, nnz(keep), 100 * (1 - mean(keep)), string(repmat(', no THR1B', 1, ~thr)), 1e9 * [rmsX(1) rmsYZ], kD(1), ...
        sprintf('%.2f ', 1e9 * rmsX), sprintf('%.3f ', kD), toc(td));
end
if isempty(rows)
    NonGravForces.fail(unit, "dataMissing", "no day of %s in %s had ACT1B, GNV1B and SCA1B with records -> check the pool", options.Month, monthDir);
end
T = cell2table(rows, 'VariableNames', {'day', 'nUsed', 'maskedFraction', 'thr1b', 'rmsX', 'rmsY', 'rmsZ', 'kDrag', 'rmsXVariants', 'kDragVariants'});
rmsV = T.rmsXVariants;                          % days x variants
if iscell(rmsV), rmsV = cell2mat(rmsV); end
M = table(vName', (sDrag / sDrag(1))', median(rmsV, 1)', (sDY ./ sDD)', ...
    'VariableNames', {'variant', 'dragMeanRatio', 'medianDailyRmsX', 'jointKDrag'});
fprintf('\ndragScaleMonth %s %s: %d days used, %d without records, %d epochs at %d s\n', options.Satellite, options.Month, height(T), numel(noData), nAll, options.Step_s);
fprintf('  %-16s %10s %16s %12s\n', 'variant', 'drag ratio', 'median rms_x', 'joint k_drag');
for vv = 1:nV
    fprintf('  %-16s %10.3f %13.2f nm %12.3f\n', vName(vv), M.dragMeanRatio(vv), 1e9 * M.medianDailyRmsX(vv), M.jointKDrag(vv));
end
fprintf('  preset median daily rms y / z: %.2f / %.2f nm/s^2; k_drag per day min / median / max: %.3f / %.3f / %.3f\n', ...
    1e9 * median(T.rmsY), 1e9 * median(T.rmsZ), min(T.kDrag), median(T.kDrag), max(T.kDrag));
R = struct('days', T, 'month', M, 'noDataDays', noData, 'seconds', toc(t0));
end
