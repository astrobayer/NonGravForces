function R = dayDiagnosis(options)
%DAYDIAGNOSIS One day of ACT1B against the model, hour by hour: where a day's residual comes from.
%
%   R = dayDiagnosis                                   % GRACE-D 2020-01-22 against D 01-23 and C 01-22
%   R = dayDiagnosis('Cases', ["D 2020-01-22" "D 2020-01-04"])
%
%   WHY (2026-11-01, owner: "Zerlege GRACE-D Tag fuer Tag", then this diagnosis first). The frozen
%   baseline's cross-track residual of GRACE-D over 2020-01 is 12.08 nm/s^2 against C's 2.51. Read
%   day by day from the per-day bias sigmas of tools/baseline_reference.json, ONE day carries it:
%   2020-01-22, the first day after the safe mode (19-21 have no ACT1B), about half a day of data at
%   ~85 nm/s^2 cross-track; without it the month is ~3.2 nm/s^2. tools/dragScaleMonth.m skipped the
%   day because its THR1B is not on GPS time (time_ref 'R', the receiver), so no firing was cut.
%   This tool asks WHAT on that day: thrusters, a re-commissioning transient, attitude activity,
%   flagged or missing data.
%
%   METHOD, per case (a satellite and a day) on the ACT1B/GNV1B epochs every Step_s seconds:
%   the preset's Model (drag.densityScale "none": the cross-track axis does not see the density
%   scale, and the tool must not depend on the Swarm month) against ACT1B with ONE bias per day and
%   scale 1, as dragScaleMonth does. Per GPS hour of the day it prints
%     n        ACT1B records in the hour (1 Hz: 3600 when complete) and model epochs used
%     rms      x / y / z of model minus ACT1B after the day's bias [nm/s^2]
%     mean y   the hour's mean cross-track residual (a bias that moves within the day)
%     noise y  ACT1B's own 1 Hz cross-track noise, the RMS of diff(y) over sqrt(2), model-free [nm/s^2]
%     thr      THR1B records with a firing in the hour, by time_ref (G = GPS, R = receiver), and
%              the summed on-time [s]; R records are placed at their own time tag, which is the
%              receiver's clock and close to GPS time but not corrected here
%     masked   rms y after cutting every firing as dragScaleMonth does ([t - 0.5 s, t + on + 7 s],
%              widened by the CRN kernel's 70.3 s half span), with R firings at their tag
%     rate     rms of the SCA1B body rate about its hour mean [urad/s] and the SCA1B gap count
%     flag     ACT1B records with any quality-flag bit set
%   and a day line with the same quantities. A case whose THR1B has both time references reports
%   how many records carry each.
%
%   IT LIVES IN tools/: it reads the DATAPOOL. release-checks.yml runs it ("what: daydiag").
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   Pool     1x1   string  -      DATAPOOL root; default "": NonGravForces.dataPool, DATAPOOL then READING_DATA_POOL
%   Cases    1xC   string  -      "<sat> <yyyy-mm-dd>"; default ["D 2020-01-22" "D 2020-01-23" "C 2020-01-22"]
%   Step_s   1x1   double  s      Model grid on the ACT1B epochs; default 10
%   Preset   1x1   string  -      default "standard"
%
% OUTPUT ARGUMENTS
%   R   1xC struct  .case .hours (table, one row per hour) .day (1x1 struct, the day line) .timeRefs
%
% ERRORS
%   NonGravForces:dataPool:dataPoolUnset  no 'Pool' and neither DATAPOOL nor READING_DATA_POOL set
%   NonGravForces:dayDiagnosis:dataMissing   no pool, or a case's ACT1B/GNV1B/SCA1B missing
%   gracefoReading:*, referenceSystem:*      propagated
%
% Created 2026-11-01 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = ""
    options.Cases (1,:) string = ["D 2020-01-22" "D 2020-01-23" "C 2020-01-22"]
    options.Step_s (1,1) double {mustBePositive, mustBeInteger} = 10
    options.Preset (1,1) string = "standard"
end
unit = "dayDiagnosis";
options.Pool = NonGravForces.dataPool(options.Pool, 'Option', "Pool");   % refuses when none is named
root = fileparts(fileparts(mfilename('fullpath')));
if strlength(options.Pool) == 0 || ~isfolder(options.Pool)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL at '%s' -> set DATAPOOL or pass 'Pool'", options.Pool);
end
eopFile = fullfile(options.Pool, "EOP", "finals2000A.all");
swFile = fullfile(options.Pool, "SpaceWeather", "SW-All.csv");
if ~isfile(eopFile) || ~isfile(swFile)
    eopFile = fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all');
    swFile = fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv');
end
eop = referenceSystem.readEop(eopFile);
sw = NonGravForces.io.spaceWeather(swFile);
o = NonGravForces.Options(options.Preset).set("drag.densityScale", "none");
halfCrn = (1407 - 1) / 2 * 0.1;
R = struct('case', {}, 'hours', {}, 'day', {}, 'timeRefs', {});
for c = options.Cases
    parts = split(c);
    sat = parts(1); d = parts(2);
    dayDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(d, 5), d);
    f = @(p) fullfile(dayDir, p + "_" + d + "_" + sat + "_04.txt");
    if ~all(arrayfun(@(p) isfile(f(p)), ["ACT1B" "GNV1B" "SCA1B"]))
        NonGravForces.fail(unit, "dataMissing", "%s: ACT1B, GNV1B or SCA1B missing in %s -> fetch the month into the pool", c, dayDir);
    end
    tc = tic;
    [epa, am, auxA] = NonGravForces.io.readACC1B(f("ACT1B"));
    [ep, r, v] = NonGravForces.io.readGNV1B(f("GNV1B"), 'Eop', eop);
    [epq, q] = NonGravForces.io.readSCA1B(f("SCA1B"));
    % THR1B through gracefoReading directly: NonGravForces.io.readTHR1B refuses a non-GPS time_ref,
    % which is exactly the record this diagnosis needs to see.
    tF = zeros(0, 1); onF = zeros(0, 1); refF = strings(0, 1); refs = "";
    if isfile(f("THR1B"))
        [tt, obs, raw] = gracefoReading.readL1bMulti(f("THR1B"));
        names = string(raw.obsNames(:)');
        iOn = find(startsWith(names, "on_time_"));
        on = obs(:, iOn) / 1000;                       % ms in the file
        fired = any(on > 0, 2);
        tr = repmat("G", size(tt, 1), 1);
        if isfield(raw, 'time_ref')
            x = raw.time_ref;
            if ischar(x), x = cellstr(x(:)); end
            tr = string(x(:));
            if isscalar(tr), tr = repmat(tr, size(tt, 1), 1); end
        end
        epT = referenceSystem.daySecondToEpoch(tt, "GPS");
        tF = epT.sec(fired); tF = tF(:); onF = max(on(fired, :), [], 2); refF = tr(fired);
        u = unique(tr);
        refs = strjoin(arrayfun(@(x) sprintf("%s:%d", x, nnz(tr == x)), u), " ");
    end
    % model on the shared grid
    [tf, loc] = ismember(ep.sec, epa.sec);
    tf = tf & mod(round(ep.sec), options.Step_s) == 0;
    e = referenceSystem.daySecondToEpoch([ep.day(tf) ep.sec(tf)], ep.scale);
    att = NonGravForces.geom.AttitudeProvider(epq, q);
    M = NonGravForces.Model(NonGravForces.geom.graceFO(sat), o, 'Eop', eop, 'SpaceWeather', sw);
    a = M.accel(e, r(:, tf), v(:, tf), att.at(e));
    t = ep.sec(tf); t = t(:)';
    res = a - am(:, loc(tf));
    res = res - mean(res, 2);                         % one bias per day and axis
    cut = false(1, numel(t));
    for j = 1:numel(tF)
        cut = cut | (t >= tF(j) - 0.5 - halfCrn & t <= tF(j) + onF(j) + 7 + halfCrn);
    end
    resM = res(:, ~cut); resM = resM - mean(resM, 2);
    tM = t(~cut);
    % model-free noise of ACT1B, 1 Hz
    ta = epa.sec(:)'; ya = am(2, :);
    % SCA1B body rate from successive quaternions (atan2 form, never acos)
    tq = epq.sec(:)'; qq = q; if size(qq, 1) ~= 4, qq = qq'; end
    dq = abs(sum(qq(:, 1:end-1) .* qq(:, 2:end), 1));
    ang = 2 * atan2(sqrt(max(0, 1 - dq.^2)), dq);
    dtq = diff(tq);
    rate = ang ./ dtq; tRate = tq(1:end-1);
    gapQ = tRate(dtq > 1.5 * median(dtq));
    if isfield(auxA, 'qualflgBits') && ~isempty(auxA.qualflgBits), flg = any(auxA.qualflgBits, 2)'; else, flg = false(size(ta)); end
    H = 24; rows = zeros(H, 15);
    for h = 0:H-1
        in = @(x) x >= 3600 * h & x < 3600 * (h + 1);
        k = in(t); km = in(tM); ka = in(ta); kq = in(tRate); kf = in(tF);
        yh = ya(ka); dy = diff(yh(:)'); dta = diff(ta(ka));
        noise = sqrt(mean(dy(dta < 1.5).^2) / 2);
        rh = rate(kq);
        rows(h + 1, :) = [h, nnz(ka), nnz(k), 1e9 * sqrt(mean(res(:, k).^2, 2))', 1e9 * mean(res(2, k)), 1e9 * noise, ...
            nnz(kf & refF == "G"), nnz(kf & refF ~= "G"), sum(onF(kf)), 1e9 * sqrt(mean(resM(2, km).^2)), ...
            1e6 * sqrt(mean((rh - mean(rh)).^2)), nnz(in(gapQ)), nnz(flg(ka))];
    end
    T = array2table(rows, 'VariableNames', {'hour', 'nAct1b', 'nModel', 'rmsX', 'rmsY', 'rmsZ', 'meanY', 'noiseY', ...
        'thrG', 'thrR', 'onTime_s', 'maskedRmsY', 'rateRms_urad_s', 'scaGaps', 'flagged'});
    day = struct('nAct1b', numel(ta), 'nModel', numel(t), 'rms', 1e9 * sqrt(mean(res.^2, 2))', ...
        'maskedRmsY', 1e9 * sqrt(mean(resM(2, :).^2)), 'maskedFraction', mean(cut), 'firings', numel(tF), ...
        'firingsR', nnz(refF ~= "G"), 'noiseY', 1e9 * sqrt(mean(diff(ya(:)').^2) / 2), 'flagged', nnz(flg));
    fprintf('\ndayDiagnosis %s %s: ACT1B %d records, model %d epochs at %d s, THR1B time_ref %s, %d firings (%d not GPS)  (%.0f s)\n', ...
        sat, d, day.nAct1b, day.nModel, options.Step_s, refs, day.firings, day.firingsR, toc(tc));
    fprintf('  day: rms x / y / z %.2f / %.2f / %.2f nm/s^2, y with every firing cut %.2f (%.1f %% cut), ACT1B 1 Hz noise y %.2f, flagged %d\n', ...
        day.rms, day.maskedRmsY, 100 * day.maskedFraction, day.noiseY, day.flagged);
    fprintf('  %4s %6s %5s %7s %7s %7s %7s %7s %5s %5s %7s %8s %8s %5s %5s\n', 'hour', 'nACT', 'nMod', 'rmsX', 'rmsY', 'rmsZ', ...
        'meanY', 'noiseY', 'thrG', 'thrR', 'on[s]', 'maskRmsY', 'rate[u]', 'gapQ', 'flag');
    for h = 1:H
        fprintf('  %4d %6d %5d %7.2f %7.2f %7.2f %7.2f %7.2f %5d %5d %7.2f %8.2f %8.3f %5d %5d\n', rows(h, :));
    end
    R(end + 1) = struct('case', c, 'hours', T, 'day', day, 'timeRefs', refs); %#ok<AGROW>
end
end
