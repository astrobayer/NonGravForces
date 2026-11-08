function R = regressionBaseline(options)
%REGRESSIONBASELINE Calibrate the reference month 2020-01 (GRACE-C and -D) and freeze or compare it.
%
%   R = regressionBaseline                         % compare with tools/baseline_reference.json
%   R = regressionBaseline('Mode', "freeze")       % write the reference (after a deliberate change)
%   R = regressionBaseline('Days', 14:16, 'Satellites', "C", 'Mode', "trial")   % a quick trial, writes nothing
%
%   THE REGRESSION BASELINE OF N-10 (CONCEPT 6.6, ORC-10: "releases change no
%   calibrated number beyond k x bootstrap sigma"). For each satellite, every
%   day of the month that the DATAPOOL holds is calibrated JOINTLY through
%   NonGravForces.app.calibrate (N-09: one common scale, a linear bias per day,
%   a bootstrap over days) against the `standard` Model along GNV1B with SCA1B
%   attitude. The numbers that matter are frozen into the committed file
%   tools/baseline_reference.json with their sigmas and the toolbox version
%   they were last verified at; a later run compares every number with it and
%   fails a number that moved by more than K of its sigma -- the scale's
%   bootstrap sigma, a bias's formal sigma.
%
%   IT LIVES IN tools/, LIKE benchmark.m, and for the same reasons: a month in
%   the standard preset takes about 1.5 h on the owner's host (measured
%   2026-10-16: 86 s per satellite-day), which no CI job allows, and the
%   result is a property of a release, not of a commit. auditPackage.py makes
%   it unavoidable all the same: it refuses a tree whose toolbox version is not
%   the one the reference file records, so a version bump without a baseline
%   run cannot pass the audit.
%
%   A DAY WITHOUT RECORDS IS A MISSING DAY, NOT AN ERROR. GRACE-D's ACT1B,
%   GNV1B and SCA1B of 2020-01-19..21 hold num_records 0 (safe mode), and
%   gracefoReading refuses such a file (gracefoReading:readL1b:noData); that
%   one identifier is caught and the day listed in noDataDays -- every other
%   reader error still stops the run. The first month run (2026-10-16) lost
%   three hours to it, after GRACE-C had already been computed.
%
%   EACH SATELLITE IS WRITTEN AS SOON AS IT IS DONE ("freeze"), into the
%   reference that already exists, so a failure on one satellite costs only
%   that satellite, and the two can be run as separate calls.
%
%   THE DAYS run REFUSES ARE RECORDED, NOT HIDDEN. GRACE-D was in safe mode
%   around 2020-01-18..22 (timeSeriesPreparation's defect catalogue, DC-1), and
%   a missing or refused day changes what the common scale averages over; the
%   reference lists the days used, and a comparison with other days fails.
%
%   THR1B IS PASSED (2026-11-04, owner: "Yes. Go ahead with it"). Until 0.18.0 the baseline gave
%   app.calibrate no 'Thrusters', so its residuals carried every firing; tools/dayDiagnosis.m
%   (PR #131) found that GRACE-D's 2020-01-22, the day after the safe mode, re-acquires attitude on
%   thrusters (112 firings, 90.6 s on-time in hour 12) and that cutting the firings takes the day's
%   cross-track residual from ~85 to 0.75 nm/s^2. Each day's THR1B is read by gracefoReading and
%   handed to timeSeriesPreparation.run's mask stage, one Arc per day (its "thr" layer, margin from
%   the preset). A day whose THR1B file is missing is refused, never run unmasked; a THR1B file with
%   num_records 0 is a day without firings. THR1B records not on GPS time (GRACE-D 01-22, time_ref
%   R, hours 6-11) are placed at their own tag: they fall before ACT1B resumes and touch no sample.
%   The reference records the choice (meta.thrusters), and a comparison across it is refused.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default "": NonGravForces.dataPool, DATAPOOL then READING_DATA_POOL
%   Month       1x1   string  -      "yyyy-mm"; default "2020-01"
%   Satellites  1xS   string  -      Default ["C" "D"]
%   Days        1xD   double  -      Days of the month to use; default [] = every day folder present
%   Preset      1x1   string  -      Model preset; default "standard"
%   NBoot       1x1   double  -      Bootstrap replicates over days; default 100
%   K           1x1   double  -      Tolerance in sigmas; default 3
%   Mode        1x1   string  -      "compare" (default), "freeze" (write the reference) or "trial"
%                                    (compute and report only)
%   Product     1x1   string  -      The observation: "ACT1B" (default) or "ACH1B", GRACE-FO's hybrid
%                                    transplant (ACU + ACS + ACM + AC0, from the .ACX archives), which
%                                    exists for GRACE-D only (N-10b, ORC-12). Only "ACT1B" is frozen and
%                                    compared; "ACH1B" runs in "trial" mode. timeSeriesPreparation
%                                    >= 1.5 prepares ACH1B exactly as ACT1B (its X-19); until then
%                                    the arcs were labelled "ACT1B", which gave bitwise the same
%                                    result on 2020-01-05/06 (measured 2026-10-20)
%   Reference   1x1   string  -      The reference file; default tools/baseline_reference.json. Another
%                                    file serves a short trial of the comparison itself
%   Thrusters   1x1   logical -      Pass each day's THR1B to calibrate's mask stage; default true
%                                    (false reproduces the baselines frozen up to 0.18.0)
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   R     1x1   struct  -      .result (per satellite: .usedDays .refusedDays .noDataDays .scale
%                              .scaleSigma .bias .biasSigma .residualRms .seconds) .compare (table of every compared number; empty unless
%                              "compare") .pass (logical) .seconds
%
% ERRORS
%   NonGravForces:dataPool:dataPoolUnset  no 'Pool' and neither DATAPOOL nor READING_DATA_POOL set
%   NonGravForces:regressionBaseline:dataMissing  no pool, no month folder, no reference to compare, or a
%                                                 used day without THR1B while 'Thrusters' is true
%   NonGravForces:regressionBaseline:inputRange   an unknown Mode
%   NonGravForces:regressionBaseline:inputThrusters  a comparison with a reference frozen under the other
%                                                 'Thrusters' choice
%
% EXAMPLE
%   R = regressionBaseline('Pool', "E:\DATAPOOL", 'Days', 15:16, 'Satellites', "C", 'NBoot', 0, 'Mode', "trial");
%
% Created 2026-10-16 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = ""
    options.Month (1,1) string = "2020-01"
    options.Satellites (1,:) string = ["C" "D"]
    options.Days double = []
    options.Preset (1,1) string = "standard"
    options.NBoot (1,1) double {mustBeInteger, mustBeNonnegative} = 100
    options.K (1,1) double {mustBePositive} = 3
    options.Mode (1,1) string = "compare"
    options.Reference (1,1) string = ""
    options.Product (1,1) string = "ACT1B"
    options.Thrusters (1,1) logical = true
end
unit = "regressionBaseline";
options.Pool = NonGravForces.dataPool(options.Pool, 'Option', "Pool");   % refuses when none is named
if ~ismember(options.Product, ["ACT1B" "ACH1B"]) || (options.Product == "ACH1B" && (any(options.Satellites ~= "D") || options.Mode ~= "trial"))
    NonGravForces.fail(unit, "inputRange", "'Product' is ""ACT1B"", or ""ACH1B"" for GRACE-D in ""trial"" mode (the hybrid transplant exists for D only, and only ACT1B is frozen) -> got %s for %s in %s", ...
        options.Product, strjoin(options.Satellites, ","), options.Mode);
end
if ~ismember(options.Mode, ["compare" "freeze" "trial"])
    NonGravForces.fail(unit, "inputRange", "'Mode' must be ""compare"", ""freeze"" or ""trial"" -> got ""%s""", options.Mode);
end
root = fileparts(fileparts(mfilename('fullpath')));
refFile = fullfile(root, 'tools', 'baseline_reference.json');
if strlength(options.Reference) > 0, refFile = options.Reference; end
monthDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(options.Month, 5));
if strlength(options.Pool) == 0 || ~isfolder(monthDir)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL month folder %s -> set DATAPOOL or pass 'Pool', and fetch the month with orbitIntegration.fetchData(""gracefoL1b"", ...)", monthDir);
end
if options.Mode == "compare" && ~isfile(refFile)
    NonGravForces.fail(unit, "dataMissing", "%s does not exist -> run regressionBaseline('Mode', ""freeze"") once to create it", refFile);
end
eop = referenceSystem.readEop(fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
days = dir(fullfile(monthDir, options.Month + "-*"));
days = string({days([days.isdir]).name});
if ~isempty(options.Days), days = days(ismember(str2double(extractAfter(days, 8)), options.Days)); end

% The standard preset's combined density scale (0.17.0, O-19): Swarm's factor of the month from the pool, the
% F30 driver fit on any day without one. A preset that asks for it is never run unscaled.
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
t0 = tic;
result = struct();
meta = struct('toolbox', string(NonGravForces.toolboxVersion("string")), 'generated', string(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm')) + " UTC", ...
    'matlab', string(version), 'month', options.Month, 'preset', options.Preset, 'nBoot', options.NBoot, 'k', options.K, ...
    'thrusters', options.Thrusters, ...
    'note', "ORC-10: every number may move by at most k of its sigma (scale: bootstrap over days; bias: formal)");
if options.Mode == "compare"
    refMeta = jsondecode(fileread(refFile)).meta;
    refThr = isfield(refMeta, 'thrusters') && logical(refMeta.thrusters);   % absent: frozen before THR1B was passed
    if refThr ~= options.Thrusters
        NonGravForces.fail(unit, "inputThrusters", "%s was frozen with 'Thrusters' %d and this run has %d -> compare under the same choice, or re-freeze", ...
            refFile, refThr, options.Thrusters);
    end
end
for sat = options.Satellites
    ts = tic;
    [acc, orb, att, used, noData, thr] = localReadMonth(monthDir, days, sat, eop, options.Product, options.Thrusters, unit);
    fprintf('regressionBaseline %s: read %d days (%d without records), THR1B %s in %.0f s\n', sat, numel(used), numel(noData), ...
        string(repmat('passed', 1, options.Thrusters)) + string(repmat('not passed', 1, ~options.Thrusters)), toc(ts));
    thrArgs = {};
    if options.Thrusters
        if isscalar(thr), thrArgs = {'Thrusters', thr{1}}; else, thrArgs = {'Thrusters', thr}; end
    end
    M = NonGravForces.Model(NonGravForces.geom.graceFO(sat), o, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, scaleArgs{:});
    tc = tic;
    res = NonGravForces.app.calibrate(acc, M, 'Orbit', orb, 'Attitude', att, thrArgs{:}, 'Fit', {'NBoot', options.NBoot, 'BiasBasis', "linear"});
    fprintf('regressionBaseline %s: calibrate (run, model, joint fit, %d replicates) in %.0f s\n', sat, options.NBoot, toc(tc));
    if isfield(res, 'used')
        ok = res.used; cals = [res.calibration{ok}];
    else
        ok = true; cals = res.calibration;          % one day: calibrate's single-arc result (found by N-10b)
    end
    c = cals(1);
    r = struct();
    r.usedDays = used(ok);
    r.refusedDays = used(~ok);
    r.noDataDays = noData;
    r.scale = diag(c.scale).';
    r.scaleSigma = sqrt(diag(reshape(diag(c.covScale), 3, 3))).';
    r.bias = cell2mat(arrayfun(@(k) reshape(cals(k).bias.coef.', 1, []), 1:numel(cals), 'UniformOutput', false).');
    r.biasSigma = cell2mat(arrayfun(@(k) reshape(cell2mat(cellfun(@(C) sqrt(diag(C)).', cals(k).bias.cov, 'UniformOutput', false)).', 1, []), ...
        1:numel(cals), 'UniformOutput', false).');
    r.residualRms = c.uncertainty.residualRms;
    r.seconds = toc(ts);
    result.(sat) = r;
    if options.Mode == "freeze"
        % Written now, into what exists: a later satellite's failure cannot cost this one.
        if isfile(refFile), refAll = jsondecode(fileread(refFile)); else, refAll = struct('meta', meta, 'result', struct()); end
        refAll.meta = meta; refAll.result.(sat) = r;
        fid = fopen(refFile, 'w'); fwrite(fid, jsonencode(refAll, PrettyPrint = true)); fclose(fid);
        fprintf('regressionBaseline %s: frozen into %s\n', sat, refFile);
    end
end

R = struct('result', result, 'compare', table(), 'pass', true, 'seconds', toc(t0));
switch options.Mode
    case "compare"
        ref = jsondecode(fileread(refFile));
        [R.compare, R.pass] = localCompare(ref.result, result, options.K);
        if R.pass
            % Verified at this version: the audit reads meta.toolbox, so a release must pass here first.
            ref.meta.toolbox = meta.toolbox; ref.meta.verified = meta.generated;
            fid = fopen(refFile, 'w'); fwrite(fid, jsonencode(ref, PrettyPrint = true)); fclose(fid);
        end
end
log = fullfile(root, 'tools', 'baseline_last.txt');
fid = fopen(log, 'a');
fprintf(fid, "%s  NonGravForces %s  %s  mode %s  pass %d  %.0f s\n", meta.generated, meta.toolbox, options.Month, options.Mode, R.pass, R.seconds);
fclose(fid);
end

% =========================================================================================
function [acc, orb, att, used, noData, thr] = localReadMonth(monthDir, days, sat, eop, product, withThr, unit)
%LOCALREADMONTH The ACT1B arcs, GNV1B orbits, one SCA1B attitude provider and (withThr) each day's THR1B Arc
%of the days present; thr is a cell, one entry per used day, [] for a THR1B file without records.
names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
acc = timeSeriesPreparation.Arc.empty; orb = struct('epoch', {}, 'r', {}, 'v', {}); used = strings(1, 0); noData = strings(1, 0);
thr = cell(1, 0);
dayS = []; secS = []; qS = [];
for d = days
    f = @(p) fullfile(monthDir, d, p + "_" + d + "_" + sat + "_04.txt");
    if ~all(arrayfun(@(p) isfile(f(p)), [product "GNV1B" "SCA1B"])), continue, end
    try
        [epC, a] = NonGravForces.io.readACC1B(f(product));
        [ep, r, v] = NonGravForces.io.readGNV1B(f("GNV1B"), 'Eop', eop);
        [e, q] = NonGravForces.io.readSCA1B(f("SCA1B"));
    catch err
        if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end
        noData(end + 1) = d; %#ok<AGROW>
        continue
    end
    if isempty(a), noData(end + 1) = d; continue, end %#ok<AGROW>
    th = [];
    if withThr
        if ~isfile(f("THR1B"))
            NonGravForces.fail(unit, "dataMissing", "%s %s has ACT1B but no THR1B in the pool -> fetch the day, or run with 'Thrusters' false", sat, d);
        end
        try
            [tt, ot, at] = gracefoReading.readL1b(f("THR1B"), "TimeFormat", "gpsns");
            th = timeSeriesPreparation.Arc(tt, ot, at);
        catch err
            if err.identifier ~= "gracefoReading:readL1b:noData", rethrow(err); end   % no firings that day
        end
    end
    thr{end + 1} = th; %#ok<AGROW>
    acc(end + 1) = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', struct("obsNames", {names}, ...
        "variableNames", {names}, "variableFieldNames", {names}, "variableUnits", units, "fileName", product + "_" + d + "_" + sat + "_04.txt", "product", product)); %#ok<AGROW>
    orb(end + 1) = struct('epoch', ep, 'r', r, 'v', v); %#ok<AGROW>
    dayS = [dayS; e.day(:)]; secS = [secS; e.sec(:)]; qS = [qS, q]; %#ok<AGROW>
    used(end + 1) = d; %#ok<AGROW>
end
att = NonGravForces.geom.AttitudeProvider(referenceSystem.daySecondToEpoch([dayS, secS], "GPS"), qS);
end

function [T, pass] = localCompare(ref, now, k)
%LOCALCOMPARE Every frozen number against the new one, in sigmas of the frozen one.
rows = {};
for sat = string(fieldnames(ref))'
    a = ref.(sat); b = now.(sat);
    sameDays = isequal(string(a.usedDays(:)'), b.usedDays);
    rows(end + 1, :) = {sat, "usedDays", double(~sameDays), 0, sameDays}; %#ok<AGROW>
    ax = ["x" "y" "z"];
    for i = [1 3]
        z = abs(b.scale(i) - a.scale(i)) / a.scaleSigma(i);
        rows(end + 1, :) = {sat, "S_" + ax(i) + ax(i), z, k, z <= k}; %#ok<AGROW>
    end
    if sameDays
        z = abs(b.bias - a.bias) ./ a.biasSigma;
        rows(end + 1, :) = {sat, "bias (worst of " + numel(z) + ")", max(z, [], 'all'), k, max(z, [], 'all') <= k}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'satellite', 'quantity', 'sigmas', 'k', 'pass'});
pass = all(T.pass);
end
