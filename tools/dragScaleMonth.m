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
%   VARIANTS against the preset: msis21, constantCd 2.4 (Huckfeldt's base case) and both;
%   DTM2020 (CNES, fitted to TU Delft's accelerometer densities; academic licence) alone and
%   with Cd 2.4. Each reports its along-track bias-only RMS and k_drag.
%
%   SWARM DENSITY SCALE (O-19, 2026-10-27): when the pool holds ESA's Swarm POD-only densities
%   of the month (<Pool>/Swarm/DNSxPOD/<yyyy>/SW_OPER_DNS<x>POD_2_<yyyy-mm>.csv, fetched by
%   setupNonGravForces('FetchPool')), env.densityScale gives each day f = sum(rho_Swarm) /
%   sum(rho_model) per satellite at Swarm's own positions. Three variants: the preset's drag times
%   Swarm B's factor (~515 km, above GRACE-FO) and times A+C's (~450 km, below) -- exact, because the
%   gas force is linear in the species densities -- and the Model option drag.densityScale "swarm",
%   which interpolates log f in height to GRACE-FO's geodetic height at every epoch.
%
%   FLUX CORRECTION (O-19 step 2, 2026-10-28): with the driver variants, when 'FluxRecord' exists, the
%   variant "flux driver" -- the preset's Model on a space-weather struct whose F10.7 and F10.7A are
%   shifted per UTC day by delta = a0 + c . z, the F30 + season + square record densityDriverFit wrote
%   for its "day trimmed flux" set: an exospheric temperature correction whose height profile is MSIS's own.
%
%   COMBINED SCALE (2026-10-29): with Swarm and F30 both in the pool, the variant
%   "drag.densityScale combined" -- env.densityScale with the driver scale as 'Fallback' -- and the count of
%   epochs it took from the fallback (none when Swarm covers every day of the month).
%
%   DRIVER DENSITY SCALE (O-19, 2026-10-28): when the pool holds CLS's adjusted F30
%   (<Pool>/SpaceWeather/F30/radio_flux_adjusted_observation.txt), two Model variants with
%   drag.densityScale "driver" -- env.densityDriverScale's prediction from F30, Ap and season, at
%   Swarm's level (Offset 0) and at GRACE-FO 1's (Offset -0.041). They need no density observation
%   of the month, so they are the extrapolable correction judged against the same ACT1B.
%
%   SPACE WEATHER AND EOP: the complete files in the pool (<Pool>/SpaceWeather/SW-All.csv,
%   <Pool>/EOP/finals2000A.all, fetched with the month by setupNonGravForces('PoolMonths')) when
%   they are there, else the committed 2019-09..2020-03 excerpts, which serve 2020-01 only. A month
%   the chosen files do not cover is refused rather than evaluated with a held driver.
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
%   FluxRecord  1x1   string  -      densityDriverFit's coefficient record for the flux variant; default the
%                                    2000-08..2026-08 all-mission record in <Pool>/Swarm
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   R     1x1   struct  -      .days (table: day, nUsed, maskedFraction, thr1b, rmsX/Y/Z
%                              of the preset, kDrag of the preset, rmsX and kDrag per variant)
%                              .month (table: variant, drag mean ratio, median daily rmsX,
%                              joint kDrag) .swarm (table: day, B and A+C factors, the Model option's
%                              mean factor and GRACE-FO's mean geodetic height; empty without Swarm
%                              data) .densityScale (the env.densityScale result, [] without)
%                              .noDataDays .seconds
%
% ERRORS
%   NonGravForces:dragScaleMonth:dataMissing  no pool or no month folder, or space weather / EOP files
%                                             that do not cover the month (81 days before it included)
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
    options.FluxRecord (1,1) string = fullfile(string(getenv("READING_DATA_POOL")), "Swarm", "densityDriverFit_record_msis00_2000-08_2026-08_champ+grace+gracefo+swarm.mat")
end
unit = "dragScaleMonth";
root = fileparts(fileparts(mfilename('fullpath')));
monthDir = fullfile(options.Pool, "GRACEFO", "L1B", "RL04", extractBefore(options.Month, 5));
if strlength(options.Pool) == 0 || ~isfolder(monthDir)
    NonGravForces.fail(unit, "dataMissing", "no DATAPOOL month folder %s -> set READING_DATA_POOL or pass 'Pool'", monthDir);
end
eopFile = fullfile(options.Pool, "EOP", "finals2000A.all");
swFile = fullfile(options.Pool, "SpaceWeather", "SW-All.csv");
if ~isfile(eopFile) || ~isfile(swFile)
    eopFile = fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all');
    swFile = fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv');
end
eop = referenceSystem.readEop(eopFile);
sw = NonGravForces.io.spaceWeather(swFile);
m0 = datetime(options.Month, 'InputFormat', 'yyyy-MM');
mjd0 = juliandate(m0 - caldays(81)) - 2400000.5;   % DTM2020's 81-day mean reaches back that far
mjd1 = juliandate(dateshift(m0, 'end', 'month') + caldays(1)) - 2400000.5;
seen = sw.mjd(all(isfinite(sw.index.kp3h), 2) & isfinite(sw.index.f107Observed));   % observed rows, not predictions
swSpan = [min(seen) max(seen)];
if swSpan(1) > mjd0 || swSpan(2) < mjd1
    NonGravForces.fail(unit, "dataMissing", "%s spans MJD %d..%d, the month needs %d..%d -> fetch the month with setupNonGravForces('PoolMonths', ""%s"")", ...
        swFile, swSpan(1), swSpan(2), mjd0, mjd1, options.Month);
end
fprintf('dragScaleMonth: space weather %s, EOP %s\n', swFile, eopFile);
dd = dir(fullfile(monthDir, options.Month + "-*"));
days = string({dd([dd.isdir]).name});
if ~isempty(options.Days), days = days(ismember(str2double(extractAfter(days, 8)), options.Days)); end

o = NonGravForces.Options(options.Preset);
vName = ["preset" "msis21" "constantCd 2.4" "msis21 + Cd 2.4" "dtm2020" "dtm2020 + Cd 2.4"];
vOpt = {o, o.set("drag.density", "msis21"), o.set("drag.gsi", "constantCd").set("drag.cd", 2.4), ...
    o.set("drag.density", "msis21").set("drag.gsi", "constantCd").set("drag.cd", 2.4), ...
    o.set("drag.density", "dtm2020"), o.set("drag.density", "dtm2020").set("drag.gsi", "constantCd").set("drag.cd", 2.4)};
nM = numel(vName);                              % variants that are a Model of their own
macro = NonGravForces.geom.graceFO(options.Satellite);
models = cell(1, nM);
for v = 1:nM
    models{v} = NonGravForces.Model(macro, vOpt{v}, 'Eop', eop, 'SpaceWeather', sw);
end
% The driver-predicted density scale (env.densityDriverScale, PR #119), when the pool holds CLS's F30:
% Swarm's level (Offset 0) and GRACE-FO 1's (-0.041). It needs no density observation of the month.
f30File = fullfile(options.Pool, "SpaceWeather", "F30", "radio_flux_adjusted_observation.txt");
haveDriver = isfile(f30File);
if haveDriver
    cls = NonGravForces.io.spaceWeather(f30File);
    oD = o.set("drag.densityScale", "driver");
    for off = [0 -0.041]
        D = NonGravForces.env.densityDriverScale(sw, cls, 'Offset', off);
        vName = [vName sprintf("densityScale driver %+.3f", off)]; %#ok<AGROW>
        models{end + 1} = NonGravForces.Model(macro, oD, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', D); %#ok<AGROW>
    end
    nM = numel(vName);
    % The FLUX correction (step 2): the same drivers predict a shift delta of NRLMSISE-00's F10.7 and
    % F10.7A instead of a factor, so the height profile of the correction is MSIS's own. Its coefficients
    % are the record densityDriverFit wrote for the flux set; the Model sees a space-weather struct whose
    % query returns the shifted fluxes.
    if isfile(options.FluxRecord)
        rec = load(options.FluxRecord, 'records');
        fr = rec.records.dayTrimmedFlux;              % matlab.lang.makeValidName("day trimmed flux")
        if isfield(fr, 'coef')
            Dz = NonGravForces.env.densityDriverScale(sw, cls);
            deltaAt = @(e) localFluxShift(Dz.drivers(e), fr);
            swF = sw;
            swF.query = @(e) localShiftedQuery(sw.query(e), deltaAt(e));
            vName = [vName "flux driver"];
            models{end + 1} = NonGravForces.Model(macro, o, 'Eop', eop, 'SpaceWeather', swF);
            nM = numel(vName);
            fprintf('dragScaleMonth: flux record %s (%d rows): a0 %.4f sfu, centre %.3f sfu\n', options.FluxRecord, fr.rows, fr.a0, fr.centre_sfu);
        else
            fprintf('dragScaleMonth: the flux record in %s was refused by its fit: the flux variant is not run\n', options.FluxRecord);
        end
    else
        fprintf('dragScaleMonth: no %s: the flux variant is not run\n', options.FluxRecord);
    end
else
    fprintf('dragScaleMonth: no %s in the pool: the driver variants are not run\n', f30File);
end
% Swarm POD densities of the month, when the pool has them: the density-scale variants
swFiles = string(fullfile(options.Pool, "Swarm", "DNSxPOD", extractBefore(options.Month, 5), ...
    "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_" + options.Month + ".csv"));
haveSwarm = all(isfile(swFiles));
S = [];
if haveSwarm
    S = NonGravForces.env.densityScale(swFiles, sw, 'Model', o.drag.density);
    % drag.densityScale "combined" (PR #120): Swarm's factor on its days, the driver scale (Offset 0) on the rest
    if haveDriver
        C = NonGravForces.env.densityScale(swFiles, sw, 'Model', o.drag.density, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
        vName = [vName "drag.densityScale combined"];
        models{end + 1} = NonGravForces.Model(macro, o.set("drag.densityScale", "combined"), 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', C);
    end
    vName = [vName "drag.densityScale swarm"];
    models{end + 1} = NonGravForces.Model(macro, o.set("drag.densityScale", "swarm"), 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
    nM = numel(vName);
    vName = [vName "swarm B scale" "swarm A+C scale"];
end
nV = numel(vName);
swRows = cell(0, 5);
nDriverEpochs = 0;                              % epochs the combined scale took from the driver fallback
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
        fprintf('dragScaleMonth %s %s: skipped, a product has no records (%s)\n', options.Satellite, d, err.message);
        noData(end + 1) = d; %#ok<AGROW>
        continue
    end
    [tf, loc] = ismember(ep.sec, epa.sec);
    tf = tf & mod(round(ep.sec), options.Step_s) == 0;
    % With the Swarm scale every variant drops the GPS day's first 20 s: they belong to the previous UTC
    % day (GPS - UTC = 18 s), whose factor the month's Swarm files do not hold on the 1st.
    if haveSwarm, tf = tf & ep.sec >= 20; end
    if nnz(tf) < 100
        fprintf('dragScaleMonth %s %s: skipped, %d GNV1B epochs on the %d s grid shared with ACT1B (GNV1B %d records, ACT1B %d, first secs %s / %s)\n', ...
            options.Satellite, d, nnz(tf), options.Step_s, numel(ep.sec), numel(epa.sec), mat2str(ep.sec(1:min(3, end))', 8), mat2str(epa.sec(1:min(3, end))', 8));
        noData(end + 1) = d; %#ok<AGROW>
        continue
    end
    e = referenceSystem.daySecondToEpoch([ep.day(tf) ep.sec(tf)], ep.scale);
    att = NonGravForces.geom.AttitudeProvider(epq, q);
    qk = att.at(e);
    act = am(:, loc(tf));
    t = ep.sec(tf);
    keep = true(1, nnz(tf)); thr = isfile(f("THR1B"));
    if thr
        try
            [~, th] = NonGravForces.io.readTHR1B(f("THR1B"));
        catch err
            % GRACE-D 2020-01-22, the day after the safe mode, carries THR1B records not on GPS time
            % (time_ref 'R'); the firings cannot be placed, so the day cannot be masked, and is skipped.
            if err.identifier ~= "NonGravForces:readTHR1B:dataTimeRef", rethrow(err); end
            fprintf('dragScaleMonth %s %s: skipped, THR1B time_ref is not GPS (the firings cannot be placed)\n', options.Satellite, d);
            noData(end + 1) = d; %#ok<AGROW>
            continue
        end
        k = any(th.fired, 1);
        tFire = reshape(th.epoch.sec(k), 1, []); onT = max(th.onTime_s(:, k), [], 1);
        cut = false(numel(t), 1);
        for j = 1:numel(tFire)                 % a loop, not a K x F matrix: a day has ~10^3 firings
            cut = cut | (t(:) >= tFire(j) - 0.5 - halfCrn & t(:) <= tFire(j) + onT(j) + 7 + halfCrn);
        end
        keep = ~cut';
    end
    rmsX = zeros(1, nV); kD = zeros(1, nV); rmsYZ = [NaN NaN];
    variants = cell(1, nV);                       % {a, drag} per variant
    for vv = 1:nM
        [a, info] = models{vv}.accel(e, r(:, tf), v(:, tf), qk);
        variants{vv} = {a, info};
    end
    if haveSwarm
        jd = find(S.day == ep.day(find(tf, 1)), 1);   % the GPS day; its first 18 s were dropped above
        aP = variants{1}{1}; iP = variants{1}{2};
        fac = [S.high(jd), S.low(jd)];
        for j = 1:2                               % the preset with its density scaled by fac(j)
            iS = iP; iS.drag = fac(j) * iP.drag;
            variants{nM + j} = {aP + (fac(j) - 1) * iP.drag, iS};
        end
        iM = variants{nM}{2};
        if haveDriver, nDriverEpochs = nDriverEpochs + nnz(C.source(e) == "driver"); end
        swRows(end + 1, :) = {d, S.high(jd), S.low(jd), mean(iM.densityScale(keep)), mean(iM.lla(3, keep))}; %#ok<AGROW>
    end
    for vv = 1:nV
        a = variants{vv}{1}; info = variants{vv}{2};
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
SW = cell2table(swRows, 'VariableNames', {'day', 'fB', 'fAC', 'fModelMean', 'hGraceFO_m'});
if haveSwarm
    fprintf('  Swarm density scale rho_Swarm / rho_%s per day, B / A+C / the Model option over the kept epochs: median %.3f / %.3f / %.3f (min..max %.3f..%.3f)\n', ...
        o.drag.density, median(SW.fB), median(SW.fAC), median(SW.fModelMean), min(SW.fModelMean), max(SW.fModelMean));
    fprintf('  mean heights A / B / C over the month: %.1f / %.1f / %.1f km, GRACE-FO (geodetic) %.1f km\n', ...
        mean(S.height_m(1, :), 'omitnan') / 1e3, mean(S.height_m(2, :), 'omitnan') / 1e3, mean(S.height_m(3, :), 'omitnan') / 1e3, mean(SW.hGraceFO_m) / 1e3);
    if haveDriver
        fprintf('  combined scale: %d of %d epochs from the driver fallback (the rest from Swarm)\n', nDriverEpochs, nAll);
    end
else
    fprintf('  no Swarm densities of %s in the pool: the density-scale variants were not run\n', options.Month);
end
R = struct('days', T, 'month', M, 'swarm', SW, 'densityScale', S, 'noDataDays', noData, 'seconds', toc(t0));
end

function delta = localFluxShift(z, fr)
%LOCALFLUXSHIFT The predicted shift (sfu) per epoch from env.densityDriverScale's 8xK drivers, with the
% square recentred at the flux record's own centre (the drivers carry the factor fit's).
z(8, :) = ((z(1, :) - fr.centre_sfu) / 100) .^ 2;
delta = fr.a0 + sum(fr.coef(:) .* z, 1);
end

function q = localShiftedQuery(q, delta)
%LOCALSHIFTEDQUERY The space-weather query with F10.7 and F10.7A shifted by delta (floored at 1 sfu, as in the fit).
q.f107 = max(reshape(q.f107, 1, []) + delta, 1);
q.f107a = max(reshape(q.f107a, 1, []) + delta, 1);
end
