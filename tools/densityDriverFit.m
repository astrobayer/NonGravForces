function R = densityDriverFit(options)
%DENSITYDRIVERFIT A driver model of the observed density scale: Swarm 2014-now, cross-validated over years.
%
%   R = densityDriverFit                                  % Swarm 2014-01..2026-08 from READING_DATA_POOL
%   R = densityDriverFit('First', "2019-01", 'Last', "2020-12")
%   R = densityDriverFit('First', "2000-08", 'Missions', ["swarm" "champ" "grace" "gracefo"])
%
%   STEP 1 OF THE EXTRAPOLABLE DENSITY CORRECTION (HANDOVER O-19, owner 2026-10-27). The
%   Swarm factor of env.densityScale exists only where Swarm measured. A correction that
%   extrapolates in time, backwards before 2014 and forwards into forecasts, must depend on
%   DRIVERS that exist for every date (F10.7 since 1947, F30 since 1951, ap since 1932) and
%   on time only through season, local time and a slow trend. This tool measures how well
%   such a regression predicts the observed factor on data it was not fitted to.
%
%   DATA: for every Swarm satellite, UTC day and orbit node (ascending, descending), the
%   observed factor f = sum(rho_Swarm) / sum(rho_model) over every 'Stride'-th valid record
%   (the ratio of sums of env.densityScale), with the pass's mean geodetic height and its
%   local solar time at the equator (circular mean of |lat| < 30 deg). Two nodes per day
%   give two local times 12 h apart; the planes drift through all local times in ~4.5
%   months (A/C) and ~4.8 months (B).
%   MISSIONS (2026-10-28, step 3): 'Missions' adds TU Delft's ACCELEROMETER densities of CHAMP
%   (2000-2010, ~300-450 km), GRACE-A (2002-2017, ~350-500 km) and GRACE-FO 1 (2018-, ~490 km),
%   ESA TOLEOS v02 from <Pool>/TOLEOS (setupNonGravForces('DensityMonths')), each its own
%   satellite in the rows. They bring solar cycle 23 -- Fbar up to ~250 sfu in 2000-2002 -- which
%   2014-2019 never reached, so the solar terms are fitted where 2024 lies instead of extrapolated.
%
%   MODEL: log f = c0 + c . (z - mean z), the drivers z
%     Fbar    81-day mean of the solar proxy ENDING the previous day (no future data)
%     dF      the previous day's proxy minus Fbar
%     ap      daily Ap of the day
%     annual and semiannual harmonics of the day of year
%     diurnal harmonic of the pass's local time
%     height  (h - 470 km) / 100 km
%     trend   years since 2020-01-01
%   with the proxy F10.7 (CelesTrak) or F30 (CLS, adjusted), nested from 3 to 11 drivers.
%   NONLINEAR SOLAR TERMS (2026-10-28): log Fbar instead of Fbar, a centred square ((Fbar - c)/100)^2 with c
%   each proxy's mean over the data set (it was 120 sfu for both until PR #118, which made F30's collinear),
%   and the cycle's PHASE, Fbar now minus Fbar a year earlier (/100; > 0 rising, < 0 declining) -- the
%   linear cycle-24 relation failed forward into cycle 25, and f is not monotonic in Fbar at high
%   activity (0.95 in 2015, 0.77 in 2024 at higher flux), which a saturation or a hysteresis explains.
%   Each fit is timeSeriesPreparation.regressDrivers (CONCEPT 3.2: no estimator here); a
%   prediction is the algebra coef(1) + (z - mean z_train) coef(2:end).
%
%   VALIDATION, the point of the tool:
%     - blocked cross-validation with one fold per calendar year (timeSeriesPreparation.
%       crossval, a buffer of 27 days between training and test blocks): the RMS of the
%       held-out log f, i.e. the relative error of a predicted factor;
%     - forward: fitted before 'Split', predicting after it (solar cycle 24 decline ->
%       cycle 25 rise), and backward: fitted after, predicting before;
%     - the constant model (the mean factor) as the baseline every driver must beat.
%   Then the full fit of each proxy's largest model with NBoot block-bootstrap sigmas.
%
%   FLUX CORRECTION (2026-10-28, step 2 of O-19): a factor on the density carries no height
%   dependence of its own, so missions at different heights read different factors (CHAMP +0.036,
%   GRACE-FO 1 -0.041 in log f). The 'flux' data sets instead correct the model THROUGH ITS OWN
%   PHYSICS: per row the shift delta (sfu) of NRLMSISE-00's F10.7 and F10.7A inputs -- an exospheric
%   temperature correction, the height profile then MSIS's own -- that makes the model's summed
%   density equal the observed one, found on a grid of shifts ('FluxGrid', default -60..+30 sfu,
%   pchip in log rho; a row outside the grid is dropped and counted). The drivers and validation are
%   those above, with delta as the response; every loss and residual of a flux set is multiplied by
%   the row's own d log rho / d delta, so it reads in log f and compares one to one with the
%   factor sets. The F30 "+season F^2" fit of every set is kept as a coefficient record in
%   <Pool>/Swarm/densityDriverFit_record_<model>_<span>[_<missions>].mat (tools/dragScaleMonth.m reads
%   the flux one of 2000-08..2026-08, all missions).
%
%   IT LIVES IN tools/, like dragScaleMonth.m: it reads the DATAPOOL and runs for minutes.
%   .github/workflows/density-fit.yml runs it on the self-hosted runner.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default getenv("READING_DATA_POOL")
%   First, Last 1x1   string  -      "yyyy-mm" span of months; default 2014-01..2026-08
%   Missions    1xM   string  -      "swarm" (A, B, C; satellites 1-3), "champ", "grace" (GRACE-A),
%                                    "gracefo" (GRACE-FO 1), numbered on in this order; default "swarm".
%                                    A month outside a mission's span is skipped, not missing
%   Model       1x1   string  -      env.density backend the factor is for; default "msis00"
%   Stride      1x1   double  -      every Stride-th valid record (30 s Swarm, 10 s TOLEOS); default 10
%   MinRecords  1x1   double  -      fewer records in a pass drop it; default 40
%   Split       1x1   string  -      "yyyy-mm-dd" of the forward/backward split; default 2020-01-01
%   NBoot       1x1   double  -      bootstrap replicates of the full fits; default 200
%   MaxVif      1x1   double  -      regressDrivers' VIF gate; default 10
%   Sets        1xS   string  -      data sets validated: "pass", "pass trimmed", "day", "day trimmed",
%                                    "pass trimmed flux", "day trimmed flux" (default "day trimmed" and
%                                    "day trimmed flux"; the raw ones are noise-dominated)
%   FluxGrid    1xG   double  sfu    shifts of F10.7 and F10.7A for the flux sets; default -60:15:30
%   Cache       1x1   logical -      keep the passes in <Pool>/Swarm/densityDriverFit_<model>_<span>_s<Stride>_n<MinRecords>[_<missions>].mat
%                                    and reuse them (default true); false recomputes and does not write
%
% OUTPUT ARGUMENTS
%   R  1x1 struct: .rows (table: mjd, satellite, ascending, logF, height_m, lt_h, n, sumObs, sumModel),
%      .satellites (the labels of the satellite numbers),
%      .drivers (names), .sets.pass, .passTrimmed, .day, .dayTrimmed (each: .cv table of model, cvRms, cvMae,
%      forwardRms, forwardMae, backwardRms, backwardMae; .cvSe; .fit, the regressDrivers result of each proxy's
%      largest model; .satelliteResidual (per proxy: rows, median, mean, RMS of log f per satellite); .nRows;
%      .squareCentre [F10.7 F30] in sfu), .missingFiles, .emptyFiles, .seconds
%
%   DATA SETS (2026-10-28, after the first run found std(log f) 0.65 per pass while the yearly medians
%   moved only between 0.5 and 1.07): per pass and per day (both nodes, the ratio of the summed sums),
%   each also TRIMMED of rows more than log 2 from their satellite's 27-day running median. The losses
%   are reported as RMS and as median absolute error, so that a few noise-dominated passes cannot decide
%   between models.
%
% ERRORS
%   NonGravForces:densityDriverFit:dataMissing  SW-All or the CLS F30 file is not in the pool,
%                                               or no Swarm pass survived
%   NonGravForces:readSwarmDensity:*, NonGravForces:density:*, timeSeriesPreparation:*  propagated,
%                                               except readSwarmDensity's data* refusals (an empty
%                                               month), listed in emptyFiles
%
% EXAMPLE
%   R = densityDriverFit('Pool', "E:\DATAPOOL", 'First', "2019-01", 'Last', "2020-12");
%
% Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Pool (1,1) string = string(getenv("READING_DATA_POOL"))
    options.First (1,1) string = "2014-01"
    options.Last (1,1) string = "2026-08"
    options.Model (1,1) string = "msis00"
    options.Stride (1,1) double {mustBeInteger, mustBePositive} = 10
    options.MinRecords (1,1) double {mustBeInteger, mustBePositive} = 40
    options.Split (1,1) string = "2020-01-01"
    options.NBoot (1,1) double {mustBeInteger, mustBePositive} = 200
    options.MaxVif (1,1) double {mustBeGreaterThanOrEqual(options.MaxVif, 1)} = 10
    options.Cache (1,1) logical = true
    options.Missions (1,:) string {mustBeMember(options.Missions, ["swarm" "champ" "grace" "gracefo"])} = "swarm"
    options.Sets (1,:) string {mustBeMember(options.Sets, ["pass" "pass trimmed" "day" "day trimmed" "pass trimmed flux" "day trimmed flux"])} = ["day trimmed" "day trimmed flux"]
    options.FluxGrid (1,:) double {mustBeFinite} = -60:15:30
end
unit = "densityDriverFit";
swFile = fullfile(options.Pool, "SpaceWeather", "SW-All.csv");
f30File = fullfile(options.Pool, "SpaceWeather", "F30", "radio_flux_adjusted_observation.txt");
if ~isfile(swFile) || ~isfile(f30File)
    NonGravForces.fail(unit, "dataMissing", "%s or %s is not in the pool -> fetch them with setupNonGravForces('FetchPool', ...)", swFile, f30File);
end
sw = NonGravForces.io.spaceWeather(swFile);
cls = NonGravForces.io.spaceWeather(f30File);
t0 = tic;

% ---- 1. the observed factor per satellite, UTC day and node
% the sources: label, file of a month, months covered (the HAPI headers' spans)
src = struct('label', {}, 'file', {}, 'first', {}, 'last', {});
toleos = struct('mission', {"champ", "grace", "gracefo"}, 'label', {"CHAMP", "GRACE-A", "GRACE-FO1"}, ...
    'dataset', {"CH_OPER_DNS_ACC_2_", "GR_OPER_DNS1ACC_2_", "GF_OPER_DNS1ACC_2_"}, ...
    'first', {"2000-07", "2002-04", "2018-05"}, 'last', {"2010-09", "2017-10", "2100-01"});
if any(options.Missions == "swarm")
    for s = ["A" "B" "C"]
        src(end + 1) = struct('label', "Swarm " + s, 'file', @(m) fullfile(options.Pool, "Swarm", "DNSxPOD", string(year(m)), ...
            "SW_OPER_DNS" + s + "POD_2_" + string(m, 'yyyy-MM') + ".csv"), 'first', "2013-12", 'last', "2100-01"); %#ok<AGROW>
    end
end
for t = toleos(ismember([toleos.mission], options.Missions))
    src(end + 1) = struct('label', t.label, 'file', @(m) fullfile(options.Pool, "TOLEOS", t.dataset, string(year(m)), ...
        t.dataset + string(m, 'yyyy-MM') + ".csv"), 'first', t.first, 'last', t.last); %#ok<AGROW>
end
months = datetime(options.First, 'InputFormat', 'yyyy-MM'):calmonths(1):datetime(options.Last, 'InputFormat', 'yyyy-MM');
G = options.FluxGrid; nG = numel(G);
rows = zeros(0, 9 + nG);            % mjd, satellite, ascending, log f, height, LT, n, sum rho_obs, sum rho_model, sums at G
missingFiles = strings(1, 0); emptyFiles = strings(1, 0);
% The passes depend only on the files and these settings: kept in the pool (v3: with the sums at the flux grid).
cache = fullfile(options.Pool, "Swarm", "densityDriverFit_" + options.Model + "_" + options.First + "_" + options.Last + ...
    "_s" + options.Stride + "_n" + options.MinRecords + "_v3g" + strjoin(string(G), "_") + localMissionTag(options.Missions) + ".mat");
if options.Cache && isfile(cache)
    c = load(cache, 'rows', 'missingFiles', 'emptyFiles');
    rows = c.rows; missingFiles = c.missingFiles; emptyFiles = c.emptyFiles;
    months = months([]);
    fprintf('densityDriverFit: %d passes from %s\n', size(rows, 1), cache);
end
for m = months
    for is = 1:numel(src)
        if string(m, 'yyyy-MM') < src(is).first || string(m, 'yyyy-MM') > src(is).last, continue, end
        file = src(is).file(m);
        if ~isfile(file), missingFiles(end + 1) = file; continue, end %#ok<AGROW>
        try
            [ep, lla, rho, aux] = NonGravForces.io.readSwarmDensity(file);
        catch err
            if ~startsWith(err.identifier, "NonGravForces:readSwarmDensity:data"), rethrow(err); end
            emptyFiles(end + 1) = file; continue %#ok<AGROW>
        end
        dd = double(ep.day(:))'; ss = double(ep.sec(:))';
        lat = lla(1, :);
        asc = [diff(lat) > 0, false];
        asc(end) = asc(max(end - 1, 1));
        use = find(aux.valid & isfinite(rho) & rho > 0);
        use = use(1:options.Stride:end);
        if isempty(use), emptyFiles(end + 1) = file; continue, end %#ok<AGROW>
        e = referenceSystem.daySecondToEpoch([dd(use)' ss(use)'], "UTC");
        % the model at the points, and with F10.7 and F10.7A shifted by each grid value (the flux sets);
        % the shift 0 is the model itself, the same drivers sw gives
        q = sw.query(e);
        rhoG = zeros(nG, numel(use));
        for kg = 1:nG
            drv = struct('f107', max(q.f107 + G(kg), 1), 'f107a', max(q.f107a + G(kg), 1), 'ap', q.apHistory);
            rhoG(kg, :) = reshape(NonGravForces.env.density(options.Model, e, lla(:, use), [], 'Drivers', drv), 1, []);
        end
        if any(G == 0)
            rhoM = rhoG(find(G == 0, 1), :);
        else
            rhoM = reshape(NonGravForces.env.density(options.Model, e, lla(:, use), sw), 1, []);
        end
        [g, ~, grp] = unique(dd(use) * 2 + asc(use));
        for k = 1:numel(g)
            jm = find(grp == k)';
            j = use(jm);
            if numel(j) < options.MinRecords, continue, end
            eq = j(abs(lat(j)) < deg2rad(30));
            if isempty(eq), continue, end
            a = aux.lst_h(eq) * pi / 12;
            lt = mod(atan2(mean(sin(a)), mean(cos(a))) * 12 / pi, 24);
            rows(end + 1, :) = [floor(g(k) / 2), is, mod(g(k), 2), log(sum(rho(j)) / sum(rhoM(jm))), ...
                mean(lla(3, j)), lt, numel(j), sum(rho(j)), sum(rhoM(jm)), sum(rhoG(:, jm), 2)']; %#ok<AGROW>
        end
    end
    fprintf('densityDriverFit: %s  %6d passes  (%.0f s)\n', string(m, 'yyyy-MM'), size(rows, 1), toc(t0));
end
if isempty(rows)
    NonGravForces.fail(unit, "dataMissing", "no pass between %s and %s in %s -> fetch the months with setupNonGravForces('SwarmMonths' or 'DensityMonths', ...)", ...
        options.First, options.Last, options.Pool);
end
rows = sortrows(rows, [1 2 3]);     % time order: crossval's blocks are blocks of time
if options.Cache && ~isfile(cache)
    save(cache, 'rows', 'missingFiles', 'emptyFiles');
end

% ---- 2. the data sets: per pass (node) and per day (both nodes, ratio of the summed sums), each raw and
% trimmed. Trimmed drops a row whose log f is more than log 2 from its satellite's 27-day running median:
% the POD-only density is noise-dominated where the density is low (ESA's product notes), and one such row
% weighs in an RMS like a hundred good ones.
[gd, ~, gi] = unique(rows(:, [1 2]), 'rows');
sumObs = accumarray(gi, rows(:, 8)); sumMod = accumarray(gi, rows(:, 9)); nDay = accumarray(gi, rows(:, 7));
hDay = accumarray(gi, rows(:, 5) .* rows(:, 7)) ./ nDay;
sumG = zeros(size(gd, 1), nG);
for kg = 1:nG, sumG(:, kg) = accumarray(gi, rows(:, 9 + kg)); end
dayRows = [gd(:, 1), gd(:, 2), nan(size(gd, 1), 1), log(sumObs ./ sumMod), hDay, nan(size(gd, 1), 1), nDay, sumObs, sumMod, sumG];
passTrim = localTrim(rows); dayTrim = localTrim(dayRows);
sets = struct('name', {"pass", "pass trimmed", "day", "day trimmed", "pass trimmed flux", "day trimmed flux"}, ...
    'rows', {rows, passTrim, dayRows, dayTrim, localFlux(passTrim, G, "pass trimmed flux"), localFlux(dayTrim, G, "day trimmed flux")}, ...
    'lt', {true, true, false, false, true, false}, 'flux', {false, false, false, false, true, true});
names = ["Fbar" "dF" "ap" "sinAnnual" "cosAnnual" "sinSemiannual" "cosSemiannual" "sinLT" "cosLT" "height" "trend" "logFbar" "FbarSq" "phase"];
fprintf('\ndensityDriverFit %s..%s: %d passes (%s), missing files %d, empty files %d\n', options.First, options.Last, size(rows, 1), ...
    strjoin([src.label] + ":" + string(arrayfun(@(k) sum(rows(:, 2) == k), 1:numel(src))), ", "), numel(missingFiles), numel(emptyFiles));
out = struct(); records = struct();
for d = find(ismember([sets.name], options.Sets))
    v = localValidate(sets(d).rows, sets(d).name, sets(d).lt, sets(d).flux, names, [src.label], sw, cls, options);
    out.(matlab.lang.makeValidName(sets(d).name)) = v;
    records.(matlab.lang.makeValidName(sets(d).name)) = v.record;
end
if options.Cache
    recordFile = fullfile(options.Pool, "Swarm", "densityDriverFit_record_" + options.Model + "_" + options.First + "_" + options.Last + localMissionTag(options.Missions) + ".mat");
    save(recordFile, 'records');
    fprintf('  coefficient records of %s written to %s\n', strjoin(string(fieldnames(records))', ", "), recordFile);
end
fprintf('  %.0f s\n', toc(t0));

R.rows = array2table(rows(:, 1:9), 'VariableNames', ["mjd" "satellite" "ascending" "logF" "height_m" "lt_h" "n" "sumObs" "sumModel"]);
R.fluxGrid = G;
R.records = records;
R.drivers = names;
R.satellites = [src.label];
R.sets = out;
R.missingFiles = missingFiles;
R.emptyFiles = emptyFiles;
R.seconds = toc(t0);
end

function tag = localMissionTag(missions)
%LOCALMISSIONTAG The cache name's mission part: empty for Swarm alone, so the caches of 2026-10-27/28 stay valid.
tag = "";
if ~isequal(missions, "swarm"), tag = "_" + strjoin(sort(missions), "+"); end
end

function out = localFlux(rows, G, label)
%LOCALFLUX The flux rows: column 4 the shift delta of F10.7 and F10.7A (sfu) at which the model's summed
% density equals the observed sum (pchip in log rho over the grid G), the last column d log rho / d delta
% there. Columns 10.. (the sums at G) are dropped; a row outside the grid, or not monotone in it, is dropped.
nG = numel(G);
x = log(rows(:, 10:9 + nG)); y = log(rows(:, 8));
delta = nan(size(rows, 1), 1); slope = delta;
for k = 1:size(rows, 1)
    xk = x(k, :);
    if ~all(isfinite(xk)) || any(diff(xk) <= 0) || y(k) < xk(1) || y(k) > xk(end), continue, end
    delta(k) = interp1(xk, G, y(k), 'pchip');
    slope(k) = (interp1(G, xk, delta(k) + 1, 'pchip', 'extrap') - interp1(G, xk, delta(k) - 1, 'pchip', 'extrap')) / 2;
end
ok = isfinite(delta) & slope > 0;
fprintf('densityDriverFit %s: %d of %d rows inside the flux grid %s sfu (median delta %.1f sfu, median d log rho / d delta %.4f per sfu)\n', ...
    label, nnz(ok), numel(ok), mat2str(G), localQuantile(delta(ok), 0.5), localQuantile(slope(ok), 0.5));
out = [rows(ok, 1:3), delta(ok), rows(ok, 5:9), slope(ok)];
end

function rows = localTrim(rows)
%LOCALTRIM Drop rows more than log 2 from their satellite's 27-day running median of log f.
keep = true(size(rows, 1), 1);
for s = unique(rows(:, 2))'
    k = find(rows(:, 2) == s);
    t = rows(k, 1) + 0.25 * (rows(k, 3) == 1);        % the two nodes of a day apart: SamplePoints are unique
    med = movmedian(rows(k, 4), [13.5 13.5], 'SamplePoints', t);
    keep(k) = abs(rows(k, 4) - med) <= log(2);
end
rows = rows(keep, :);
end

function V = localValidate(rows, label, useLT, isFlux, names, satLabels, sw, cls, options)
%LOCALVALIDATE Drivers, blocked CV, forward/backward prediction and the full fits of one data set; printed.
% A flux set (isFlux) carries delta in column 4 and d log rho / d delta in column 10: its losses and
% residuals are multiplied by that slope, so they read in log f like the factor sets'.
mjd = rows(:, 1); y = rows(:, 4);
Fbar = nan(size(mjd)); dF = Fbar; ap = Fbar; F30bar = Fbar; dF30 = Fbar; ph107 = Fbar; ph30 = Fbar;
[ok, i] = ismember(mjd - 1, sw.mjd);
Fbar(ok) = sw.index.f107ObservedLast81(i(ok));
dF(ok) = sw.index.f107Observed(i(ok)) - Fbar(ok);
% the cycle's phase: the 81-day mean now minus a year earlier (> 0 rising, < 0 declining), past data only
[ok, i] = ismember(mjd - 366, sw.mjd);
ph107(ok) = Fbar(ok) - sw.index.f107ObservedLast81(i(ok));
[ok, i] = ismember(mjd, sw.mjd);
ap(ok) = sw.index.apDaily(i(ok));
f30 = cls.index.f30(:);
tr30 = movmean(f30, [80 0], 'omitnan');           % trailing: the reader's f30Mean81 is centred
[ok, i] = ismember(mjd - 1, cls.mjd(:));
F30bar(ok) = tr30(i(ok));
dF30(ok) = f30(i(ok)) - F30bar(ok);
[ok, i] = ismember(mjd - 366, cls.mjd(:));
ph30(ok) = F30bar(ok) - tr30(i(ok));
t = datetime(mjd + 2400000.5, 'ConvertFrom', 'juliandate');
ph = 2 * pi * (day(t, 'dayofyear') - 1) / 365.25;
lt = rows(:, 6) * pi / 12;
if ~useLT, lt = zeros(size(mjd)); end               % columns present, never selected
common = [sin(ph) cos(ph) sin(2 * ph) cos(2 * ph) sin(lt) cos(lt) (rows(:, 5) - 470e3) / 100e3 (mjd - 58849) / 365.25];
% NONLINEAR SOLAR TERMS (2026-10-28): the cycle-24 linear relation failed forward into cycle 25, and f is
% not monotonic in Fbar at high activity -- log Fbar, a centred square and the cycle's phase test why.
% The square is centred at EACH proxy's own mean over the data set (2026-10-28, PR #118): centred at 120 sfu
% for both, F30 -- which runs ~20 % lower than F10.7 -- made Fbar and its square collinear (VIF 9.5) and its
% forward fit was refused. A fixed constant per data set, so training and test rows see the same column.
ok107 = isfinite(Fbar) & isfinite(y); ok30 = isfinite(F30bar) & isfinite(y);
c107 = mean(Fbar(ok107)); c30 = mean(F30bar(ok30));
nl = @(F, P, c) [log(F), ((F - c) / 100) .^ 2, P / 100];
X = [y, Fbar, dF, ap, common, nl(Fbar, ph107, c107), F30bar, dF30, ap, common, nl(F30bar, ph30, c30)];   % 2..15 F10.7, 16..29 F30
wcol = 0; respName = "logF";
if isFlux, X = [X, rows(:, 10)]; wcol = size(X, 2); respName = "delta_sfu"; end                                                     % 30: the slope
keep = all(isfinite(X), 2);
X = X(keep, :); rows = rows(keep, :); mjd = mjd(keep);
if useLT
    nested = {1:3, 1:7, 1:10, 1:11}; nestName = ["solar+ap" "+season" "+LT+height" "+trend"];
else
    nested = {1:3, 1:7, [1:7 10], [1:7 10 11]}; nestName = ["solar+ap" "+season" "+height" "+trend"];
end
nested = [nested, {[12 2:7], [1:7 13], [1:7 14], [1:7 13 14], [12 2:7 14]}];
nestName = [nestName, "+season logF" "+season F^2" "+season phase" "+season F^2+phase" "+season logF+phase"];
spec = struct('name', "constant", 'cols', zeros(1, 0), 'n', 0);
for proxy = ["F10.7" "F30"]
    off = 1 + (proxy == "F30") * 14;
    for n = 1:numel(nested)
        spec(end + 1) = struct('name', proxy + " " + nestName(n), 'cols', off + nested{n}, 'n', n); %#ok<AGROW>
    end
end
S = numel(spec);
yr = year(datetime(mjd + 2400000.5, 'ConvertFrom', 'juliandate'));
folds = numel(unique(yr));
% The buffer is one solar rotation (27 days) of rows: the automatic Politis-White length would take the
% trend and the 81-day mean, perfectly persistent by construction, as the dependence range.
gap = ceil(27 * size(X, 1) / numel(unique(mjd)));
[cvL, cv] = timeSeriesPreparation.crossval(X, @(tr, te) localLoss(tr, te, spec, options.MaxVif, wcol), Folds = folds, Gap = gap);
split = juliandate(datetime(options.Split, 'InputFormat', 'yyyy-MM-dd')) - 2400000.5;
before = mjd < split;
fwd = nan(1, 2 * S); bwd = nan(1, 2 * S);
if any(before) && any(~before)
    tr = X; tr(~before, :) = NaN; te = X; te(before, :) = NaN;
    fwd = localLoss(tr, te, spec, options.MaxVif, wcol);
    tr = X; tr(before, :) = NaN; te = X; te(~before, :) = NaN;
    bwd = localLoss(tr, te, spec, options.MaxVif, wcol);
end
if isFlux
    fprintf('\n== %s: %d rows, delta median %.1f sfu, std %.1f, MAD %.1f; slope median %.4f log rho per sfu\n', label, size(X, 1), ...
        localQuantile(X(:, 1), 0.5), std(X(:, 1)), localQuantile(abs(X(:, 1) - localQuantile(X(:, 1), 0.5)), 0.5), localQuantile(X(:, wcol), 0.5));
else
    fprintf('\n== %s: %d rows, log f median %.3f (f %.3f), std %.3f, MAD %.3f\n', label, size(X, 1), localQuantile(X(:, 1), 0.5), ...
        exp(localQuantile(X(:, 1), 0.5)), std(X(:, 1)), localQuantile(abs(X(:, 1) - localQuantile(X(:, 1), 0.5)), 0.5));
end
fprintf('  the squares are centred at each proxy''s mean: F10.7 %.1f sfu, F30 %.1f sfu\n', c107, c30);
if isFlux
    fprintf('  per year: delta (sfu) at the 5 / 50 / 95 %% quantile, rows\n');
    for yy = unique(yr)'
        v = X(yr == yy, 1);
        fprintf('    %d  %6.1f  %6.1f  %6.1f  %6d\n', yy, localQuantile(v, 0.05), localQuantile(v, 0.5), localQuantile(v, 0.95), numel(v));
    end
else
    fprintf('  per year: f at the 5 / 50 / 95 %% quantile of log f, rows\n');
    for yy = unique(yr)'
        v = X(yr == yy, 1);
        fprintf('    %d  %.3f  %.3f  %.3f  %6d\n', yy, exp(localQuantile(v, 0.05)), exp(localQuantile(v, 0.5)), exp(localQuantile(v, 0.95)), numel(v));
    end
end
fprintf('  %-20s %8s %8s %8s %8s %8s %8s   (held-out log f; RMS and median |error|)\n', "model", "cvRMS", "cvMAE", "fwdRMS", "fwdMAE", "bwdRMS", "bwdMAE");
for s = 1:S
    fprintf('  %-20s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', spec(s).name, cvL(s), cvL(S + s), fwd(s), fwd(S + s), bwd(s), bwd(S + s));
end
fits = struct(); satTables = struct();
for proxy = ["F10.7" "F30"]
    % the full fit of the proxy's model with the smallest yearly-CV RMS
    own = find(startsWith([spec.name], proxy + " "));
    [~, b] = min(cvL(own));
    s = own(b);
    fprintf('  %s: smallest CV RMS is "%s"\n', proxy, spec(s).name);
    try
        r = timeSeriesPreparation.regressDrivers(X(:, 1), X(:, spec(s).cols), DriverNames = names(nested{spec(s).n}), ...
            MaxVif = options.MaxVif, NBoot = options.NBoot, Stage = "densityDriverFit." + label + "." + proxy);
        fprintf('  %s full fit: r2 %.3f; coefficient +- block-bootstrap sigma (%s per unit), VIF\n', proxy, r.r2, respName);
        fprintf('    %-14s %10.4f +- %.4f\n', "intercept", r.coef(1), r.sigma(1));
        for k = 1:numel(r.driverNames)
            fprintf('    %-14s %10.4g +- %.2g   VIF %.1f\n', r.driverNames(k), r.coef(k + 1), r.sigma(k + 1), r.vif(k));
        end
        % MISSION OFFSETS (2026-10-28): the full fit's residual per satellite. A mission whose median residual
        % stands apart from the others carries a calibration offset (POD-only Swarm against TU Delft's
        % accelerometer densities) that the solar terms would otherwise absorb.
        Z = X(:, spec(s).cols);
        res = X(:, 1) - (r.coef(1) + (Z - mean(Z, 1)) * r.coef(2:end));
        if isFlux, res = res .* X(:, wcol); end                % in log f: comparable with the factor sets
        % THE COEFFICIENT RECORD (2026-10-28) that env.densityDriverScale transcribes: the absolute form
        % log f = a0 + sum(c .* z), a0 = intercept - c . mean(z), every digit, and the fit's own prediction
        % at two reference days -- an oracle the unit's tests compare with, built by THIS code's drivers.
        zMean = mean(Z, 1);
        a0 = r.coef(1) - zMean * r.coef(2:end);
        fprintf('    coefficient record (%s, centre of the square %.17g sfu; response %s):\n', proxy, (proxy == "F10.7") * c107 + (proxy == "F30") * c30, respName);
        fprintf('      a0 = %.17g\n', a0);
        for k = 1:numel(r.driverNames)
            fprintf('      %-14s c = %.17g   mean = %.17g\n', r.driverNames(k), r.coef(k + 1), zMean(k));
        end
        for refDay = [54466 58849]
            j = find(mjd == refDay, 1);
            if isempty(j), continue, end
            fprintf('      reference MJD %d: z = %s  %s = %.17g\n', refDay, mat2str(Z(j, :), 17), respName, a0 + Z(j, :) * r.coef(2:end));
        end
        fprintf('    residual per satellite: rows, median, mean, RMS (log f)\n');
        sat = rows(:, 2);
        satRes = nan(numel(satLabels), 4);
        for k = 1:numel(satLabels)
            e = res(sat == k);
            if isempty(e), continue, end
            satRes(k, :) = [numel(e), localQuantile(e, 0.5), mean(e), sqrt(mean(e .^ 2))];
            fprintf('      %-10s %6d %+8.4f %+8.4f %8.4f\n', satLabels(k), satRes(k, :));
        end
        satTables.(matlab.lang.makeValidName(proxy)) = array2table(satRes, 'VariableNames', ["rows" "median" "mean" "rms"], 'RowNames', cellstr(satLabels));
    catch err
        if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
        r = struct('collinear', string(err.message));
        fprintf('  %s full fit refused: %s\n', proxy, r.collinear);
    end
    fits.(matlab.lang.makeValidName(proxy)) = r;
end
% THE RECORD of the fixed model F30 + season + square (the one env.densityDriverScale transcribes), whatever
% won the CV above: tools/dragScaleMonth.m reads the flux set's from the pool.
sRec = find([spec.name] == "F30 +season F^2", 1);
Zr = X(:, spec(sRec).cols);
try
    rr = timeSeriesPreparation.regressDrivers(X(:, 1), Zr, DriverNames = names(nested{spec(sRec).n}), MaxVif = options.MaxVif, NBoot = 2);
    V.record = struct('set', label, 'model', spec(sRec).name, 'response', respName, ...
        'driverNames', rr.driverNames, 'a0', rr.coef(1) - mean(Zr, 1) * rr.coef(2:end), 'coef', rr.coef(2:end), 'centre_sfu', c30, 'rows', size(X, 1));
    fprintf('  record %s / %s: a0 %.17g, coef %s, centre %.17g sfu\n', label, V.record.model, V.record.a0, mat2str(V.record.coef', 17), c30);
catch err
    if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
    V.record = struct('set', label, 'model', spec(sRec).name, 'response', respName, 'collinear', string(err.message));
    fprintf('  record %s / %s refused: %s\n', label, spec(sRec).name, err.message);
end
V.cv = table([spec.name]', cvL(1:S)', cvL(S + 1:end)', fwd(1:S)', fwd(S + 1:end)', bwd(1:S)', bwd(S + 1:end)', ...
    'VariableNames', ["model" "cvRms" "cvMae" "forwardRms" "forwardMae" "backwardRms" "backwardMae"]);
V.cvSe = cv.se;
V.fit = fits;
V.nRows = size(X, 1);
V.squareCentre = [c107 c30];
V.satelliteResidual = satTables;
end

function q = localQuantile(x, p)
%LOCALQUANTILE The p-quantile by sorting, linear between order statistics (prctile is the Statistics Toolbox).
x = sort(x(isfinite(x)));
if isempty(x), q = NaN; return, end
h = 1 + (numel(x) - 1) * p;
q = x(floor(h)) + (h - floor(h)) * (x(min(ceil(h), end)) - x(floor(h)));
end

function L = localLoss(tr, te, spec, maxVif, wcol)
%LOCALLOSS RMS and median absolute error of the held-out log f for each model, fitted on the training rows.
% wcol > 0: the response is a flux shift and column wcol its d log rho / d delta, which turns the error into log f.
S = numel(spec);
L = nan(1, 2 * S);
itr = all(isfinite(tr), 2); ite = all(isfinite(te), 2);
if ~any(itr) || ~any(ite), return, end
for s = 1:S
    ytr = tr(itr, 1);
    if isempty(spec(s).cols)
        yhat = repmat(mean(ytr), nnz(ite), 1);
    else
        Ztr = tr(itr, spec(s).cols);
        try
            r = timeSeriesPreparation.regressDrivers(ytr, Ztr, MaxVif = maxVif, NBoot = 2);   % the minimum regressDrivers takes; CV reads coef only
        catch err
            if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
            continue
        end
        yhat = r.coef(1) + (te(ite, spec(s).cols) - mean(Ztr, 1)) * r.coef(2:end);
    end
    e = te(ite, 1) - yhat;
    if wcol > 0, e = e .* te(ite, wcol); end
    L(s) = sqrt(mean(e .^ 2));
    L(S + s) = localQuantile(abs(e), 0.5);
end
end
