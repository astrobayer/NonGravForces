function R = densityDriverFit(options)
%DENSITYDRIVERFIT A driver model of the observed density scale: Swarm 2014-now, cross-validated over years.
%
%   R = densityDriverFit                                  % Swarm 2014-01..2026-08 from READING_DATA_POOL
%   R = densityDriverFit('First', "2019-01", 'Last', "2020-12")
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
%   NONLINEAR SOLAR TERMS (2026-10-28): log Fbar instead of Fbar, a centred square ((Fbar - 120)/100)^2,
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
%   IT LIVES IN tools/, like dragScaleMonth.m: it reads the DATAPOOL and runs for minutes.
%   .github/workflows/density-fit.yml runs it on the self-hosted runner.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   Pool        1x1   string  -      DATAPOOL root; default getenv("READING_DATA_POOL")
%   First, Last 1x1   string  -      "yyyy-mm" span of Swarm months; default 2014-01..2026-08
%   Model       1x1   string  -      env.density backend the factor is for; default "msis00"
%   Stride      1x1   double  -      every Stride-th valid 30 s record; default 10
%   MinRecords  1x1   double  -      fewer records in a pass drop it; default 40
%   Split       1x1   string  -      "yyyy-mm-dd" of the forward/backward split; default 2020-01-01
%   NBoot       1x1   double  -      bootstrap replicates of the full fits; default 200
%   MaxVif      1x1   double  -      regressDrivers' VIF gate; default 10
%   Sets        1xS   string  -      data sets validated: "pass", "pass trimmed", "day", "day trimmed"
%                                    (default the two trimmed ones; the raw ones are noise-dominated)
%   Cache       1x1   logical -      keep the passes in <Pool>/Swarm/densityDriverFit_<model>_<span>_s<Stride>_n<MinRecords>.mat
%                                    and reuse them (default true); false recomputes and does not write
%
% OUTPUT ARGUMENTS
%   R  1x1 struct: .rows (table: mjd, satellite, ascending, logF, height_m, lt_h, n, sumObs, sumModel),
%      .drivers (names), .sets.pass, .passTrimmed, .day, .dayTrimmed (each: .cv table of model, cvRms, cvMae,
%      forwardRms, forwardMae, backwardRms, backwardMae; .cvSe; .fit, the regressDrivers result of each proxy's
%      largest model; .nRows), .missingFiles, .emptyFiles, .seconds
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
    options.Sets (1,:) string {mustBeMember(options.Sets, ["pass" "pass trimmed" "day" "day trimmed"])} = ["pass trimmed" "day trimmed"]
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
sats = ["A" "B" "C"];
months = datetime(options.First, 'InputFormat', 'yyyy-MM'):calmonths(1):datetime(options.Last, 'InputFormat', 'yyyy-MM');
rows = zeros(0, 9);                 % mjd, satellite, ascending, log f, height, LT, n, sum rho_obs, sum rho_model
missingFiles = strings(1, 0); emptyFiles = strings(1, 0);
% The passes depend only on the files and these settings: kept in the pool (v2: with the two sums).
cache = fullfile(options.Pool, "Swarm", "densityDriverFit_" + options.Model + "_" + options.First + "_" + options.Last + ...
    "_s" + options.Stride + "_n" + options.MinRecords + "_v2.mat");
if options.Cache && isfile(cache)
    c = load(cache, 'rows', 'missingFiles', 'emptyFiles');
    rows = c.rows; missingFiles = c.missingFiles; emptyFiles = c.emptyFiles;
    months = months([]);
    fprintf('densityDriverFit: %d passes from %s\n', size(rows, 1), cache);
end
for m = months
    for is = 1:numel(sats)
        file = fullfile(options.Pool, "Swarm", "DNSxPOD", string(year(m)), "SW_OPER_DNS" + sats(is) + "POD_2_" + string(m, 'yyyy-MM') + ".csv");
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
        rhoM = reshape(NonGravForces.env.density(options.Model, e, lla(:, use), sw), 1, []);
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
                mean(lla(3, j)), lt, numel(j), sum(rho(j)), sum(rhoM(jm))]; %#ok<AGROW>
        end
    end
    fprintf('densityDriverFit: %s  %6d passes  (%.0f s)\n', string(m, 'yyyy-MM'), size(rows, 1), toc(t0));
end
if isempty(rows)
    NonGravForces.fail(unit, "dataMissing", "no Swarm pass between %s and %s in %s -> fetch the months with setupNonGravForces('SwarmMonths', ...)", ...
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
dayRows = [gd(:, 1), gd(:, 2), nan(size(gd, 1), 1), log(sumObs ./ sumMod), hDay, nan(size(gd, 1), 1), nDay, sumObs, sumMod];
sets = struct('name', {"pass", "pass trimmed", "day", "day trimmed"}, 'rows', {rows, localTrim(rows), dayRows, localTrim(dayRows)}, ...
    'lt', {true, true, false, false});
names = ["Fbar" "dF" "ap" "sinAnnual" "cosAnnual" "sinSemiannual" "cosSemiannual" "sinLT" "cosLT" "height" "trend" "logFbar" "FbarSq" "phase"];
fprintf('\ndensityDriverFit %s..%s: %d passes (%s), missing files %d, empty files %d\n', options.First, options.Last, size(rows, 1), ...
    strjoin(sats + ":" + string(arrayfun(@(k) sum(rows(:, 2) == k), 1:3)), " "), numel(missingFiles), numel(emptyFiles));
out = struct();
for d = find(ismember([sets.name], options.Sets))
    out.(matlab.lang.makeValidName(sets(d).name)) = localValidate(sets(d).rows, sets(d).name, sets(d).lt, names, sw, cls, options);
end
fprintf('  %.0f s\n', toc(t0));

R.rows = array2table(rows, 'VariableNames', ["mjd" "satellite" "ascending" "logF" "height_m" "lt_h" "n" "sumObs" "sumModel"]);
R.drivers = names;
R.sets = out;
R.missingFiles = missingFiles;
R.emptyFiles = emptyFiles;
R.seconds = toc(t0);
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

function V = localValidate(rows, label, useLT, names, sw, cls, options)
%LOCALVALIDATE Drivers, blocked CV, forward/backward prediction and the full fits of one data set; printed.
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
nl = @(F, P) [log(F), ((F - 120) / 100) .^ 2, P / 100];
X = [y, Fbar, dF, ap, common, nl(Fbar, ph107), F30bar, dF30, ap, common, nl(F30bar, ph30)];   % 2..15 F10.7, 16..29 F30
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
[cvL, cv] = timeSeriesPreparation.crossval(X, @(tr, te) localLoss(tr, te, spec, options.MaxVif), Folds = folds, Gap = gap);
split = juliandate(datetime(options.Split, 'InputFormat', 'yyyy-MM-dd')) - 2400000.5;
before = mjd < split;
fwd = nan(1, 2 * S); bwd = nan(1, 2 * S);
if any(before) && any(~before)
    tr = X; tr(~before, :) = NaN; te = X; te(before, :) = NaN;
    fwd = localLoss(tr, te, spec, options.MaxVif);
    tr = X; tr(before, :) = NaN; te = X; te(~before, :) = NaN;
    bwd = localLoss(tr, te, spec, options.MaxVif);
end
fprintf('\n== %s: %d rows, log f median %.3f (f %.3f), std %.3f, MAD %.3f\n', label, size(X, 1), localQuantile(X(:, 1), 0.5), ...
    exp(localQuantile(X(:, 1), 0.5)), std(X(:, 1)), localQuantile(abs(X(:, 1) - localQuantile(X(:, 1), 0.5)), 0.5));
fprintf('  per year: f at the 5 / 50 / 95 %% quantile of log f, rows\n');
for yy = unique(yr)'
    v = X(yr == yy, 1);
    fprintf('    %d  %.3f  %.3f  %.3f  %6d\n', yy, exp(localQuantile(v, 0.05)), exp(localQuantile(v, 0.5)), exp(localQuantile(v, 0.95)), numel(v));
end
fprintf('  %-20s %8s %8s %8s %8s %8s %8s   (held-out log f; RMS and median |error|)\n', "model", "cvRMS", "cvMAE", "fwdRMS", "fwdMAE", "bwdRMS", "bwdMAE");
for s = 1:S
    fprintf('  %-20s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', spec(s).name, cvL(s), cvL(S + s), fwd(s), fwd(S + s), bwd(s), bwd(S + s));
end
fits = struct();
for proxy = ["F10.7" "F30"]
    % the full fit of the proxy's model with the smallest yearly-CV RMS
    own = find(startsWith([spec.name], proxy + " "));
    [~, b] = min(cvL(own));
    s = own(b);
    fprintf('  %s: smallest CV RMS is "%s"\n', proxy, spec(s).name);
    try
        r = timeSeriesPreparation.regressDrivers(X(:, 1), X(:, spec(s).cols), DriverNames = names(nested{spec(s).n}), ...
            MaxVif = options.MaxVif, NBoot = options.NBoot, Stage = "densityDriverFit." + label + "." + proxy);
        fprintf('  %s full fit: r2 %.3f; coefficient +- block-bootstrap sigma (log f per unit), VIF\n', proxy, r.r2);
        fprintf('    %-14s %10.4f +- %.4f\n', "intercept", r.coef(1), r.sigma(1));
        for k = 1:numel(r.driverNames)
            fprintf('    %-14s %10.4g +- %.2g   VIF %.1f\n', r.driverNames(k), r.coef(k + 1), r.sigma(k + 1), r.vif(k));
        end
    catch err
        if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
        r = struct('collinear', string(err.message));
        fprintf('  %s full fit refused: %s\n', proxy, r.collinear);
    end
    fits.(matlab.lang.makeValidName(proxy)) = r;
end
V.cv = table([spec.name]', cvL(1:S)', cvL(S + 1:end)', fwd(1:S)', fwd(S + 1:end)', bwd(1:S)', bwd(S + 1:end)', ...
    'VariableNames', ["model" "cvRms" "cvMae" "forwardRms" "forwardMae" "backwardRms" "backwardMae"]);
V.cvSe = cv.se;
V.fit = fits;
V.nRows = size(X, 1);
end

function q = localQuantile(x, p)
%LOCALQUANTILE The p-quantile by sorting, linear between order statistics (prctile is the Statistics Toolbox).
x = sort(x(isfinite(x)));
if isempty(x), q = NaN; return, end
h = 1 + (numel(x) - 1) * p;
q = x(floor(h)) + (h - floor(h)) * (x(min(ceil(h), end)) - x(floor(h)));
end

function L = localLoss(tr, te, spec, maxVif)
%LOCALLOSS RMS and median absolute error of the held-out log f for each model, fitted on the training rows.
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
    L(s) = sqrt(mean(e .^ 2));
    L(S + s) = localQuantile(abs(e), 0.5);
end
end
