function R = densityDriverFit(options)
%DENSITYDRIVERFIT A driver model of the observed density scale: Swarm 2014-now, cross-validated over years.
%
%   R = densityDriverFit                                  % Swarm 2014-01..2026-09 from READING_DATA_POOL
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
%   First, Last 1x1   string  -      "yyyy-mm" span of Swarm months; default 2014-01..2026-09
%   Model       1x1   string  -      env.density backend the factor is for; default "msis00"
%   Stride      1x1   double  -      every Stride-th valid 30 s record; default 10
%   MinRecords  1x1   double  -      fewer records in a pass drop it; default 40
%   Split       1x1   string  -      "yyyy-mm-dd" of the forward/backward split; default 2020-01-01
%   NBoot       1x1   double  -      bootstrap replicates of the full fits; default 200
%   MaxVif      1x1   double  -      regressDrivers' VIF gate; default 10
%
% OUTPUT ARGUMENTS
%   R  1x1 struct: .rows (table: mjd, satellite, ascending, logF, height_m, lt_h, n),
%      .drivers (names), .cv (table: model, cvRms, cvSe, forwardRms, backwardRms),
%      .fit (struct per proxy: regressDrivers result of the largest model), .missingFiles,
%      .emptyFiles, .seconds
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
    options.Last (1,1) string = "2026-09"
    options.Model (1,1) string = "msis00"
    options.Stride (1,1) double {mustBeInteger, mustBePositive} = 10
    options.MinRecords (1,1) double {mustBeInteger, mustBePositive} = 40
    options.Split (1,1) string = "2020-01-01"
    options.NBoot (1,1) double {mustBeInteger, mustBePositive} = 200
    options.MaxVif (1,1) double {mustBeGreaterThanOrEqual(options.MaxVif, 1)} = 10
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
rows = zeros(0, 7);                 % mjd, satellite, ascending, log f, height, LT, n
missingFiles = strings(1, 0); emptyFiles = strings(1, 0);
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
                mean(lla(3, j)), lt, numel(j)]; %#ok<AGROW>
        end
    end
    fprintf('densityDriverFit: %s  %6d passes  (%.0f s)\n', string(m, 'yyyy-MM'), size(rows, 1), toc(t0));
end
if isempty(rows)
    NonGravForces.fail(unit, "dataMissing", "no Swarm pass between %s and %s in %s -> fetch the months with setupNonGravForces('SwarmMonths', ...)", ...
        options.First, options.Last, options.Pool);
end
rows = sortrows(rows, [1 2 3]);     % time order: crossval's blocks are blocks of time

% ---- 2. the drivers, from data that exist on and before each day
mjd = rows(:, 1); y = rows(:, 4);
Fbar = nan(size(mjd)); dF = Fbar; ap = Fbar; F30bar = Fbar; dF30 = Fbar;
[ok, i] = ismember(mjd - 1, sw.mjd);
Fbar(ok) = sw.index.f107ObservedLast81(i(ok));
dF(ok) = sw.index.f107Observed(i(ok)) - Fbar(ok);
[ok, i] = ismember(mjd, sw.mjd);
ap(ok) = sw.index.apDaily(i(ok));
f30 = cls.index.f30(:);
tr30 = movmean(f30, [80 0], 'omitnan');           % trailing: the reader's f30Mean81 is centred
[ok, i] = ismember(mjd - 1, cls.mjd(:));
F30bar(ok) = tr30(i(ok));
dF30(ok) = f30(i(ok)) - F30bar(ok);
t = datetime(mjd + 2400000.5, 'ConvertFrom', 'juliandate');
ph = 2 * pi * (day(t, 'dayofyear') - 1) / 365.25;
lt = rows(:, 6) * pi / 12;
common = [sin(ph) cos(ph) sin(2 * ph) cos(2 * ph) sin(lt) cos(lt) (rows(:, 5) - 470e3) / 100e3 (mjd - 58849) / 365.25];
names = ["Fbar" "dF" "ap" "sinAnnual" "cosAnnual" "sinSemiannual" "cosSemiannual" "sinLT" "cosLT" "height" "trend"];
X = [y, Fbar, dF, ap, common, F30bar, dF30, ap, common];   % columns 2..12 F10.7, 13..23 F30
keep = all(isfinite(X), 2);
X = X(keep, :); rows = rows(keep, :); mjd = mjd(keep);
nested = {1:3, 1:7, 1:10, 1:11};
nestName = ["solar+ap" "+season" "+LT+height" "+trend"];
spec = struct('name', "constant", 'cols', zeros(1, 0));
for proxy = ["F10.7" "F30"]
    off = 1 + (proxy == "F30") * 11;
    for n = 1:numel(nested)
        spec(end + 1) = struct('name', proxy + " " + nestName(n), 'cols', off + nested{n}); %#ok<AGROW>
    end
end
S = numel(spec);

% ---- 3. validation
yr = year(datetime(mjd + 2400000.5, 'ConvertFrom', 'juliandate'));
folds = numel(unique(yr));
% The buffer is one solar rotation (27 days) of passes: the automatic Politis-White length would
% take the trend and the 81-day mean, perfectly persistent by construction, as the dependence range.
gap = ceil(27 * size(X, 1) / numel(unique(mjd)));
[cvRms, cv] = timeSeriesPreparation.crossval(X, @(tr, te) localLoss(tr, te, spec, options.MaxVif), Folds = folds, Gap = gap);
split = juliandate(datetime(options.Split, 'InputFormat', 'yyyy-MM-dd')) - 2400000.5;
before = mjd < split;
fwd = nan(1, S); bwd = nan(1, S);
if any(before) && any(~before)
    tr = X; tr(~before, :) = NaN; te = X; te(before, :) = NaN;
    fwd = localLoss(tr, te, spec, options.MaxVif);
    tr = X; tr(before, :) = NaN; te = X; te(~before, :) = NaN;
    bwd = localLoss(tr, te, spec, options.MaxVif);
end
cvTable = table([spec.name]', cvRms', cv.se', fwd', bwd', 'VariableNames', ["model" "cvRms" "cvSe" "forwardRms" "backwardRms"]);

% ---- 4. full fits of each proxy's largest model, with bootstrap sigmas
fit = struct();
for proxy = ["F10.7" "F30"]
    s = find([spec.name] == proxy + " " + nestName(end));
    try
        fit.(matlab.lang.makeValidName(proxy)) = timeSeriesPreparation.regressDrivers(X(:, 1), X(:, spec(s).cols), ...
            DriverNames = names, MaxVif = options.MaxVif, NBoot = options.NBoot, Stage = "densityDriverFit." + proxy);
    catch err
        if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
        fit.(matlab.lang.makeValidName(proxy)) = struct('collinear', string(err.message));
    end
end

% ---- 5. report
fprintf('\ndensityDriverFit %s..%s: %d passes (%s), %d folds (one per year), log f mean %.3f (f = %.3f), std %.3f\n', ...
    options.First, options.Last, size(rows, 1), strjoin(sats + ":" + string(arrayfun(@(k) sum(rows(:, 2) == k), 1:3)), " "), ...
    folds, mean(X(:, 1)), exp(mean(X(:, 1))), std(X(:, 1)));
fprintf('  %-22s %9s %9s %11s %11s   (RMS of held-out log f, i.e. relative error of f)\n', "model", "cv", "+-se", "forward", "backward");
for s = 1:S
    fprintf('  %-22s %9.4f %9.4f %11.4f %11.4f\n', spec(s).name, cvRms(s), cv.se(s), fwd(s), bwd(s));
end
for proxy = string(fieldnames(fit))'
    r = fit.(proxy);
    if isfield(r, 'collinear'), fprintf('\n  %s full fit refused: %s\n', proxy, r.collinear); continue, end
    fprintf('\n  %s full fit: r2 %.3f, %d rows; coefficient +- block-bootstrap sigma (log f per unit), VIF\n', proxy, r.r2, r.nRows);
    fprintf('    %-14s %10.4f +- %.4f\n', "intercept", r.coef(1), r.sigma(1));
    for k = 1:numel(names)
        fprintf('    %-14s %10.4g +- %.2g   VIF %.1f\n', names(k), r.coef(k + 1), r.sigma(k + 1), r.vif(k));
    end
end
fprintf('\n  per year, observed mean f by satellite (A B C):\n');
for yy = unique(yr)'
    sel = yr == yy;
    fprintf('    %d  %s\n', yy, strjoin(compose("%.3f", arrayfun(@(k) exp(mean(X(sel & rows(:, 2) == k, 1))), 1:3)), "  "));
end
fprintf('  missing files %d, empty files %d, %.0f s\n', numel(missingFiles), numel(emptyFiles), toc(t0));

R.rows = array2table(rows, 'VariableNames', ["mjd" "satellite" "ascending" "logF" "height_m" "lt_h" "n"]);
R.drivers = names;
R.cv = cvTable;
R.fit = fit;
R.missingFiles = missingFiles;
R.emptyFiles = emptyFiles;
R.seconds = toc(t0);
end

function L = localLoss(tr, te, spec, maxVif)
%LOCALLOSS RMS of the held-out log f for each model, fitted on the training rows.
L = nan(1, numel(spec));
itr = all(isfinite(tr), 2); ite = all(isfinite(te), 2);
if ~any(itr) || ~any(ite), return, end
for s = 1:numel(spec)
    ytr = tr(itr, 1);
    if isempty(spec(s).cols)
        yhat = repmat(mean(ytr), nnz(ite), 1);
    else
        Ztr = tr(itr, spec(s).cols);
        try
            r = timeSeriesPreparation.regressDrivers(ytr, Ztr, MaxVif = maxVif, NBoot = 1);
        catch err
            if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
            continue
        end
        yhat = r.coef(1) + (te(ite, spec(s).cols) - mean(Ztr, 1)) * r.coef(2:end);
    end
    L(s) = sqrt(mean((te(ite, 1) - yhat) .^ 2));
end
end
