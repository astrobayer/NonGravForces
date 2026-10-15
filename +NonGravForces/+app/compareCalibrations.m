function [T, info] = compareCalibrations(obs, ref, specs, options)
%COMPARECALIBRATIONS Rank candidate calibration models by blocked cross-validation, per axis.
%
%   T = NonGravForces.app.compareCalibrations(obs, ref, specs)
%   [T, info] = NonGravForces.app.compareCalibrations(obs, ref, specs, 'Names', names, 'Folds', 5, 'Fit', {'Scale', "full"})
%
%   Model comparison for the accelerometer calibration (N-03b; the v3 plan's
%   WP-14 `compare`). Each entry of `specs` is one candidate model, written as
%   the name-value pairs NonGravForces.app.fitCalibration takes -- a bias
%   basis, a scale pattern, a time-varying scale, K2 -- and every candidate
%   is scored by the same yardstick: how well, fitted on part of the arc, it
%   PREDICTS the observation on the part it did not see.
%
%   BLOCKED, BECAUSE THE ROWS ARE DEPENDENT. timeSeriesPreparation.crossval
%   cuts the arc into 'Folds' contiguous blocks and, for each, fits on the
%   other blocks less a buffer of 'Gap' rows on either side, then scores the
%   held-out block (hv-block CV, Racine 2000). The buffer defaults to the
%   larger of crossval's automatic block length (Politis-White on the
%   observation) and the CRN window in rows: a filtered sample is a weighted
%   sum over 140.7 s of raw samples, so without that buffer a test row and
%   its training neighbours share raw data and the score rewards fitting
%   noise. ONE crossval call scores every candidate -- the loss function fits
%   all of them on the training rows -- so every candidate meets exactly the
%   same folds.
%
%   WHY NOT AIC OR BIC. Both count the n rows as n independent pieces of
%   information. The CRN-filtered residual is correlated over ~140 s, so the
%   effective sample size is far below n, the penalty term is far too small,
%   and both select models that are too flexible. They are reported, from
%   the full-arc fit, as reference columns whose variable descriptions say
%   IC_unreliable -- the v3 plan's requirement -- and never used to rank.
%
%   Each fit runs with 'NBoot' 0 (ranking needs the prediction, not its
%   uncertainty) and with 'Span' set to the whole arc's usable span, so a
%   B-spline basis is defined over a held-out edge block rather than zero
%   there. The 'Kernel' belongs to the observation and is common to all
%   candidates: pass it in 'Fit', never in a spec.
%
%   THE MINIMUM-DATA GATE IS SCALED ON A FOLD. app.fitCalibration refuses an
%   axis that keeps less than 'MinWeightFraction' (0.5) of its weight, and a
%   held-out block plus its buffers is weight the gate would count as lost --
%   with 2 folds, always more than half. So on a fold the fraction is the
%   candidate's own times the share of finite rows the fold keeps for
%   training: the gate still judges the DATA, not the cross-validation design.
%   The full-arc fit uses the fraction unchanged.
%
%   A LOCAL BASIS NEEDS KNOTS WIDER THAN A HELD-OUT BLOCK. A hat or B-spline
%   function whose support inside the arc falls wholly within a held-out
%   window has no training row, so the fit on that fold is singular -- not a
%   bad score, no score. At an edge fold that window is the block plus one
%   buffer plus half the kernel span, and a knot spacing above it also clears
%   every interior fold (supports of 2h and 4h exceed block + 2 buffers). So a
%   candidate whose bias or scale knot spacing does not exceed
%   (ceil(N/Folds) + Gap) rows + Tf/2 is REFUSED before any fit, rather than
%   failing mid-fold or, worse, skipping the folds it cannot predict and
%   scoring better for it.
%
% INPUT ARGUMENTS
%   Name   Size  Type                        Units  Description
%   ----   ----  ----                        -----  -----------
%   obs    1x1   timeSeriesPreparation.Arc   m/s^2  Observation, as for app.fitCalibration
%   ref    1x1   timeSeriesPreparation.Arc   m/s^2  Model on the kernel's raw grid, as for app.fitCalibration
%   specs  1xC   cell                        -      Candidate models: each a cell of app.fitCalibration
%                                                   name-value pairs ({} is the default model); C >= 1
%
% NAME-VALUE ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   Names  1xC   string  -      A name per candidate; default "spec1" ... "specC"
%   Fit    1xN   cell    -      Name-value pairs common to every candidate (a spec's own pair wins);
%                               'Kernel' belongs here; default {}
%   Folds  1x1   double  -      Contiguous blocks, integer >= 2; default 5
%   Gap    1x1   double  rows   Buffer on each side of a test block; default [] = the larger of
%                               crossval's automatic block length and the kernel span in rows
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   T     Cx8   table   -      One row per candidate: name (string), cvMse (1x3, m^2/s^4, mean over
%                              the folds of the test-row mean squared prediction error), cvMseSe
%                              (1x3, its standard error over the folds), rank (1x3, 1 = lowest
%                              cvMse on that axis), nParams (1x3), nRows (1x3, rows of the full-arc
%                              fit), aic and bic (1x3, from the full-arc fit; IC_unreliable).
%                              T.Properties.UserData: .best (1x3 string, the rank-1 name per axis)
%                              .folds .gap .note
%   info  1x1   struct  -      .foldLoss (1xC cell of Folds x 3) .testRows (Folds x 2) .gap
%                              .calibrations (1xC cell of the full-arc app.Calibration)
%
% ACCURACY
%   The score is exact arithmetic on the fits it composes (crossval partitions
%   exactly; each fit has app.fitCalibration's accuracy). As a statistic it is
%   one realisation: cvMseSe treats the fold losses as independent, which the
%   buffer makes approximately true (crossval's own caveat), so two
%   candidates within about two cvMseSe of each other are not distinguished
%   by this arc. AIC = n log(RSS/n) + 2k and BIC = n log(RSS/n) + k log n are
%   exact formulas applied to an n that overstates the information.
%
% ERRORS
%   NonGravForces:compareCalibrations:depMissing   timeSeriesPreparation is not on the path
%   NonGravForces:compareCalibrations:inputFolds   a candidate's knot spacing does not exceed the
%                                                  held-out window, so some fold cannot fit it
%   NonGravForces:compareCalibrations:inputSpecs   specs is not a non-empty cell of name-value cells,
%                                                  a spec sets 'Kernel', 'NBoot' or 'Span', or
%                                                  'Names' does not have one name per spec
%   NonGravForces:fitCalibration:*                 a candidate the estimator refuses
%   timeSeriesPreparation:crossval:*               folds that leave no training or test row
%   timeSeriesPreparation:blockLength:tooShort     the automatic gap on too short an arc
%
% REFERENCES
%   Racine J. (2000) J. Econometrics 99:39 (hv-block cross-validation);
%   Burman P., Chow E., Nolan D. (1994) Biometrika 81:351; Akaike H. (1974)
%   IEEE Trans. Autom. Control 19:716; Schwarz G. (1978) Ann. Stat. 6:461;
%   timeSeriesPreparation docs/CONCEPT.md 3.2 and 6.1 (row 7); the v3 plan
%   (docs/inputs/HANDOVER_AccelCal_v3.md) WP-14.
%
% EXAMPLE
%   [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
%   specs = {{'BiasBasis', "constant"}, {'BiasBasis', "linear"}, {'BiasBasis', "pspline", 'BiasKnotSpacing_s', 1800}};
%   T = NonGravForces.app.compareCalibrations(obs, ref, specs, 'Names', ["const" "linear" "spline"]);
%   T(:, ["name" "cvMse" "rank"]), T.Properties.UserData.best
%
% See also NonGravForces.app.fitCalibration, timeSeriesPreparation.crossval, NonGravForces.app.calibrate.
%
% KEYWORDS calibration model comparison cross-validation blocked hv-block AIC BIC accelerometer
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    obs (1,1)
    ref (1,1)
    specs
    options.Names string = string.empty
    options.Fit cell = {}
    options.Folds (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.Folds, 2)} = 5
    options.Gap double {mustBeScalarOrEmpty, mustBeInteger, mustBeNonnegative} = []
end
unit = "compareCalibrations";
if isempty(which('timeSeriesPreparation.crossval'))
    NonGravForces.fail(unit, "depMissing", "the folds come from timeSeriesPreparation.crossval -> run setupTimeSeriesPreparation");
end
if ~iscell(specs) || isempty(specs) || ~all(cellfun(@(c) iscell(c) && mod(numel(c), 2) == 0, specs))
    NonGravForces.fail(unit, "inputSpecs", "specs must be a non-empty cell of name-value cells, one per candidate -> e.g. {{'BiasBasis', ""constant""}, {}}");
end
C = numel(specs);
for c = 1:C
    keys = lower(string(specs{c}(1:2:end)));
    if any(ismember(keys, ["kernel" "nboot" "span"]))
        NonGravForces.fail(unit, "inputSpecs", "spec %d sets 'Kernel', 'NBoot' or 'Span'; the kernel is common (pass it in 'Fit') and the other two are set by the comparison -> remove them from the spec", c);
    end
end
names = options.Names;
if isempty(names), names = "spec" + string(1:C); end
if numel(names) ~= C
    NonGravForces.fail(unit, "inputSpecs", "'Names' must give one name per spec -> got %d names for %d specs", numel(names), C);
end
names = reshape(names, 1, C);

% The span every fit's bases must cover: the whole arc's usable epochs, as fitCalibration takes them.
live = find(obs.weights > 0 & ~obs.filled & any(isfinite(obs.obs), 2));
if isempty(live), live = 1:numel(obs.t); end
span = [obs.t(live(1)), obs.t(live(end))];
args = cell(1, C);
for c = 1:C
    args{c} = [localMerge(options.Fit, specs{c}), {'NBoot', 0, 'Span', span}];
end

% The buffer: the dependence of the observation, and at least the filter's reach.
gap = options.Gap;
if isempty(gap)
    gap = max(ceil(timeSeriesPreparation.internal.blockLength(obs.obs)));
    K = localKernel(options.Fit);
    if ~isempty(K)
        dt = double(median(diff(obs.t))) / 1e9;
        gap = max(gap, ceil(K.Tf / dt));
    end
end

% Every local basis must outlast the widest held-out window (see the help).
dt = double(median(diff(obs.t))) / 1e9;
half = 0;
Kc = localKernel(options.Fit);
if ~isempty(Kc), half = Kc.Tf / 2; end
window = (ceil(numel(obs.t) / options.Folds) + gap) * dt + half;
for c = 1:C
    for b = [struct('basis', "BiasBasis", 'h', "BiasKnotSpacing_s", 'b0', "linear", 'h0', 10800), ...
             struct('basis', "ScaleBasis", 'h', "ScaleKnotSpacing_s", 'b0', "constant", 'h0', 43200)]
        kind = string(localGet(args{c}, b.basis, b.b0));
        h = localGet(args{c}, b.h, b.h0);
        if ismember(kind, ["hat" "pspline"]) && h <= window
            NonGravForces.fail(unit, "inputFolds", "candidate %s: its %s %s with %g s knots has functions that a held-out window of %.0f s contains whole, so some fold cannot fit it -> raise %s above %.0f s or use more 'Folds'", ...
                names(c), kind, lower(extractBefore(b.basis, "Basis")), h, window, b.h, window);
        end
    end
end

loss = @(xTrain, xTest) localScore(obs, ref, args, xTrain, xTest);
[~, cv] = timeSeriesPreparation.crossval(obs.obs, loss, Folds = options.Folds, Gap = gap);

% Per candidate: the fold losses, their mean and standard error, and the full-arc fit for the ICs.
cvMse = nan(C, 3); cvSe = nan(C, 3); nPar = nan(C, 3); nRows = nan(C, 3); aic = nan(C, 3); bic = nan(C, 3);
foldLoss = cell(1, C); cals = cell(1, C);
for c = 1:C
    F = cv.foldLoss(:, (c - 1) * 3 + (1:3));
    foldLoss{c} = F;
    for i = 1:3
        v = F(isfinite(F(:, i)), i);
        if ~isempty(v), cvMse(c, i) = mean(v); end
        if numel(v) >= 2, cvSe(c, i) = std(v) / sqrt(numel(v)); end
    end
    [cals{c}, fi] = NonGravForces.app.fitCalibration(obs, ref, args{c}{:});
    for i = 1:3
        r = fi.residual.obs(fi.rows(:, i), i);
        n = numel(r); k = numel(fi.names{i});
        nPar(c, i) = k; nRows(c, i) = n;
        aic(c, i) = n * log(sum(r.^2) / n) + 2 * k;
        bic(c, i) = n * log(sum(r.^2) / n) + k * log(n);
    end
end
rnk = nan(C, 3);
best = strings(1, 3);
for i = 1:3
    [~, order] = sort(cvMse(:, i));                 % NaN sorts last
    rnk(order, i) = (1:C)';
    best(i) = names(order(1));
end

T = table(names', cvMse, cvSe, rnk, nPar, nRows, aic, bic, ...
    'VariableNames', {'name', 'cvMse', 'cvMseSe', 'rank', 'nParams', 'nRows', 'aic', 'bic'});
T.Properties.VariableUnits = {'', 'm^2/s^4', 'm^2/s^4', '', '', '', '', ''};
T.Properties.VariableDescriptions = {'candidate', 'blocked-CV mean squared prediction error per axis', ...
    'standard error of cvMse over the folds', '1 = lowest cvMse on that axis', 'parameters per axis', ...
    'rows of the full-arc fit', 'IC_unreliable: AIC of the full-arc fit, n overstates the information', ...
    'IC_unreliable: BIC of the full-arc fit, n overstates the information'};
T.Properties.UserData = struct('best', best, 'folds', options.Folds, 'gap', gap, ...
    'note', "ranked by blocked CV; aic and bic are IC_unreliable reference columns (correlated residual)");
info = struct('foldLoss', {foldLoss}, 'testRows', cv.testRows, 'gap', gap, 'calibrations', {cals});
end

% =========================================================================================
function v = localScore(obs, ref, args, xTrain, xTest)
%LOCALSCORE One fold: fit every candidate on the training rows, score it on the test rows.
C = numel(args);
v = nan(1, 3 * C);
obsTrain = obs.derive(obs.t, xTrain, "NonGravForces.app.compareCalibrations.train");
kept = nnz(any(isfinite(xTrain), 2)) / max(1, nnz(any(isfinite(obs.obs), 2)));
for c = 1:C
    a = localSet(args{c}, 'MinWeightFraction', kept * localGet(args{c}, 'MinWeightFraction', 0.5));
    [~, fi] = NonGravForces.app.fitCalibration(obsTrain, ref, a{:});
    r = xTest - fi.prediction.obs;
    for i = 1:3
        e = r(isfinite(r(:, i)), i);
        if ~isempty(e), v((c - 1) * 3 + i) = mean(e.^2); end
    end
end
end

function out = localMerge(common, own)
%LOCALMERGE Name-value pairs: the common ones, overridden by the candidate's own.
out = own;
ownKeys = lower(string(own(1:2:end)));
for k = 1:2:numel(common)
    if ~ismember(lower(string(common{k})), ownKeys)
        out = [out, common(k:k + 1)]; %#ok<AGROW>
    end
end
end

function args = localSet(args, name, value)
%LOCALSET Set a name-value pair (case-insensitive), replacing it if present.
for k = 1:2:numel(args)
    if strcmpi(string(args{k}), name), args{k + 1} = value; return; end
end
args = [args, {name, value}];
end

function v = localGet(args, name, default)
%LOCALGET The value of a name-value pair (case-insensitive), or the default.
v = default;
for k = 1:2:numel(args)
    if strcmpi(string(args{k}), name), v = args{k + 1}; end
end
end

function K = localKernel(common)
%LOCALKERNEL The kernel the fits will use: 'Kernel' from the common pairs, else the L1B ACC kernel.
K = timeSeriesPreparation.crnKernel("ACC");
for k = 1:2:numel(common)
    if strcmpi(string(common{k}), "Kernel")
        kv = common{k + 1};
        if isempty(kv), K = [];
        elseif isstruct(kv), K = kv;
        end
    end
end
end
