function [B, info] = forceBudget(models, epoch, r, v, q, options)
%FORCEBUDGET Per-effect per-axis acceleration budget over an arc, in chunks, for one or several models.
%
%   B = NonGravForces.app.forceBudget(model, epoch, r, v, q)
%   B = NonGravForces.app.forceBudget({fast, standard}, epoch, r, v, q, 'Labels', ["fast" "standard"])
%   [B, info] = NonGravForces.app.forceBudget(___, 'Chunk', 50, 'Frame', "gcrs")
%
%   THE QUESTION THIS ANSWERS is not "what is the acceleration" - Model.accel
%   answers that - but "which effect carries the arc, on which axis, and what
%   does a modelling choice cost". So the table carries a share of the total
%   beside every RMS, and several models may be reduced side by side in one
%   call, which is the comparison a budget is normally wanted for (turn ERP off,
%   turn self-shadowing off, change the preset, and read the difference).
%
%   WHAT IT ADDS OVER Model.budget. Model.budget evaluates every epoch of the
%   arc in ONE call and reduces it; a day at 1 Hz is 86400 epochs and the
%   per-panel state of all of them at once. This unit evaluates in chunks and
%   accumulates only the sums of squares and the running maxima, so the memory
%   is independent of the arc length, and it labels and stacks the models.
%   It is a reduction over Model.accel and computes no physics of its own
%   (principle P2 applied one layer up): with 'Chunk' at or above the epoch
%   count the result is BITWISE what Model.budget returns, and that identity is
%   the reference test.
%
%   THE CHUNKED SUM IS NOT BITWISE the unchunked one, and the help says so
%   rather than leaving a test to discover it: summing K squares in one pass and
%   summing them chunk by chunk differ in association. The maxima are exact
%   either way, and so are the shares formed from them.
%
% INPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   models   1x1|1xM Model   -      One NonGravForces.Model, or a cell array of them. Every
%                            or cell model must accept the same epochs; they may differ in preset,
%                                   options, macro model and mass
%   epoch    1x1     struct  -      referenceSystem epoch struct with K entries
%   r        3xK     double  m      Position, GCRS
%   v        3xK     double  m/s    Velocity, GCRS
%   q        4xK     double  -      BODY->GCRS quaternions, or [] to take them from each
%                                   model's attitude provider
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   Frame    1x1   string  -      "body" (default) or "gcrs"
%   Chunk    1x1   double  -      Epochs per evaluation, >= 1; default the first model's
%                                 Options.series.chunk. The result does not depend on it
%                                 beyond the rounding named under ACCURACY
%   Labels   1xM   string  -      One label per model, default the model's preset name with
%                                 an index appended when two share it
%
% OUTPUT ARGUMENTS
%   Name     Size      Type    Units  Description
%   ----     ----      ----    -----  -----------
%   B        (4M)x10   table   -      Four rows per model - drag, radiation, thermal, total -
%                                     with .model .effect .rms_x .rms_y .rms_z .max_x .max_y
%                                     .max_z .rms_norm .share (rms_norm over the model's total
%                                     rms_norm; 1 on the total row, and 0 throughout for a null
%                                     arc rather than 0/0). VariableUnits are m/s^2 on every
%                                     numeric column but share
%   info     1x1       struct  -      .K, .chunks, .chunkSize, .frame, .labels, .provenance
%
% ACCURACY
%   No arithmetic of its own beyond the reduction: every number comes from
%   Model.accel. With one chunk the table is bitwise Model.budget's. With
%   several chunks the RMS differs by the association of the sum alone, at the
%   1e-16 level relative; the maxima and the shares are exact. Shares are formed
%   from the rms_norm column, so they sum to more than 1 whenever two effects
%   partly cancel in the total - that is a property of the arc, not an error,
%   and the total row is reported so the reader can see it.
%
% ERRORS
%   NonGravForces:forceBudget:inputModels   models is not a Model or a cell array of Models
%   NonGravForces:forceBudget:inputShape    r, v or q do not match the epoch count, or Labels is not 1xM
%   NonGravForces:forceBudget:inputRange    Chunk < 1 or not finite
%   NonGravForces:Model:*                   propagated from the evaluation (missing Eop, attitude, weather)
%
% REFERENCES
%   HANDOVER.md Part G LP-7 (WP-30, moved there from LP-6 on 2026-09-08);
%   Woeske F. et al. (2019) Adv. Space Res. 63:1318 (the GRACE budget whose
%   shape this reproduces); Siemes C. et al. (2023) Adv. Space Res.
%
% EXAMPLE
%   S = NonGravForces.env.densityScale(swarmFiles, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));   % the standard preset's density scale
%   M = NonGravForces.Model(macro, NonGravForces.Options("standard"), 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
%   B = NonGravForces.app.forceBudget(M, ep, r, v, []);
%   B(B.effect ~= "total", {'effect', 'rms_norm', 'share'})
%
% See also NonGravForces.Model, NonGravForces.app.calibrate, NonGravForces.app.densityInversion.
%
% KEYWORDS budget force accounting per-effect RMS arc comparison
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
arguments
    models
    epoch (1,1) struct
    r double
    v double
    q double = []
    options.Frame (1,1) string {mustBeMember(options.Frame, ["body","gcrs"])} = "body"
    options.Chunk double = []
    options.Labels string = string.empty
end
unit = "forceBudget";
if isa(models, 'NonGravForces.Model')
    mdl = num2cell(reshape(models, 1, []));
elseif iscell(models) && ~isempty(models) && all(cellfun(@(m) isa(m, 'NonGravForces.Model'), models(:)))
    mdl = reshape(models, 1, []);
else
    NonGravForces.fail(unit, "inputModels", "models must be a NonGravForces.Model or a non-empty cell array of them -> got %s", class(models));
end
M = numel(mdl);
K = numel(epoch.sec);
if ~(ismatrix(r) && isequal(size(r), [3 K]) && isequal(size(v), [3 K]))
    NonGravForces.fail(unit, "inputShape", "r and v must be 3x%d for this epoch -> got %s and %s", K, mat2str(size(r)), mat2str(size(v)));
end
if ~isempty(q) && ~isequal(size(q), [4 K])
    NonGravForces.fail(unit, "inputShape", "q must be 4x%d or [] -> got %s", K, mat2str(size(q)));
end
if isempty(options.Chunk)
    chunk = mdl{1}.options.series.chunk;
else
    chunk = options.Chunk;
end
if ~isscalar(chunk) || ~isfinite(chunk) || chunk < 1
    NonGravForces.fail(unit, "inputRange", "'Chunk' must be a finite scalar >= 1 -> got %s", mat2str(chunk));
end
chunk = max(1, round(chunk));
if isempty(options.Labels)
    labels = localLabels(mdl);
else
    labels = reshape(options.Labels, 1, []);
    if ~isequal(size(labels), [1 M])
        NonGravForces.fail(unit, "inputShape", "'Labels' must hold one label per model (%d) -> got %d", M, numel(labels));
    end
end

eff = ["drag"; "radiation"; "thermal"; "total"];
nE = numel(eff);
first = 1:chunk:K;
last = min(first + chunk - 1, K);
nC = numel(first);
sumsq = zeros(nE, 3, M);
mx = zeros(nE, 3, M);
for m = 1:M
    for c = 1:nC
        idx = first(c):last(c);
        ep = referenceSystem.daySecondToEpoch([epoch.day(idx) epoch.sec(idx)], epoch.scale);
        qc = []; if ~isempty(q), qc = q(:, idx); end
        [a, ai] = mdl{m}.accel(ep, r(:, idx), v(:, idx), qc, 'Frame', options.Frame);
        vals = {ai.drag, ai.radiation, ai.thermal, a};
        for i = 1:nE
            sumsq(i, :, m) = sumsq(i, :, m) + sum(vals{i}.^2, 2)';
            mx(i, :, m) = max(mx(i, :, m), max(abs(vals{i}), [], 2)');
        end
    end
end
rms = sqrt(sumsq / K);                                  % mean is sum/K, so one chunk is bitwise Model.budget
rmsNorm = reshape(vecnorm(rms, 2, 2), nE, M);
modelCol = strings(nE * M, 1); effCol = strings(nE * M, 1);
val = zeros(nE * M, 7); share = zeros(nE * M, 1);
for m = 1:M
    rows = (m - 1) * nE + (1:nE);
    modelCol(rows) = labels(m);
    effCol(rows) = eff;
    val(rows, :) = [rms(:, :, m), mx(:, :, m), rmsNorm(:, m)];
    tot = rmsNorm(nE, m);                               % the total row is the last effect
    if tot > 0
        share(rows) = rmsNorm(:, m) / tot;
    end                                                 % a null arc leaves every share at 0 rather than 0/0
end
B = table(modelCol, effCol, val(:, 1), val(:, 2), val(:, 3), val(:, 4), val(:, 5), val(:, 6), val(:, 7), share, ...
    'VariableNames', {'model', 'effect', 'rms_x', 'rms_y', 'rms_z', 'max_x', 'max_y', 'max_z', 'rms_norm', 'share'});
B.Properties.VariableUnits = ["", "", repmat("m/s^2", 1, 7), ""];
B.Properties.Description = sprintf('NonGravForces.app.forceBudget, %s frame, %d epochs in %d chunk(s), %d model(s)', ...
    options.Frame, K, nC, M);
info = struct('K', K, 'chunks', nC, 'chunkSize', chunk, 'frame', options.Frame, 'labels', labels, ...
    'provenance', NonGravForces.provenance("app.forceBudget", 'Options', struct('Frame', options.Frame, 'Chunk', chunk, 'models', M)));
end

function labels = localLabels(mdl)
%LOCALLABELS The preset of each model, made unique by an index when two share one.
M = numel(mdl);
labels = strings(1, M);
for m = 1:M
    labels(m) = mdl{m}.options.preset;
end
dup = labels;
for m = 1:M
    if nnz(dup == dup(m)) > 1
        labels(m) = dup(m) + "#" + m;
    end
end
end
