function [cal, info] = fitCalibration(obs, ref, options)
%FITCALIBRATION Estimate an accelerometer calibration from a filtered observation arc and a raw-rate model arc, through lsqvce.fit.
%
%   cal = NonGravForces.app.fitCalibration(obs, ref)
%   [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'Scale', "full", 'BiasBasis', "pspline", 'NBoot', 0)
%
%   The estimator of the calibration (N-02b; CONCEPT 6.2). Per axis i, with H
%   the CRN kernel of timeSeriesPreparation.crnKernel:
%
%       y_i(t) = H*b_i(t) + sum_j S_ij (H*a_j)(t) [+ K2_i (H*a_i^2)(t), demeaned] + e_i(t)
%
%   y is the observation `obs` (e.g. ACT1B at 1 Hz), a the model `ref` on the
%   kernel's raw grid (10 Hz), b_i(t) = B(t) c_i the bias curve on the basis
%   NonGravForces.app.Calibration.basis evaluates. EVERY deterministic column
%   is filtered by H exactly as the observation was -- the bias basis, the
%   model and the K2 regressor -- by timeSeriesPreparation.crnApply at the
%   observation epochs, so the fitted curve is the one the instrument carried
%   before the L1B filter (CONCEPT 6.2). The K2 column is H applied to a_i^2,
%   because the nonlinearity acts on the instrument's raw signal, before the
%   filter; (H a_i)^2 differs from it by the filter's smoothing of a square.
%   Its mean over the fitted rows goes into the bias, so K2 is reported with a
%   bias that carries -K2 * mean(H a_i^2).
%
%   THE THREE AXES ARE THREE PROBLEMS (CONCEPT 6.6), each one lsqvce.fit call
%   on declaration kind A: one system, one day, one arc per contiguous run of
%   usable rows (a gap or a NaN window of the filter ends an arc). An entry of
%   S that is not estimated is held at its a-priori value and moved to the
%   observation side. The columns are equilibrated to unit RMS before the
%   solve and the estimates scaled back, so a bias of 1e-7 m/s^2 and a scale
%   of 1 are solved at the same precision.
%
%   SCALE. "diag" estimates S_ii, "full" all nine, a 3x3 logical mask any
%   pattern. The cross-track scale S_yy is held at its a-priori value by
%   default ('FixCrossTrack', true; CONCEPT 6.2: the cross-track signal is too
%   small to calibrate a scale against). A row that estimates more than one
%   model column is GATED: when the variance inflation factor of any of its
%   model columns (timeSeriesPreparation.observability) exceeds 'VifMax', the
%   call refuses rather than report a shear the data cannot separate (v3 I.5-3).
%
%   NOISE. "white" (default) solves under unit weights; the formal sigmas
%   then assume white noise, which the CRN-filtered residual is not. That is
%   why 'NBoot' (default 200) adds a moving-block bootstrap of every parameter
%   through timeSeriesPreparation.bootstrap, refitting each replicate with
%   lsqvce.fit: the covariance the Calibration carries is the BOOTSTRAP one
%   when NBoot > 0, and uncertainty.inflation reports bootstrap/formal sigma.
%   "estimate" hands the noise to lsqvce's joint variance-component estimator
%   (Estimate "noise", the law of each axis estimated, AR whitener); a struct
%   is declared as the lsqvce covariance spec of every axis as given (e.g. an
%   'arma' law from timeSeriesPreparation.noise). The white default is a
%   decision, not a convenience: the filtered residual is band-limited at
%   0.035 Hz, so a whitening law estimated from it is ill-conditioned above
%   the band, while the block bootstrap needs no law at all.
%
%   EIV BOUND (v3 D4). The model carries the dominant error, so a regression
%   scale is attenuated toward zero by lambda = var(H a_j) / (var(H a_j) +
%   sigma_e^2) (Fuller 1987, sect. 1.1). With the residual variance of axis j,
%   divided by S_jj^2, as an UPPER bound on sigma_e^2, uncertainty.eiv.lambda
%   is a lower bound on the attenuation factor and scaleBiasBound = 1 - lambda
%   the worst-case relative bias of every scale on model axis j. It is a
%   bound, not a correction.
%
% INPUT ARGUMENTS
%   Name  Size  Type                        Units  Description
%   ----  ----  ----                        -----  -----------
%   obs   1x1   timeSeriesPreparation.Arc   m/s^2  Observation, three axes, units m/s2, already
%                                                  through the L1B chain (e.g. ACT1B); its
%                                                  epochs lie on ref's grid
%   ref   1x1   timeSeriesPreparation.Arc   m/s^2  Model accelerations, three axes, units m/s2, on
%                                                  the kernel's raw grid around obs (10 Hz for
%                                                  "ACC"); with 'Kernel' [] on obs's epochs
%
% NAME-VALUE ARGUMENTS
%   Name               Size  Type    Units  Description
%   ----               ----  ----    -----  -----------
%   Kernel             1x1   string  -      "crn" (default): crnKernel("ACC"); a crnKernel struct;
%                                           or [] for an unfiltered observation
%   BiasBasis          1x1   string  -      "constant" | "linear" (default) | "hat" | "pspline"
%   BiasKnotSpacing_s  1x1   double  s      Knot spacing of "hat" and "pspline"; default 10800.
%                                           Below the kernel span Tf it warns
%   Scale              3x3   string  -      "diag" (default) | "full" | a 3x3 logical mask
%   FixCrossTrack      1x1   logical -      Hold S_yy at its a-priori value; default true
%   ScaleApriori       3x3   double  -      Values of the entries not estimated; default eye(3)
%   Quadratic          1x1   logical -      Estimate K2 per axis; default false
%   Noise              1x1   string  -      "white" (default) | "estimate" | an lsqvce
%                                           covariance spec struct
%   NBoot              1x1   double  -      Bootstrap replicates, 0 for none; default 200
%   Seed               1x1   double  -      Bootstrap seed; default 0
%   VifMax             1x1   double  -      VIF gate for a row with several model columns;
%                                           default 10
%   MinWeightFraction  1x1   double  -      Minimum-data gate of observability; default 0.5
%
% OUTPUT ARGUMENTS
%   Name  Size  Type                            Units  Description
%   ----  ----  ----                            -----  -----------
%   cal   1x1   NonGravForces.app.Calibration   -      Convention "observation"; covScale and bias
%                                                      cov from the bootstrap when NBoot > 0, else
%                                                      formal; .uncertainty (.sigmaScaleFormal,
%                                                      .sigmaScaleBootstrap, .inflation 3x3,
%                                                      .residualRms 1x3 m/s^2, .eiv, .noise, .nBoot,
%                                                      .stopReason); .gates (.vif, .cond,
%                                                      .maxAbsCorr, .weightFraction, .vifMax,
%                                                      .verdict)
%   info  1x1   struct                          -      .residual (Arc at obs epochs, NaN where not
%                                                      fitted) .names (1x3 cell of column names)
%                                                      .rows (Nx3 logical) .Qxx (1x3 cell, the
%                                                      unit-weight cofactor of each axis, in the
%                                                      unscaled parameters) .results (1x3 cell of
%                                                      lsqvce.Results) .observability (1x3 cell)
%
% ACCURACY
%   The estimate is lsqvce's QR solve on equilibrated columns: rounding of
%   order cond * eps with cond the equilibrated condition number the gates
%   report. Each filtered column carries the rounding of one Nf-term dot
%   product (crnApply). The formal sigmas assume the declared noise; the
%   bootstrap sigmas are Monte-Carlo estimates (timeSeriesPreparation.
%   bootstrap's accuracy) with the block length of the Politis-White rule on
%   the residual of the axis. The
%   bias is a curve before the filter; Calibration.apply evaluates it
%   unfiltered (its help names that approximation).
%
% ERRORS
%   NonGravForces:fitCalibration:depMissing        lsqvce or timeSeriesPreparation is not on the path
%   NonGravForces:fitCalibration:inputArc          obs or ref not a three-axis Arc in m/s2
%   NonGravForces:fitCalibration:inputScale        a mask that is neither "diag", "full" nor 3x3,
%                                                  or a non-finite 'ScaleApriori'
%   NonGravForces:fitCalibration:inputGrid         'Kernel' [] and an obs epoch ref does not carry
%   NonGravForces:fitCalibration:inputKernel       a kernel that is neither "crn", a struct nor []
%   NonGravForces:fitCalibration:numericCollinear  the VIF gate refused a row with several model
%                                                  columns
%   NonGravForces:Calibration:inputBasis           knots the basis refuses
%   timeSeriesPreparation:*                        crnApply, observability, bootstrap
%   lsqvce:*                                       the solve
%
% WARNINGS
%   NonGravForces:fitCalibration:knotsBelowFilterSpan  knot spacing below the kernel span Tf: the
%                                                      filtered basis columns are nearly dependent
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md 6.2 (the observation equation, the
%   gates, the EIV bound) and 6.6 (three axes, three problems); the v3 plan
%   (docs/inputs/HANDOVER_AccelCal_v3.md) D4, I.5-3, WP-08, WP-09, WP-13;
%   Fuller W. (1987) Measurement Error Models, Wiley, sect. 1.1 (attenuation);
%   Kuensch H. (1989) Ann. Statist. 17:1217 (moving-block bootstrap);
%   Kutner M. et al. (2005) Applied Linear Statistical Models, ch. 10 (VIF);
%   Klinger B., Mayer-Guerr T. (2016) Adv. Space Res. 58:1597 (K2 term).
%
% EXAMPLE
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
%   cal = NonGravForces.app.fitCalibration(obs, ref, 'FixCrossTrack', false, 'NBoot', 50);
%   [diag(cal.scale), diag(truth.scale)]
%
% See also NonGravForces.app.Calibration, NonGravForces.app.syntheticPair, lsqvce.fit, timeSeriesPreparation.crnApply.
%
% KEYWORDS calibration accelerometer scale bias estimator lsqvce CRN bootstrap EIV VIF
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    obs (1,1)
    ref (1,1)
    options.Kernel = "crn"
    options.BiasBasis (1,1) string {mustBeMember(options.BiasBasis, ["constant" "linear" "hat" "pspline"])} = "linear"
    options.BiasKnotSpacing_s (1,1) double {mustBePositive, mustBeFinite} = 10800
    options.Scale = "diag"
    options.FixCrossTrack (1,1) logical = true
    options.ScaleApriori (3,3) double = eye(3)
    options.Quadratic (1,1) logical = false
    options.Noise = "white"
    options.NBoot (1,1) double {mustBeInteger, mustBeNonnegative} = 200
    options.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.VifMax (1,1) double {mustBePositive} = 10
    options.MinWeightFraction (1,1) double {mustBePositive, mustBeLessThanOrEqual(options.MinWeightFraction, 1)} = 0.5
end
unit = "fitCalibration";
if isempty(which('lsqvce.fit')) || isempty(which('timeSeriesPreparation.crnApply'))
    NonGravForces.fail(unit, "depMissing", "the estimate goes through lsqvce.fit and the filter through timeSeriesPreparation.crnApply -> run lsqvce_setup and setupTimeSeriesPreparation");
end
for a = {obs, ref}
    arc = a{1};
    if ~isa(arc, 'timeSeriesPreparation.Arc') || size(arc.obs, 2) ~= 3 || any(arc.units ~= "m/s2")
        NonGravForces.fail(unit, "inputArc", "obs and ref must each be one timeSeriesPreparation.Arc with three acceleration columns in m/s2 -> got a %s", class(arc));
    end
end
mask = localMask(options.Scale, unit);
if options.FixCrossTrack, mask(2, 2) = false; end
S0 = options.ScaleApriori;
if ~all(isfinite(S0(:)))
    NonGravForces.fail(unit, "inputScale", "'ScaleApriori' must be finite -> it supplies every entry not estimated");
end
K = localKernel(options.Kernel, unit);

% --- the bias basis, on the grid the filter reads -----------------------------------------
% The knots span every raw sample a filter window reads, half a kernel beyond the first and the
% last observation epoch, so the B-splines are a partition of unity wherever H looks.
t1 = obs.t(1); tN = obs.t(end);
tA = t1; tB = tN;
if ~isempty(K)
    half = int64(K.Nh) * int64(round(1e9 / K.fs));
    tA = t1 - half; tB = tN + half;
end
bias = struct('basis', options.BiasBasis, 't0', t1, 'knots', int64([]), 'coef', [], 'cov', {{[], [], []}});
h = int64(round(options.BiasKnotSpacing_s * 1e9));
switch options.BiasBasis
    case "hat"
        bias.knots = tA + int64(0:ceil(double(tB - tA) / double(h))) * h;
    case "pspline"
        bias.knots = tA + int64(-3:ceil(double(tB - tA) / double(h)) + 3) * h;
end
if ismember(options.BiasBasis, ["hat" "pspline"]) && ~isempty(K) && options.BiasKnotSpacing_s < K.Tf
    warning("NonGravForces:fitCalibration:knotsBelowFilterSpan", ...
        "knot spacing %g s is below the kernel span %g s: the filtered basis columns are nearly dependent -> raise 'BiasKnotSpacing_s'", options.BiasKnotSpacing_s, K.Tf);
end

% --- the filtered columns at the observation epochs ----------------------------------------
if isempty(K)
    [found, at] = ismember(obs.t, ref.t);
    if ~all(found)
        NonGravForces.fail(unit, "inputGrid", "with 'Kernel' [] every observation epoch must be a model epoch -> %d of %d are not", nnz(~found), numel(found));
    end
    Ha = ref.obs(at, :);
    HB = NonGravForces.app.Calibration.basis(bias, obs.t);
    Hq = ref.obs(at, :).^2;
else
    Ha = timeSeriesPreparation.crnApply(ref, K, OutputEpochs = obs.t).obs;
    Braw = NonGravForces.app.Calibration.basis(bias, ref.t);
    HB = timeSeriesPreparation.crnApply(localArc(ref.t, Braw, "b"), K, OutputEpochs = obs.t).obs;
    Hq = [];
    if options.Quadratic
        Hq = timeSeriesPreparation.crnApply(localArc(ref.t, ref.obs.^2, "q"), K, OutputEpochs = obs.t).obs;
    end
end
P = size(HB, 2);
bNames = "b_" + string(1:P);
axisNames = ["x" "y" "z"];

% --- three problems ----------------------------------------------------------------------
N = numel(obs.t);
usable = obs.weights > 0 & ~obs.filled;
scale = S0; covScale = zeros(9); coefB = zeros(3, P); covB = cell(1, 3);
quad = nan(3, 1);
sigF = nan(3); sigB = nan(3); resRms = nan(1, 3); stopReason = strings(1, 3);
res = nan(N, 3); rowsUsed = false(N, 3);
info = struct('residual', [], 'names', {cell(1, 3)}, 'rows', [], 'Qxx', {cell(1, 3)}, ...
    'results', {cell(1, 3)}, 'observability', {cell(1, 3)});
gates = struct('vif', {cell(1, 3)}, 'cond', nan(1, 3), 'maxAbsCorr', nan(1, 3), ...
    'weightFraction', nan(1, 3), 'vifMax', options.VifMax, 'verdict', strings(1, 3));
blockLen = nan(1, 3);
for i = 1:3
    est = find(mask(i, :)); fix = find(~mask(i, :));
    if options.Quadratic
        X = [HB, Ha(:, est), Hq(:, i)];
        names = [bNames, "S_" + axisNames(i) + axisNames(est), "K2_" + axisNames(i)];
    else
        X = [HB, Ha(:, est)];
        names = [bNames, "S_" + axisNames(i) + axisNames(est)];
    end
    l = obs.obs(:, i) - Ha(:, fix) * S0(i, fix).';
    ok = usable & isfinite(l) & all(isfinite(X), 2);
    if options.Quadratic
        X(:, end) = X(:, end) - mean(X(ok, end));
    end
    rowsUsed(:, i) = ok;

    % gate: identifiability before the fit (timeSeriesPreparation's metrics)
    yArc = localArc(obs.t, obs.obs(:, i), "y");
    Xg = X; Xg(~ok, :) = NaN;
    rep = timeSeriesPreparation.observability(yArc, Xg, Names = names, MinWeightFraction = options.MinWeightFraction);
    info.observability{i} = rep;
    gates.vif{i} = rep.vif(:, 1).'; gates.cond(i) = rep.cond(1);
    gates.maxAbsCorr(i) = rep.maxAbsCorr(1); gates.weightFraction(i) = rep.weightFraction(1);
    modelVif = rep.vif(P + (1:numel(est)), 1);
    if numel(est) > 1 && any(~(modelVif <= options.VifMax))
        NonGravForces.fail(unit, "numericCollinear", "axis %s estimates %d model columns and their VIF reaches %.3g > %g: the shear is not separable on this arc -> use 'Scale' ""diag"" or a longer arc", ...
            axisNames(i), numel(est), max(modelVif), options.VifMax);
    end
    gates.verdict(i) = "pass";
    if any(~(modelVif <= options.VifMax)), gates.verdict(i) = "warn"; end

    % the solve: lsqvce.fit, kind A, one arc per contiguous run
    [x, Qxx, Cxx, R] = localSolve(X, l, ok, obs.t, options.Noise);
    info.results{i} = R; info.Qxx{i} = Qxx; info.names{i} = names;
    stopReason(i) = string(R.stopReason);
    v = nan(N, 1); v(ok) = l(ok) - X(ok, :) * x;
    res(:, i) = v;
    resRms(i) = sqrt(mean(v(ok).^2));
    if options.NBoot > 0
        % The block length follows the dependence of the ERRORS, so it is the Politis-White rule on
        % the residual (bootstrap's own automatic rule, read off a two-replicate call), not on the
        % regressors, whose constant bias column has no autocorrelation to measure.
        [~, bl] = timeSeriesPreparation.bootstrap(v, @(r) mean(r, 1, "omitnan"), NBoot = 2, Stage = "fitCalibration.block");
        z = [l, X]; z(~ok, :) = NaN;
        [~, bi] = timeSeriesPreparation.bootstrap(z, @localBootStat, NBoot = options.NBoot, BlockLength = bl.blockLength, ...
            Seed = options.Seed, Stage = "fitCalibration." + axisNames(i));
        rows = all(isfinite(bi.replicates), 2);
        Cuse = cov(bi.replicates(rows, :));
        blockLen(i) = bi.blockLength;
    else
        Cuse = Cxx;
    end

    coefB(i, :) = x(1:P).';
    covB{i} = Cuse(1:P, 1:P);
    sIdx = P + (1:numel(est));
    scale(i, est) = x(sIdx).';
    lin = (est - 1) * 3 + i;                         % positions of S(i, est) in S(:)
    covScale(lin, lin) = Cuse(sIdx, sIdx);
    sigF(i, est) = sqrt(diag(Cxx(sIdx, sIdx))).';
    if options.NBoot > 0, sigB(i, est) = sqrt(diag(Cuse(sIdx, sIdx))).'; end
    if options.Quadratic, quad(i) = x(end); end
end

% --- the EIV bound per model axis --------------------------------------------------------
lambda = nan(1, 3);
for j = 1:3
    ok = rowsUsed(:, j);
    vs = var(Ha(ok, j));
    lambda(j) = vs / (vs + var(res(ok, j)) / scale(j, j)^2);
end

bias.coef = coefB; bias.cov = covB;
uncertainty = struct('sigmaScaleFormal', sigF, 'sigmaScaleBootstrap', sigB, 'inflation', sigB ./ sigF, ...
    'residualRms', resRms, 'eiv', struct('lambda', lambda, 'scaleBiasBound', 1 - lambda, ...
    'note', "lower bound on the attenuation factor, from the residual variance as an upper bound on the model error; a bound, not a correction"), ...
    'noise', localNoiseName(options.Noise), 'nBoot', options.NBoot, 'blockLength', blockLen, 'stopReason', stopReason);
prov = NonGravForces.provenance("app.fitCalibration", 'Options', localHashable(options), ...
    'Note', "obs " + obs.product + ", ref " + ref.product + ", kernel " + localKernelName(K));
cal = NonGravForces.app.Calibration(struct('scale', scale, 'covScale', covScale, 'bias', bias, ...
    'quadratic', quad, 'uncertainty', uncertainty, 'gates', gates, 'convention', "observation", 'provenance', prov));
info.residual = obs.derive(obs.t, res, "NonGravForces.app.fitCalibration", struct('basis', options.BiasBasis));
info.rows = rowsUsed;
end

% =========================================================================================
function [x, Qxx, Cxx, R] = localSolve(X, l, ok, t, noise)
%LOCALSOLVE One axis through lsqvce.fit, kind A: one system, one day, one arc per run of rows.
s = sqrt(mean(X(ok, :).^2, 1)); s(s == 0) = 1;
sl = sqrt(mean(l(ok).^2)); if sl == 0, sl = 1; end
Xe = X ./ s; le = l / sl;
idx = find(ok);
cut = [0; find(diff(idx) > 1); numel(idx)];
nArc = numel(cut) - 1;
D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
tSec = double(t - t(1)) / 1e9;
for k = 1:nArc
    r = idx(cut(k) + 1:cut(k + 1));
    D.A{k} = Xe(r, :); D.B{k} = zeros(numel(r), 0); D.l{k} = le(r); D.t{k} = tSec(r);
end
p = size(X, 2);
data = struct();
data.commonParams = [repmat({'common'}, p, 1), num2cell((1:p).')];
data.days = 1;
sys = struct('buildDay', @(day, obsType) D, 'covariance', struct('acc', localCovSpec(noise)));
data.systems = sys;
if isstruct(noise) || noise ~= "estimate"
    R = lsqvce.fit(data, lsqvce.Options(Verbose = false));
else
    o = lsqvce.Options(Estimate = "noise", Pool = "stream", ArcFactors = false, Verbose = false);
    o.Noise = struct('Whitener', "ar");
    R = lsqvce.fit(data, o);
end
% Back to the unscaled parameters. Q is the solve's own cofactor (sol.Cxx = sol.s0sq * Q), read
% from the accumulator rather than divided out of Cxx, which a noise-free arc leaves 0 or NaN.
x = R.x(:) * sl ./ s(:);
u = sl^2 ./ (s(:) * s(:).');
Qxx = R.acc.commonCofactor() .* u;
Cxx = R.sol.Cxx .* u;
end

function c = localCovSpec(noise)
%LOCALCOVSPEC The lsqvce covariance spec of the one observation type.
if isstruct(noise)
    c = noise;
elseif noise == "estimate"
    c = struct('type', 'diag', 'sigma', 1, 'estimate', "law");
else
    c = struct('type', 'diag', 'sigma', 1);
end
end

function x = localBootStat(z)
%LOCALBOOTSTAT One bootstrap replicate: the same equilibrated solve, through lsqvce.fit (kind B).
z = z(all(isfinite(z), 2), :);
l = z(:, 1); X = z(:, 2:end);
s = sqrt(mean(X.^2, 1)); s(s == 0) = 1;
sl = sqrt(mean(l.^2)); if sl == 0, sl = 1; end
R = lsqvce.fit(struct('A', X ./ s, 'l', l / sl), lsqvce.Options(Verbose = false));
x = (R.x(:) * sl ./ s(:)).';
end

function mask = localMask(scale, unit)
%LOCALMASK The estimated entries of S as a 3x3 logical.
if (ischar(scale) || isstring(scale)) && isscalar(string(scale))
    switch string(scale)
        case "diag", mask = logical(eye(3));
        case "full", mask = true(3);
        otherwise
            NonGravForces.fail(unit, "inputScale", "'Scale' must be ""diag"", ""full"" or a 3x3 logical mask -> got ""%s""", string(scale));
    end
elseif (islogical(scale) || isnumeric(scale)) && isequal(size(scale), [3 3])
    mask = logical(scale);
else
    NonGravForces.fail(unit, "inputScale", "'Scale' must be ""diag"", ""full"" or a 3x3 logical mask -> got a %s %s", mat2str(size(scale)), class(scale));
end
end

function K = localKernel(kernel, unit)
%LOCALKERNEL Resolve 'Kernel': the L1B ACC kernel, a given crnKernel struct, or none.
if isempty(kernel)
    K = [];
elseif (ischar(kernel) || isstring(kernel)) && string(kernel) == "crn"
    K = timeSeriesPreparation.crnKernel("ACC");
elseif isstruct(kernel) && isscalar(kernel) && all(isfield(kernel, {'F', 'fs', 'Nf', 'Nh', 'Tf', 'preset'}))
    K = kernel;
else
    NonGravForces.fail(unit, "inputKernel", "'Kernel' must be ""crn"", a timeSeriesPreparation.crnKernel struct or [] -> got a %s", class(kernel));
end
end

function arc = localArc(t, X, prefix)
%LOCALARC A design arc of named columns, the shape timeSeriesPreparation.Arc builds from.
names = cellstr(prefix + "_" + string(1:size(X, 2)));
u = struct();
for k = 1:numel(names), u.(names{k}) = '-'; end
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", u, "fileName", "DESIGN", "product", "DESIGN");
arc = timeSeriesPreparation.Arc(t, X, aux);
end

function s = localNoiseName(noise)
if isstruct(noise), s = "declared " + string(noise(1).type); else, s = string(noise); end
end

function s = localKernelName(K)
if isempty(K), s = "none"; else, s = string(K.preset); end
end

function o = localHashable(options)
%LOCALHASHABLE The options as provenance hashes them: a kernel struct by its preset only.
o = options;
if isstruct(o.Kernel), o.Kernel = "struct " + string(o.Kernel.preset); end
if isstruct(o.Noise), o.Noise = "declared " + string(o.Noise(1).type); end
end
