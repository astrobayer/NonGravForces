function [cal, info] = calibrate(yMeas, aModel, epoch, options)
%CALIBRATE Accelerometer bias, scale and jumps from a measured/modelled pair, by weighted least squares.
%
%   cal = NonGravForces.app.calibrate(yMeas, aModel, epoch)
%   cal = NonGravForces.app.calibrate(yMeas, aModel, epoch, 'Bias', "hat", 'KnotInterval_s', 3600)
%   [cal, info] = NonGravForces.app.calibrate(___, 'Scale', "full", 'Events', T, 'Robust', "huber")
%
%   THE MODEL BEING FITTED is the one sensor.Electrostatic implements and
%   qc.residuals applies, so a calibration estimated here goes straight into
%   either without a change of convention:
%
%       yMeas(:,k) = bias(:,k) + S aModel(:,k) + e(:,k)
%       bias(i,k)  = sum_j c(i,j) phi_j(t_k)  +  sum_m d(i,m) H(t_k - t_m)
%
%   with phi a basis in time and H the unit step at each declared jump.
%
%   IT SEPARATES BY AXIS, and that is the whole reason no general estimator is
%   needed here. Row i of yMeas depends on row i of the bias coefficients and row
%   i of S and on nothing else, so the problem is THREE independent least squares
%   of a few dozen unknowns rather than one of a few hundred. Each is solved by
%   QR on the design matrix -- never by forming the normal equations, which
%   squares the condition number of a problem that is already badly scaled.
%
%   THE COLUMNS ARE EQUILIBRATED BEFORE THE SOLVE, and this is not cosmetic: the
%   basis columns are of order 1 and the scale column is an acceleration of order
%   1e-7, so the raw design matrix has a condition number of 1e7 from the units
%   alone. Each column is divided by its norm, the scaled system is solved, and
%   the solution and its covariance are transformed back. What `cal.condition`
%   reports is the condition of the EQUILIBRATED matrix, which is the one that
%   says something about the geometry of the problem rather than about metres per
%   second squared.
%
%   WHAT IT DOES NOT DO, and where the limits are: no variance components, no
%   spectral noise model, no correlated weighting. The weights are per axis and
%   per sample, diagonal. Accelerometer noise is emphatically NOT white -- it
%   rises as 1/f below a few mHz -- so a bias basis fine enough to absorb the
%   drift is doing the work a coloured covariance would otherwise do, and the
%   formal sigmas below are consequently OPTIMISTIC. They are reported because a
%   relative comparison between axes and between arcs is still informative; they
%   are not a claim about absolute accuracy. Whiten with qc.residuals first if an
%   honest sigma matters.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units   Description
%   ----     ----  ----    -----   -----------
%   yMeas    3xK   double  m/s^2   Measured accelerations, sensor frame
%   aModel   3xK   double  m/s^2   Modelled non-gravitational accelerations, same frame
%   epoch    1x1   struct  -       referenceSystem epoch struct with K entries; only differences
%                                  of its seconds are used, so any scale is fine
%
% NAME-VALUE ARGUMENTS
%   Name             Size  Type     Units  Description
%   ----             ----  ----     -----  -----------
%   Bias             1x1   string   -      "constant" (default) | "linear" | "hat". "hat" is a
%                                          continuous piecewise-linear basis on a uniform knot
%                                          grid: local support, a partition of unity, and no end
%                                          conditions to get wrong. A cubic B-spline basis is NOT
%                                          offered, because a drifting bias is not required to
%                                          have continuous curvature and the smoothness would be
%                                          bought with three times the code
%   KnotInterval_s   1x1   double   s      Knot spacing for "hat", > 0; default 3600
%   Scale            1x1   string   -      "diag" (default) | "full" | "none". "diag" estimates
%                                          S(i,i), "full" the whole row i of S (axis couplings)
%   Events           Nx1   double   s      Seconds from the first epoch at which the bias steps;
%                                          one extra unknown per event per axis. Take them from
%                                          qc.eventTable. Default none
%   Weights          3xK   double   -      Per-sample weights (1/sigma^2 up to a common factor),
%                                          or 1x3 for one per axis; default all ones
%   Mask             3xK   logical  -      false where a sample must not be used; default all true
%   Robust           1x1   string   -      "none" (default) | "huber". IRLS with the Huber weight
%   HuberK           1x1   double   -      Huber threshold in units of the robust sigma; default 1.345
%   MaxIterations    1x1   double   -      IRLS cap, >= 1; default 10
%   MaxCondition     1x1   double   -      Refuse above this condition of the equilibrated design
%                                          matrix; default 1e10
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units   Description
%   ----   ----  ----    -----   -----------
%   cal    1x1   struct  -       .bias 3xK evaluated bias, .biasCoefficients 3xnB,
%                                .jump 3xNe step sizes, .scale 3x3 (identity where not
%                                estimated), .sigma0 3x1 a-posteriori sigma per axis,
%                                .sigma 3x(nB+ns+Ne) formal parameter sigmas,
%                                .covariance 1x3 cell of parameter covariances,
%                                .names 1xP parameter names, .residual 3xK,
%                                .dof 3x1, .condition 3x1, .iterations 3x1, .used 3xK logical
%   info   1x1   struct  -       .K .nBasis .nScale .nEvents .basis .provenance
%
% ACCURACY
%   Planted parameters are recovered exactly in the noise-free case: the model is
%   linear, so a design matrix of full rank returns the truth to the conditioning
%   of the QR and nothing else. With noise the recovery is unbiased and the
%   formal sigmas scale as 1/sqrt(K) -- both asserted in tests/tNgfApp. The
%   caveat under WHAT IT DOES NOT DO applies to the sigmas, not to the estimates.
%
% ERRORS
%   NonGravForces:calibrate:inputShape      yMeas, aModel, Weights or Mask do not match the epoch count
%   NonGravForces:calibrate:inputRange      KnotInterval_s <= 0, MaxIterations < 1, a weight < 0,
%                                           or an event outside the arc
%   NonGravForces:calibrate:numericNaN      an unmasked sample is NaN or Inf
%   NonGravForces:calibrate:numericSingular the equilibrated design matrix exceeds MaxCondition,
%                                           or an axis has fewer usable samples than unknowns
%
% REFERENCES
%   Klinger B., Mayer-Guerr T. (2016) Adv. Space Res. 58:1597 (bias and scale
%   parametrisation for GRACE); Behzadpour S. et al. (2021) JGR Solid Earth
%   126:e2020JB021297 (daily scale on a fine knot grid); Bjoerck A. (1996)
%   Numerical Methods for Least Squares Problems, sect. 2.4 (QR, not the normal
%   equations) and 7.3 (column equilibration); Huber P. (1981) Robust Statistics.
%   HANDOVER.md Part G LP-7 (WP-19).
%
% EXAMPLE
%   [~, ai] = M.accel(ep, r, v, q);
%   cal = NonGravForces.app.calibrate(measured, ai.drag + ai.radiation + ai.thermal, ep, ...
%             'Bias', "hat", 'KnotInterval_s', 3600, 'Scale', "diag");
%   cal.scale                        % the along-track entry is the published ~1.05 for GRACE-C
%   s = NonGravForces.sensor.Electrostatic(cal.bias(:,1), cal.scale, [0;0;0]);
%
% See also NonGravForces.sensor.Electrostatic, NonGravForces.qc.residuals, NonGravForces.qc.eventTable, NonGravForces.app.densityInversion.
%
% KEYWORDS calibration bias scale factor least squares jumps robust accelerometer
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    yMeas double
    aModel double
    epoch (1,1) struct
    options.Bias (1,1) string {mustBeMember(options.Bias, ["constant","linear","hat"])} = "constant"
    options.KnotInterval_s (1,1) double = 3600
    options.Scale (1,1) string {mustBeMember(options.Scale, ["none","diag","full"])} = "diag"
    options.Events double = zeros(0, 1)
    options.Weights double = []
    options.Mask logical = logical.empty
    options.Robust (1,1) string {mustBeMember(options.Robust, ["none","huber"])} = "none"
    options.HuberK (1,1) double = 1.345
    options.MaxIterations (1,1) double = 10
    options.MaxCondition (1,1) double = 1e10
end
unit = "calibrate";
K = numel(epoch.sec);
if ~isequal(size(yMeas), [3 K]) || ~isequal(size(aModel), [3 K])
    NonGravForces.fail(unit, "inputShape", "yMeas and aModel must be 3x%d for this epoch -> got %s and %s", K, mat2str(size(yMeas)), mat2str(size(aModel)));
end
mask = options.Mask;
if isempty(mask), mask = true(3, K); end
if ~isequal(size(mask), [3 K])
    NonGravForces.fail(unit, "inputShape", "'Mask' must be 3x%d or [] -> got %s", K, mat2str(size(mask)));
end
W = options.Weights;
if isempty(W), W = ones(3, K); elseif isequal(size(W), [1 3]) || isequal(size(W), [3 1]), W = repmat(W(:), 1, K); end
if ~isequal(size(W), [3 K])
    NonGravForces.fail(unit, "inputShape", "'Weights' must be 3x%d, 1x3 or [] -> got %s", K, mat2str(size(options.Weights)));
end
if any(W(:) < 0) || any(~isfinite(W(:)))
    NonGravForces.fail(unit, "inputRange", "weights must be finite and non-negative -> a negative weight is not a down-weighting, it flips the residual's sign");
end
if any(~isfinite(yMeas(mask))) || any(~isfinite(aModel(mask)))
    NonGravForces.fail(unit, "numericNaN", "an UNMASKED sample of yMeas or aModel is NaN/Inf -> mask the gaps out instead of letting them into the normal matrix");
end
if ~(options.KnotInterval_s > 0) || ~isfinite(options.KnotInterval_s)
    NonGravForces.fail(unit, "inputRange", "'KnotInterval_s' must be finite and positive -> it is the spacing of the hat knots in seconds");
end
if options.MaxIterations < 1 || ~isfinite(options.MaxIterations)
    NonGravForces.fail(unit, "inputRange", "'MaxIterations' must be a finite number >= 1 -> got %g", options.MaxIterations);
end

% Seconds from the first epoch. Only differences enter, so the time scale of the
% epoch struct is irrelevant and no conversion is done for it.
t = (epoch.day(:)' - epoch.day(1)) * 86400 + (epoch.sec(:)' - epoch.sec(1));
span = max(t) - min(t);
ev = options.Events(:)';
if ~isempty(ev) && (any(ev <= min(t)) || any(ev >= max(t)))
    NonGravForces.fail(unit, "inputRange", "every event must lie strictly inside the arc [%g, %g] s -> a step at an endpoint is a constant, which the bias basis already carries", min(t), max(t));
end

% --- bias basis -------------------------------------------------------------
switch options.Bias
    case "constant"
        Phi = ones(K, 1);
        bNames = "bias0";
    case "linear"
        % Normalised to [-1, 1]: a raw t in seconds would put 1e4 beside 1 and
        % cost four digits before the solve even starts.
        tau = zeros(1, K);
        if span > 0, tau = 2 * (t - min(t)) / span - 1; end
        Phi = [ones(K, 1), tau(:)];
        bNames = ["bias0", "biasRate"];
    otherwise
        knots = min(t):options.KnotInterval_s:max(t);
        if knots(end) < max(t), knots(end+1) = knots(end) + options.KnotInterval_s; end
        d = options.KnotInterval_s;
        Phi = max(0, 1 - abs(t(:) - knots) / d);        % K x nKnots hat functions
        bNames = "biasKnot" + string(1:numel(knots));
end
nB = size(Phi, 2);
% --- event steps ------------------------------------------------------------
Step = zeros(K, numel(ev));
for m = 1:numel(ev)
    Step(:, m) = double(t(:) >= ev(m));
end
eNames = "jump" + string(1:numel(ev));
nE = numel(ev);
% --- scale columns are per axis, so they are built inside the loop -----------
switch options.Scale
    case "none",  nS = 0; sNames = string.empty(1, 0);
    case "diag",  nS = 1; sNames = "scale_ii";
    otherwise,    nS = 3; sNames = "scale_i" + string(1:3);   % generic: row i of S, per axis
end
P = nB + nE + nS;

cal = struct();
cal.biasCoefficients = zeros(3, nB);
cal.jump = zeros(3, nE);
cal.scale = eye(3);
cal.sigma0 = zeros(3, 1);
cal.sigma = zeros(3, P);
cal.covariance = cell(1, 3);
cal.residual = nan(3, K);
cal.dof = zeros(3, 1);
cal.condition = zeros(3, 1);
cal.iterations = ones(3, 1);
cal.used = mask;
axisNames = ["x", "y", "z"];
for i = 1:3
    switch options.Scale
        case "none",  Ai = [Phi, Step];
        case "diag",  Ai = [Phi, Step, aModel(i, :)'];
        otherwise,    Ai = [Phi, Step, aModel'];
    end
    use = mask(i, :)' & W(i, :)' > 0;
    if nnz(use) < P
        NonGravForces.fail(unit, "numericSingular", "axis %s has %d usable samples for %d unknowns -> shorten the bias basis or widen the mask", axisNames(i), nnz(use), P);
    end
    A = Ai(use, :);
    y = yMeas(i, use)';
    w = W(i, use)';
    [x, s0, Cxx, cnd, iter, res] = localSolve(A, y, w, options, unit, axisNames(i));
    cal.biasCoefficients(i, :) = x(1:nB)';
    if nE > 0, cal.jump(i, :) = x(nB + (1:nE))'; end
    switch options.Scale
        case "diag", cal.scale(i, i) = x(end);
        case "full", cal.scale(i, :) = x(end-2:end)';
    end
    cal.sigma0(i) = s0;
    cal.sigma(i, :) = sqrt(max(diag(Cxx), 0))';
    cal.covariance{i} = Cxx;
    cal.residual(i, use) = res';
    cal.dof(i) = nnz(use) - P;
    cal.condition(i) = cnd;
    cal.iterations(i) = iter;
    cal.used(i, :) = use';
end
cal.bias = (Phi * cal.biasCoefficients')' + (Step * cal.jump')';
cal.names = [bNames, eNames, sNames];
info = struct('K', K, 'nBasis', nB, 'nScale', nS, 'nEvents', nE, 'basis', options.Bias, ...
    'provenance', NonGravForces.provenance("app.calibrate", 'Options', ...
        struct('Bias', options.Bias, 'Scale', options.Scale, 'KnotInterval_s', options.KnotInterval_s, ...
               'Robust', options.Robust, 'nEvents', nE)));
end

function [x, sigma0, Cxx, cnd, iter, res] = localSolve(A, y, w, options, unit, axisName)
%LOCALSOLVE Weighted least squares by QR on the equilibrated design matrix, optionally IRLS.
%
% NEVER THE NORMAL EQUATIONS. A'A squares the condition number, and this problem
% starts at 1e7 from the units alone; QR on A keeps it where it is.
iter = 1;
wk = w;
for it = 1:max(1, round(options.MaxIterations))
    sw = sqrt(wk);
    Aw = A .* sw;
    yw = y .* sw;
    % Column equilibration: each column divided by its norm, so the solve sees a
    % matrix whose conditioning is geometry rather than units. Undone below.
    cn = vecnorm(Aw, 2, 1);
    if any(cn == 0)
        NonGravForces.fail(unit, "numericSingular", "axis %s has an all-zero column in the design matrix -> a basis function that no usable sample touches cannot be estimated; coarsen KnotInterval_s", axisName);
    end
    As = Aw ./ cn;
    cnd = cond(As);
    if ~isfinite(cnd) || cnd > options.MaxCondition
        NonGravForces.fail(unit, "numericSingular", "axis %s: the equilibrated design matrix has condition %.3g above 'MaxCondition' %.3g -> the bias basis and the scale are not separable on this arc; coarsen the basis, shorten the arc, or fix the scale", axisName, cnd, options.MaxCondition);
    end
    xs = As \ yw;
    x = xs ./ cn(:);
    res = y - A * x;
    if options.Robust == "none"
        break
    end
    % Huber IRLS on the robust scale. The median absolute deviation is used rather
    % than the standard deviation, or the outliers this is meant to survive would
    % set the threshold that is supposed to catch them.
    r = res .* sqrt(w);
    s = 1.4826 * median(abs(r - median(r)));
    if s <= 0
        break
    end
    u = abs(r) / s;
    wNew = w .* min(1, options.HuberK ./ max(u, realmin));
    iter = it;
    if it > 1 && max(abs(wNew - wk)) <= 1e-10 * max(w)
        break
    end
    wk = wNew;
end
dof = numel(y) - size(A, 2);
vtv = sum(w .* res.^2);
sigma0 = sqrt(vtv / max(dof, 1));
% Covariance of the SCALED solution, then transformed back by the same diagonal.
n = size(As, 2);
Rs = triu(qr(As, 0));
Rs = Rs(1:n, 1:n);
Ninv = Rs \ (Rs' \ eye(n));            % (R'R)^-1 = R^-1 R^-T, never inv(A'A)
Cxx = sigma0^2 * (Ninv ./ (cn(:) * cn(:)'));
end
