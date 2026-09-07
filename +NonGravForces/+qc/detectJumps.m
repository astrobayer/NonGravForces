function [events, info] = detectJumps(r, options)
%DETECTJUMPS Change points in level or slope of whitened residuals by PELT with a BIC penalty; eclipse and mass priors are assigned.
%
%   [events, info] = NonGravForces.qc.detectJumps(r, 'Epoch', t)
%   [events, info] = NonGravForces.qc.detectJumps(r, 'Epoch', t, 'Model', "slope", 'Priors', ecl, 'Penalty', beta)
%
%   Minimises, per channel, sum over segments of the least-squares cost plus a
%   penalty beta per change point (Killick, Fearnhead & Eckley 2012, Pruned
%   Exact Linear Time). The cost of a segment is its residual sum of squares
%   about a constant ("level") or a straight line ("slope"). The default
%   penalty is the BIC value (p + 1) sigma^2 log K with p = 1 (level) or 2
%   (slope) parameters per segment plus one for the location of the change
%   point (Killick et al. 2012 and the changepoint package convention; with
%   p alone, measured 2026-09-07, 1000 white samples with one 4-sigma step
%   gave four change points, three of them two-sample noise excursions) and
%   sigma the robust scale of the first differences (1.4826 MAD / sqrt(2)); an
%   explicit 'Penalty' replaces it (in squared residual units). PELT is exact: its segmentation equals the exhaustive
%   optimum, which the test suite checks against tools/pyxcheck/pelt_check.py.
%
%   Each change point becomes a "jump" row: tStart = tEnd = first sample of the
%   new segment, amplitude = level of the new segment minus level of the old
%   one at the boundary (for "slope" the fitted lines evaluated there),
%   sigma = sigma sqrt(1/n1 + 1/n2). A jump within PriorWindow_s of a prior
%   row (eclipse entry/exit, MAS1B tank event) takes that source and
%   confidence 1; the others are "unassigned" with confidence
%   erf(|amplitude|/sigma/sqrt(2)). Jumps are never cut from the data: they
%   become Heaviside columns in the calibration (theory §10).
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   r         KxM   double  units  Whitened residuals, M = 1 (axis "all") or 3 (x, y, z); NaN samples are
%                                  removed before segmentation and the indices mapped back
%
% NAME-VALUE ARGUMENTS
%   Name           Size   Type      Units    Description
%   ----           ----   ----      -----    -----------
%   Epoch          Kx1    datetime  UTC      Sample epochs (required)
%   Model          1x1    string    -        "level" (default) | "slope"
%   Penalty        1x1    double    units^2  Penalty per change point, default BIC (see above)
%   MinSegment     1x1    double    samples  Shortest segment, default 2 (level) or 3 (slope)
%   Priors         Nx8    table     -        qc.eventTable rows to assign, default none
%   PriorWindow_s  1x1    double    s        Assignment window, default 60
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   events    Nx8   table   -      qc.eventTable rows, class "jump", sorted by tStart
%   info      1x1   struct  -      .changePoints (1xM cell of sample indices: first sample of each new
%                                  segment), .penalty (1xM), .sigma (1xM), .cost (1xM, optimal objective),
%                                  .segments (1xM cell of [first last level|intercept slope]), .provenance
%
% ACCURACY
%   Exact optimisation for the stated cost and penalty (PELT); the statistical
%   claims are those of the BIC: a step of height h in white noise of scale
%   sigma is found when h^2 n1 n2/(n1+n2) > beta = 2 sigma^2 log K, i.e. above
%   ~1.7 sigma for a segment pair of 5 samples at K = 1000 (measured: planted steps of 4 sigma
%   recovered to +-1 sample, amplitude within 3 sigma). Coloured 1/f noise
%   without highpass whitening produces spurious level changes: whiten first
%   (qc.residuals "highpass" or "difference").
%
% ERRORS
%   NonGravForces:detectJumps:inputShape     r not KxM with M = 1 or 3, or Epoch not Kx1
%   NonGravForces:detectJumps:inputEpoch     Epoch missing or not UTC
%   NonGravForces:detectJumps:inputRange     Penalty <= 0, MinSegment < 2 (or < 3 for slope), PriorWindow_s <= 0
%   NonGravForces:detectJumps:inputPriors    Priors not an eventTable-shaped table
%   NonGravForces:detectJumps:dataEmpty      fewer than 2 MinSegment finite samples in a channel
%
% REFERENCES
%   Killick R., Fearnhead P., Eckley I.A. (2012) Optimal detection of changepoints
%   with a linear computational cost, JASA 107:1590-1598; Schwarz G. (1978) Ann.
%   Stat. 6:461 (BIC); theory §10.
%
% EXAMPLE
%   r = NonGravForces.qc.residuals(a, aModel, 'Whitening', "highpass", 'Window', 601);
%   [ev, info] = NonGravForces.qc.detectJumps(r, 'Epoch', t, 'Priors', eclipseEvents);
%
% See also NonGravForces.qc.residuals, NonGravForces.qc.eventTable, NonGravForces.qc.screen.
%
% KEYWORDS jump change point PELT BIC quality control
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    r double
    options.Epoch = []
    options.Model (1,1) string {mustBeMember(options.Model, ["level","slope"])} = "level"
    options.Penalty double = []
    options.MinSegment double = []
    options.Priors = []
    options.PriorWindow_s (1,1) double = 60
end
unit = "detectJumps";
[K, M] = size(r);
if ~ismatrix(r) || ~(M == 1 || M == 3)
    NonGravForces.fail(unit, "inputShape", "r must be Kx1 or Kx3 -> got %s (samples down the columns)", mat2str(size(r)));
end
t = options.Epoch;
if isempty(t) || ~isdatetime(t)
    NonGravForces.fail(unit, "inputEpoch", "'Epoch' (Kx1 datetime, UTC) is required -> pass the sample epochs of r");
end
if ~isequal(size(t), [K 1])
    NonGravForces.fail(unit, "inputShape", "Epoch must be Kx1 = %dx1 -> got %s", K, mat2str(size(t)));
end
if isempty(t.TimeZone) || ~strcmp(t.TimeZone, 'UTC')
    NonGravForces.fail(unit, "inputEpoch", "Epoch must carry TimeZone UTC -> set t.TimeZone = 'UTC'");
end
slope = options.Model == "slope";
p = 1 + slope;
minSeg = options.MinSegment;
if isempty(minSeg), minSeg = 2 + slope; end
if ~isscalar(minSeg) || minSeg < 2 + slope || minSeg ~= fix(minSeg)
    NonGravForces.fail(unit, "inputRange", "MinSegment must be an integer >= %d for %s -> got %s", 2 + slope, options.Model, mat2str(minSeg));
end
if ~isempty(options.Penalty) && ~(isscalar(options.Penalty) && options.Penalty > 0)
    NonGravForces.fail(unit, "inputRange", "Penalty must be a positive scalar -> got %s", mat2str(options.Penalty));
end
if options.PriorWindow_s <= 0
    NonGravForces.fail(unit, "inputRange", "PriorWindow_s must be > 0 -> got %g", options.PriorWindow_s);
end
priors = options.Priors;
if isempty(priors), priors = NonGravForces.qc.eventTable(); end
schema = NonGravForces.qc.eventTable().Properties.VariableNames;
if ~istable(priors) || ~isequal(priors.Properties.VariableNames, schema)
    NonGravForces.fail(unit, "inputPriors", "Priors must be a qc.eventTable -> use qc.eventTable(...)");
end
axes = "all"; if M == 3, axes = ["x","y","z"]; end
cps = cell(1, M); segs = cell(1, M); pen = nan(1, M); sig = nan(1, M); cost = nan(1, M);
tS = datetime.empty(0, 1); tS.TimeZone = 'UTC'; ax = string.empty(0, 1); am = zeros(0, 1); sg = zeros(0, 1);
for m = 1:M
    idx = find(~isnan(r(:, m)));
    x = r(idx, m); n = numel(x);
    if n < 2 * minSeg
        NonGravForces.fail(unit, "dataEmpty", "channel %d has %d finite samples, fewer than 2 MinSegment = %d -> pass a longer series", m, n, 2 * minSeg);
    end
    d = diff(x); s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
    sig(m) = s;
    if isempty(options.Penalty)
        beta = (p + 1) * max(s, realmin)^2 * log(n);    % BIC: p segment parameters plus the change-point location
    else
        beta = options.Penalty;
    end
    pen(m) = beta;
    [cp, Fopt, seg] = localPelt(x, beta, minSeg, slope);
    cost(m) = Fopt;
    cps{m} = idx(cp)';                              % first sample (original index) of each new segment
    segs{m} = [idx(seg(:, 1)) idx(seg(:, 2)) seg(:, 3:end)];
    for q = 1:numel(cp)
        a = seg(q, :); b = seg(q + 1, :);           % old and new segment [first last intercept (slope)]
        n1 = a(2) - a(1) + 1; n2 = b(2) - b(1) + 1;
        if slope
            tb = b(1);                              % boundary sample (local index)
            levelOld = a(3) + a(4) * tb; levelNew = b(3) + b(4) * tb;
        else
            levelOld = a(3); levelNew = b(3);
        end
        tS(end + 1, 1) = t(idx(cp(q))); ax(end + 1, 1) = axes(m); %#ok<AGROW>
        am(end + 1, 1) = levelNew - levelOld; sg(end + 1, 1) = s * sqrt(1 / n1 + 1 / n2); %#ok<AGROW>
    end
end
n = numel(tS);
source = repmat("unassigned", n, 1);
conf = erf(abs(am) ./ max(sg, realmin) / sqrt(2));
for e = 1:n
    dt = abs(seconds(priors.tStart - tS(e)));
    [dmin, q] = min(dt);
    if ~isempty(dmin) && dmin <= options.PriorWindow_s
        source(e) = priors.source(q); conf(e) = 1;
    end
end
events = NonGravForces.qc.eventTable(tS, tS, repmat("jump", n, 1), ax, am, sg, source, conf);
events = sortrows(events, 'tStart');
info.changePoints = cps;
info.penalty = pen;
info.sigma = sig;
info.cost = cost;
info.segments = segs;
info.model = options.Model;
info.provenance = NonGravForces.provenance("qc.detectJumps", 'Options', struct('Model', options.Model, 'Penalty', options.Penalty, 'MinSegment', minSeg, 'PriorWindow_s', options.PriorWindow_s, 'priors', height(priors)));
end

function [cp, Fopt, seg] = localPelt(x, beta, minSeg, slope)
%LOCALPELT Exact penalised segmentation (Killick 2012): F(t) = min_tau F(tau) + C(tau+1..t) + beta.
n = numel(x);
tt = (1:n)';
S1 = [0; cumsum(x)]; S2 = [0; cumsum(x.^2)];
if slope
    St = [0; cumsum(tt)]; Stt = [0; cumsum(tt.^2)]; Sxt = [0; cumsum(x .* tt)];
end
F = [-beta; inf(n, 1)];                              % F(tau) stored at index tau+1
prev = zeros(n + 1, 1);
R = 0;                                               % candidate last change points tau (0-based), row
for t = minSeg:n
    adm = t - R >= minSeg;                           % admissible: the new segment tau+1..t has at least minSeg samples
    cand = reshape(R(adm), [], 1);                   % admissible candidates as a column
    if isempty(cand), continue; end
    a = cand + 1; nn = t - cand;                     % segment starts (1-based) and lengths, columns
    sx = S1(t + 1) - S1(a); sxx = S2(t + 1) - S2(a);
    if slope
        st = St(t + 1) - St(a); stt = Stt(t + 1) - Stt(a); sxt = Sxt(t + 1) - Sxt(a);
        det = nn .* stt - st.^2;
        bb = (nn .* sxt - st .* sx) ./ det;          % slope
        aa = (sx - bb .* st) ./ nn;                  % intercept
        c = sxx - 2 * aa .* sx - 2 * bb .* sxt + nn .* aa.^2 + 2 * aa .* bb .* st + bb.^2 .* stt;
    else
        c = sxx - sx.^2 ./ nn;
    end
    c = max(c, 0);
    vals = F(a) + c + beta;                          % F(tau) = F(a), one value per candidate (column)
    [Fmin, k] = min(vals);
    F(t + 1) = Fmin;
    prev(t + 1) = cand(k);
    R = [reshape(cand(vals - beta <= Fmin), 1, []), R(~adm), t];   % prune admissible tau with F(tau) + C > F(t); keep the too-recent ones
end
% backtrack
bounds = n; tau = prev(n + 1);
while tau > 0
    bounds(end + 1) = tau; %#ok<AGROW>
    tau = prev(tau + 1);
end
bounds = flip(bounds);                               % last sample of each segment
first = [1, bounds(1:end - 1) + 1];
cp = first(2:end);
Fopt = F(n + 1);
nseg = numel(bounds);
seg = zeros(nseg, 3 + slope);
for q = 1:nseg
    a = first(q); b = bounds(q); xs = x(a:b);
    if slope
        A = [ones(b - a + 1, 1) (a:b)']; coef = A \ xs; seg(q, :) = [a b coef'];
    else
        seg(q, :) = [a b mean(xs)];
    end
end
end
