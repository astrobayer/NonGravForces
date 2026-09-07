function [events, info] = detectImpulses(r, options)
%DETECTIMPULSES Matched-filter detection of pulses with a known shape (the CRN impulse response) in whitened residuals; thruster priors are assigned first.
%
%   [events, info] = NonGravForces.qc.detectImpulses(r, 'Epoch', t, 'Response', h)
%   [events, info] = NonGravForces.qc.detectImpulses(r, 'Epoch', t, 'Response', h, 'Priors', thr, 'Threshold', 5)
%
%   For each channel the residual is correlated with the ZERO-MEAN, unit-norm
%   template hz = (h - mean h)/|h - mean h| (L samples, the pulse starting at
%   its first sample):
%
%       c(t) = sum_i r(t+i-1) hz(i),   score(t) = c(t) / sigma_c,
%
%   sigma_c the robust (1.4826 MAD) scale of c. Removing the template mean makes
%   the filter blind to a constant offset or a level jump inside the window
%   (measured 2026-09-07: a 4-sigma step produced a plateau of false pulses with
%   the raw template); only the step's edge can still register, as one
%   candidate at most. The last L-1 samples are not evaluable (the window
%   would run past the end) and never host a detection. Local maxima of
%   |score| above Threshold, at least L
%   samples apart (the strongest wins), are pulses with amplitude a = c/(hz.h)
%   (the height multiplying h), sigma_a = sigma_c/(hz.h),
%   tStart at the peak sample and tEnd L-1 samples later. A pulse within
%   PriorWindow_s of a prior row (e.g. THR1B firings from io.readTHR1B) takes
%   the prior's source and confidence 1 (HANDOVER Part J item 9: priors are
%   assigned, not rediscovered); the others are "unassigned" (twang or impact
%   candidates) with confidence erf(|score|/sqrt(2)) - the two-sided Gaussian
%   probability that noise alone does not reach the score. NaN samples of r
%   are treated as zero in the correlation and counted in info.
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   r         KxM   double  units  Whitened residuals (qc.residuals), M = 1 (axis "all") or 3 (x, y, z)
%
% NAME-VALUE ARGUMENTS
%   Name           Size   Type      Units    Description
%   ----           ----   ----      -----    -----------
%   Epoch          Kx1    datetime  UTC      Sample epochs (required; eventTable times)
%   Response       Lx1    double    any      Pulse template, L >= 2, not constant (required: the CRN impulse response at the
%                                            sample rate; no default shape is assumed)
%   Priors         Nx8    table     -        qc.eventTable rows (class impulse) to assign, default none
%   Threshold      1x1    double    sigma    Detection threshold on |score|, default 5 (Gaussian false alarm
%                                            6e-7 per sample, ~0.05 per day at 1 Hz)
%   PriorWindow_s  1x1    double    s        Assignment window around a prior epoch, default 5
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   events    Nx8   table   -      qc.eventTable rows, class "impulse", sorted by tStart
%   info      1x1   struct  -      .score (KxM), .amplitude (KxM, c/(hz.h)), .sigmaCorr (1xM), .templateGain (hz.h),
%                                  .nanCount (1xM), .unmatchedPriors (rows of Priors without a detection),
%                                  .provenance
%
% ACCURACY
%   Epoch of a planted pulse recovered to +-1 sample and amplitude to 3 sigma_a
%   for pulses above Threshold (measured on synthetic CRN-shaped pulses in the
%   test suite). Two pulses closer than L samples merge into the stronger one
%   by design (impulses shorter than the CRN width are unresolvable in L1B,
%   theory §10). A template that does not match the true pulse shape lowers
%   the score by the cosine between the shapes, never raises it.
%
% ERRORS
%   NonGravForces:detectImpulses:inputShape      r not KxM with M = 1 or 3, or Epoch not Kx1 datetime
%   NonGravForces:detectImpulses:inputEpoch      Epoch missing or not UTC
%   NonGravForces:detectImpulses:inputResponse   Response missing, shorter than 2 samples, longer than K, or constant
%   NonGravForces:detectImpulses:inputPriors     Priors not an eventTable-shaped table
%   NonGravForces:detectImpulses:inputRange      Threshold or PriorWindow_s <= 0
%
% REFERENCES
%   Theory §10; Bandikova T. et al. (2019) Adv. Space Res. 64:623 (CRN impulse
%   responses per thruster pair); Meyer U. et al. (2012) GRACE ACC twangs;
%   Turin G.L. (1960) An introduction to matched filters, IRE Trans. IT-6.
%
% EXAMPLE
%   thr = NonGravForces.io.readTHR1B(thrFile);
%   [ev, info] = NonGravForces.qc.detectImpulses(r, 'Epoch', t, 'Response', crn, 'Priors', thr);
%   ev(ev.source == "unassigned", :)                 % twang / impact candidates
%
% See also NonGravForces.qc.residuals, NonGravForces.qc.eventTable, NonGravForces.io.readTHR1B, NonGravForces.qc.screen.
%
% KEYWORDS impulse matched filter thruster twang quality control
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    r double
    options.Epoch = []
    options.Response double = []
    options.Priors = []
    options.Threshold (1,1) double = 5
    options.PriorWindow_s (1,1) double = 5
end
unit = "detectImpulses";
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
h = options.Response(:);
if numel(h) < 2 || numel(h) > K || any(~isfinite(h)) || ~any(h - h(1))
    NonGravForces.fail(unit, "inputResponse", "'Response' must be a finite, non-constant template of 2..K samples -> pass the CRN impulse response at the sample rate");
end
if options.Threshold <= 0 || options.PriorWindow_s <= 0
    NonGravForces.fail(unit, "inputRange", "Threshold and PriorWindow_s must be > 0 -> got %g, %g", options.Threshold, options.PriorWindow_s);
end
priors = options.Priors;
if isempty(priors), priors = NonGravForces.qc.eventTable(); end
schema = NonGravForces.qc.eventTable().Properties.VariableNames;
if ~istable(priors) || ~isequal(priors.Properties.VariableNames, schema)
    NonGravForces.fail(unit, "inputPriors", "Priors must be a qc.eventTable -> use io.readTHR1B or qc.eventTable(...)");
end
L = numel(h);
hz = h - mean(h);
hn = hz / norm(hz);                                 % zero-mean, unit-norm matched filter
hnorm = dot(hn, h);                                 % gain: a pulse a*h gives c = a * hnorm
axes = "all"; if M == 3, axes = ["x","y","z"]; end
score = zeros(K, M); amp = zeros(K, M); sigC = nan(1, M); nanCount = zeros(1, M);
tS = datetime.empty(0, 1); tS.TimeZone = 'UTC'; tE = tS; ax = string.empty(0, 1); am = zeros(0, 1); sg = zeros(0, 1); sc = zeros(0, 1);
for m = 1:M
    x = r(:, m); bad = isnan(x); nanCount(m) = nnz(bad); x(bad) = 0;
    c = zeros(K, 1);
    for i = 1:L                                     % c(t) = sum_i x(t+i-1) hz(i), zero beyond the end
        c(1:K - i + 1) = c(1:K - i + 1) + hn(i) * x(i:K);
    end
    c(K - L + 2:K) = 0;                             % a window running past the end sees the zero padding as a step: not evaluable
    fin = c(~bad);
    s = 1.4826 * median(abs(fin - median(fin)));
    if ~(s > 0), s = NaN; end
    sigC(m) = s;
    score(:, m) = c / s;
    amp(:, m) = c / hnorm;
    % local maxima of |score| above the threshold, greedy suppression within L samples
    a = abs(score(:, m)); a(bad) = 0;
    cand = find(a >= options.Threshold);
    [~, order] = sort(a(cand), 'descend'); cand = cand(order);
    taken = false(K, 1);
    for q = cand'
        if taken(q), continue; end
        lo = max(1, q - L + 1); hi = min(K, q + L - 1);
        taken(lo:hi) = true;
        tS(end + 1, 1) = t(q); tE(end + 1, 1) = t(min(K, q + L - 1)); %#ok<AGROW>
        ax(end + 1, 1) = axes(m); am(end + 1, 1) = amp(q, m); sg(end + 1, 1) = s / hnorm; sc(end + 1, 1) = a(q); %#ok<AGROW>
    end
end
n = numel(tS);
source = repmat("unassigned", n, 1);
conf = erf(sc / sqrt(2));
matched = false(height(priors), 1);
for e = 1:n
    dt = abs(seconds(priors.tStart - tS(e)));
    [dmin, p] = min(dt);
    if ~isempty(dmin) && dmin <= options.PriorWindow_s
        source(e) = priors.source(p); conf(e) = 1; matched(p) = true;
    end
end
events = NonGravForces.qc.eventTable(tS, tE, repmat("impulse", n, 1), ax, am, sg, source, conf);
events = sortrows(events, 'tStart');
info.score = score;
info.amplitude = amp;
info.sigmaCorr = sigC;
info.templateGain = hnorm;
info.templateLength = L;
info.nanCount = nanCount;
info.unmatchedPriors = priors(~matched, :);
info.provenance = NonGravForces.provenance("qc.detectImpulses", 'Options', struct('Threshold', options.Threshold, 'PriorWindow_s', options.PriorWindow_s, 'L', L, 'priors', height(priors)));
end
