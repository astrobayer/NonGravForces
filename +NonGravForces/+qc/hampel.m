function [mask, info] = hampel(x, options)
%HAMPEL Robust point-outlier detector: moving median and scaled MAD, per column.
%
%   mask = NonGravForces.qc.hampel(x)
%   [mask, info] = NonGravForces.qc.hampel(x, 'Window', 61, 'Threshold', 4)
%
%   Flags sample i of column j when |x(i,j) - med(i,j)| > t * 1.4826 * MAD(i,j),
%   where med and MAD are computed over a centred window of W samples
%   (truncated at the edges). Median and MAD are unaffected by up to ~50 %
%   contaminated samples in a window, which is why this - and not a moving
%   mean/std - is the first pass on accelerometer residuals containing
%   thruster spikes and twangs (see the ACCURACY block).
%
%   NaN POLICY. NaN samples are never flagged and never enter the window
%   statistics; a window with fewer than 3 finite samples flags nothing.
%
% INPUT ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   x          KxM   double  any    K samples of M channels (columns). A row
%                                   vector is treated as one channel.
%
% NAME-VALUE ARGUMENTS
%   Window     1x1   double  samples Odd window length >= 3, default 61
%                                   (1 min at 1 Hz: longer than a CRN-filtered
%                                   impulse, shorter than an orbit).
%   Threshold  1x1   double  -      Multiple of the robust sigma, default 4
%                                   (Gaussian false-alarm rate 6e-5/sample).
%   MinSigma   1x1   double  units  Floor on the robust sigma, default 0:
%                                   set to the instrument noise level so that
%                                   a flat window cannot flag quantisation.
%
% OUTPUT ARGUMENTS
%   Name       Size  Type     Units  Description
%   ----       ----  ----     -----  -----------
%   mask       KxM   logical  -      true = outlier
%   info       1x1   struct   -      .median (KxM) .sigma (KxM, robust)
%                                    .score (KxM, |x-med|/sigma) .nFlagged
%                                    (1xM) .provenance
%
% ACCURACY
%   The scale factor 1.4826 makes MAD consistent with sigma for Gaussian
%   noise; for the coloured 1/f accelerometer noise it overestimates sigma
%   by up to ~10 % inside a 61-sample window, which is conservative.
%   Impulses wider than Window/2 are NOT outliers to this detector by design
%   - they go to detectImpulses.
%
% ERRORS
%   NonGravForces:hampel:inputShape    x is not a 2-D real double array
%   NonGravForces:hampel:inputWindow   Window even, < 3, or non-integer
%
% REFERENCES
%   Hampel (1974) J. Am. Stat. Assoc. 69:383; Pearson (2002) EURASIP JASP.
%
% EXAMPLE
%   x = randn(1000,1); x(500) = 12;
%   mask = NonGravForces.qc.hampel(x); find(mask)   % -> 500
%
% See also NonGravForces.qc.screen, NonGravForces.qc.detectImpulses, movmedian.
%
% KEYWORDS outlier robust median MAD quality control
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    x double
    options.Window (1,1) double = 61
    options.Threshold (1,1) double {mustBePositive} = 4
    options.MinSigma (1,1) double {mustBeNonnegative} = 0
end
if ~ismatrix(x) || ~isreal(x)
    NonGravForces.fail("hampel", "inputShape", "x must be a real KxM matrix -> got %s", mat2str(size(x)));
end
W = options.Window;
if W < 3 || mod(W, 2) ~= 1 || W ~= round(W)
    NonGravForces.fail("hampel", "inputWindow", "Window must be an odd integer >= 3 -> got %g", W);
end
if isrow(x), x = x(:); end
[K, M] = size(x);
h = (W - 1) / 2;
med = nan(K, M); sig = nan(K, M);
for j = 1:M
    xj = x(:, j);
    for i = 1:K
        lo = max(1, i - h); hi = min(K, i + h);
        w = xj(lo:hi); w = w(isfinite(w));
        if numel(w) < 3, continue; end
        m = median(w);
        med(i, j) = m;
        sig(i, j) = 1.4826 * median(abs(w - m));
    end
end
sig = max(sig, options.MinSigma);
score = abs(x - med) ./ sig;                     % NaN where sig is 0 or x NaN
mask = score > options.Threshold;                % NaN > t is false
mask(~isfinite(x)) = false;
info.median = med; info.sigma = sig; info.score = score;
info.nFlagged = sum(mask, 1);
info.provenance = NonGravForces.provenance("qc.hampel", 'Options', options);
end
