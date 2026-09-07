function [r, info] = residuals(y_meas, y_model, options)
%RESIDUALS Residuals of measured minus modelled observations under a preliminary calibration, optionally whitened.
%
%   r = NonGravForces.qc.residuals(y_meas, y_model)
%   [r, info] = NonGravForces.qc.residuals(y_meas, y_model, 'Bias', b, 'Scale', S, 'Whitening', "highpass", 'Window', 601)
%
%   r_raw = y_meas - (b + y_model * S')      per channel (columns), b 1xM, S MxM
%
%   with the preliminary calibration of the caller (defaults: b = 0, S = I,
%   i.e. the plain difference). Whitening prepares the residual for the
%   detectors, which assume white noise:
%     "none"       r = r_raw
%     "difference" r = diff(r_raw) padded with NaN in the first row: flattens
%                  the 1/f accelerometer noise, doubles the white-noise variance,
%                  turns a level jump into an impulse of one sample
%     "highpass"   r = r_raw - movmean(r_raw, Window, 'omitnan'): removes bias
%                  drift slower than Window samples and keeps jumps and pulses
%                  shorter than Window/2 intact (the choice for detectJumps)
%   NaN samples propagate (they never enter the moving mean) and are never
%   invented; the detectors handle them.
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   y_meas    KxM   double  any    K samples of M channels (columns), e.g. ACT1B accelerations
%   y_model   KxM   double  any    Modelled observable at the same samples (sensor.observe)
%
% NAME-VALUE ARGUMENTS
%   Name       Size     Type    Units    Description
%   ----       ----     ----    -----    -----------
%   Bias       1xM|KxM  double  units    Preliminary bias, default zeros
%   Scale      1xM|MxM  double  -        Preliminary scale (diagonal as 1xM, or full MxM), default identity
%   Whitening  1x1      string  -        "none" (default) | "difference" | "highpass"
%   Window     1x1      double  samples  Odd length >= 3 for "highpass" (required by it)
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   r         KxM   double  units  Whitened residuals
%   info      1x1   struct  -      .raw (KxM, before whitening), .sigmaRobust (1xM, 1.4826 MAD of r,
%                                  NaN-free), .whitening, .window, .provenance
%
% ACCURACY
%   Arithmetic only. "highpass" attenuates a step of height h to h (1 - 1/2)
%   at its centre over Window/2 samples on either side; "difference" is exact
%   for a step (one sample of height h). The robust sigma is consistent with
%   the standard deviation for Gaussian noise; for 1/f noise it overestimates
%   by up to 10 % inside short windows (conservative for the detectors).
%
% ERRORS
%   NonGravForces:residuals:inputShape      y_meas and y_model not KxM of the same size, Bias/Scale of the wrong size
%   NonGravForces:residuals:inputWindow     highpass without an odd Window >= 3
%
% REFERENCES
%   Theory §10; Klinger B., Mayer-Gürr T. (2016) Adv. Space Res. 58:1597 (ACC
%   calibration model: bias, scale); HANDOVER Part J item 9 (priors, not rediscovery).
%
% EXAMPLE
%   [ep, a] = NonGravForces.io.readACC1B(file);
%   r = NonGravForces.qc.residuals(a', aModel', 'Scale', [1.05 1 1], 'Whitening', "highpass", 'Window', 601);
%
% See also NonGravForces.qc.screen, NonGravForces.qc.hampel, NonGravForces.qc.detectImpulses, NonGravForces.qc.detectJumps.
%
% KEYWORDS residuals whitening preliminary calibration quality control
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    y_meas double
    y_model double
    options.Bias double = []
    options.Scale double = []
    options.Whitening (1,1) string {mustBeMember(options.Whitening, ["none","difference","highpass"])} = "none"
    options.Window double = []
end
unit = "residuals";
if ~ismatrix(y_meas) || ~isequal(size(y_meas), size(y_model))
    NonGravForces.fail(unit, "inputShape", "y_meas and y_model must be KxM of the same size -> got %s and %s (samples down the columns)", mat2str(size(y_meas)), mat2str(size(y_model)));
end
[K, M] = size(y_meas);
b = options.Bias;
if isempty(b), b = zeros(1, M); end
if ~(isequal(size(b), [1 M]) || isequal(size(b), [K M]))
    NonGravForces.fail(unit, "inputShape", "Bias must be 1xM or KxM with M = %d -> got %s", M, mat2str(size(b)));
end
S = options.Scale;
if isempty(S), S = eye(M); elseif isequal(size(S), [1 M]), S = diag(S); end
if ~isequal(size(S), [M M])
    NonGravForces.fail(unit, "inputShape", "Scale must be 1xM or MxM with M = %d -> got %s", M, mat2str(size(options.Scale)));
end
raw = y_meas - (b + y_model * S');
switch options.Whitening
    case "none"
        r = raw;
    case "difference"
        r = [nan(1, M); diff(raw, 1, 1)];
    case "highpass"
        W = options.Window;
        if isempty(W) || ~isscalar(W) || W < 3 || mod(W, 2) ~= 1 || W ~= fix(W)
            NonGravForces.fail(unit, "inputWindow", "highpass needs an odd Window >= 3 samples -> pass 'Window', 601 for a 10-min drift scale at 1 Hz");
        end
        r = raw - movmean(raw, W, 1, 'omitnan');
end
sig = nan(1, M);
for m = 1:M
    x = r(~isnan(r(:, m)), m);
    if numel(x) >= 3, sig(m) = 1.4826 * median(abs(x - median(x))); end
end
info.raw = raw;
info.sigmaRobust = sig;
info.whitening = options.Whitening;
info.window = options.Window;
info.provenance = NonGravForces.provenance("qc.residuals", 'Options', struct('Whitening', options.Whitening, 'Window', options.Window, 'biasGiven', ~isempty(options.Bias), 'scaleGiven', ~isempty(options.Scale)));
end
