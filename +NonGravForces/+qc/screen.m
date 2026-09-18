function [mask, events, info] = screen(y_meas, y_model, epoch, options)
%SCREEN Quality-control orchestrator: residuals -> outliers -> impulses -> jumps -> mask and event table.
%
%   [mask, events, info] = NonGravForces.qc.screen(y_meas, y_model, t)
%   [mask, events, info] = NonGravForces.qc.screen(y_meas, y_model, t, 'Bias', b, 'Scale', S, 'Response', crn, 'Thrusters', thr, 'Priors', ecl)
%
%   1. r = qc.residuals(y_meas, y_model, 'Bias', 'Scale', 'Whitening', 'Window')
%   2. outliers: qc.hampel(r, 'Window', HampelWindow, 'Threshold', HampelThreshold)
%   3. impulses (skipped, and said so in info, when no 'Response' is given: no
%      pulse shape is assumed): qc.detectImpulses runs on the residual with
%      only the ISOLATED outliers (runs of at most 2 flagged samples) replaced
%      by the local median. Hampel flags the peaks of a strong pulse as well,
%      in runs of several samples; replacing those destroys the pulse
%      (measured 2026-09-07: a 6-sigma pulse turned into four level jumps),
%      while a spike of one or two samples must not reach the matched filter
%      (it scores up to 0.26 x its height in sigma). Detected pulses claim
%      their samples (the outlier flags inside them are withdrawn) and are
%      blanked before the jump search; flagged runs that no pulse claims stay
%      outliers.
%   4. jumps: qc.detectJumps with 'JumpModel', 'JumpPenalty', 'JumpMinSegment'
%      and the 'Priors' (eclipse, MAS1B) on the residual with outliers
%      replaced by the local median and pulses blanked
%   mask marks outliers and the samples spanned by impulses; jumps are
%   reported, never masked (they become calibration columns). events is the
%   concatenated qc.eventTable of the three detectors, sorted by tStart.
%
% INPUT ARGUMENTS
%   Name      Size  Type      Units  Description
%   ----      ----  ----      -----  -----------
%   y_meas    KxM   double    units  Measured observable, M = 1 or 3 channels (columns)
%   y_model   KxM   double    units  Modelled observable
%   epoch     Kx1   datetime  UTC    Sample epochs
%
% NAME-VALUE ARGUMENTS
%   Name              Size     Type    Units    Description
%   ----              ----     ----    -----    -----------
%   Bias              1xM|KxM  double  units    Preliminary bias, default zeros (as qc.residuals)
%   Scale             1xM|MxM  double  -        Preliminary scale, diagonal as 1xM or full MxM,
%                                               default identity (as qc.residuals)
%   Whitening         1x1      string  -        "none" (default) | "difference" | "highpass"
%   Window            1x1      double  samples  Odd length >= 3, required by "highpass"
%   HampelWindow      1x1      double  samples  default 61
%   HampelThreshold   1x1      double  sigma    default 4
%   Response          Lx1      double  any      Pulse template for detectImpulses; default [] = impulse step skipped
%   Thrusters         Nx8      table   -        Impulse priors (io.readTHR1B), default none
%   ImpulseThreshold  1x1      double  sigma    default 5
%   Priors            Nx8      table   -        Jump priors (eclipse, MAS1B events), default none
%   JumpModel         1x1      string  -        "level" (default) | "slope"
%   JumpPenalty       1x1      double  units^2  default BIC (detectJumps)
%   JumpMinSegment    1x1      double  samples  default 30 (the CRN width: shorter level changes are
%                                              pulses or outliers, not jumps; the BIC alone admits
%                                              3-sigma noise segments of a few samples)
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   mask      KxM   logical  -      true = do not use the sample (outlier or inside an impulse)
%   events    Nx8   table    -      qc.eventTable of all detections
%   info      1x1   struct   -      .residual (KxM, whitened, untouched), .residualClean (outliers replaced,
%                                   pulses NaN), .outliers (KxM logical, after pulses claimed theirs),
%                                   .impulses .jumps (detector infos, or "skipped"), .isolatedOutliers
%                                   (KxM logical, the runs of <= 2 samples replaced before the matched filter),
%                                   .counts (struct outlier/impulse/jump), .provenance
%
% ACCURACY
%   Composition of the four units; see their ACCURACY blocks. The order is
%   fixed and the reason for each step is stated above. A level step may
%   register once as an impulse candidate at its edge (the zero-mean template
%   responds to the edge); the jump is then located at the end of that window,
%   up to L samples late. Measured on the synthetic day (2 pulses, 1 step,
%   1 spike in 3000 white samples): spike masked, both pulses masked over
%   their width, the step reported once (6 samples late behind the edge
%   candidate), no other jump.
%
% ERRORS
%   NonGravForces:screen:inputShape     y_meas/y_model not KxM with M = 1 or 3, or epoch not Kx1 datetime
%   NonGravForces:screen:inputEpoch     epoch not UTC
%   NonGravForces:residuals:*, NonGravForces:hampel:*, NonGravForces:detectImpulses:*, NonGravForces:detectJumps:*  propagated
%
% REFERENCES
%   Theory §10; HANDOVER D.3 workflow and Part J item 9.
%
% EXAMPLE
%   [ep, a] = NonGravForces.io.readACC1B(file); t = ...;   % Kx1 datetime UTC
%   [mask, ev] = NonGravForces.qc.screen(a', aModel', t, 'Whitening', "highpass", 'Window', 601, 'Response', crn, 'Thrusters', thr);
%
% See also NonGravForces.qc.residuals, NonGravForces.qc.hampel, NonGravForces.qc.detectImpulses, NonGravForces.qc.detectJumps, NonGravForces.qc.eventTable.
%
% KEYWORDS screen quality control mask events orchestrator
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    y_meas double
    y_model double
    epoch
    options.Bias double = []
    options.Scale double = []
    options.Whitening (1,1) string = "none"
    options.Window double = []
    options.HampelWindow (1,1) double = 61
    options.HampelThreshold (1,1) double = 4
    options.Response double = []
    options.Thrusters = []
    options.ImpulseThreshold (1,1) double = 5
    options.Priors = []
    options.JumpModel (1,1) string = "level"
    options.JumpPenalty double = []
    options.JumpMinSegment (1,1) double = 30
end
unit = "screen";
[K, M] = size(y_meas);
if ~ismatrix(y_meas) || ~isequal(size(y_meas), size(y_model)) || ~(M == 1 || M == 3)
    NonGravForces.fail(unit, "inputShape", "y_meas and y_model must be Kx1 or Kx3 of the same size -> got %s and %s", mat2str(size(y_meas)), mat2str(size(y_model)));
end
if ~isdatetime(epoch) || ~isequal(size(epoch), [K 1])
    NonGravForces.fail(unit, "inputShape", "epoch must be a Kx1 datetime -> got %s", class(epoch));
end
if isempty(epoch.TimeZone) || ~strcmp(epoch.TimeZone, 'UTC')
    NonGravForces.fail(unit, "inputEpoch", "epoch must carry TimeZone UTC -> set epoch.TimeZone = 'UTC'");
end
[r, rinfo] = NonGravForces.qc.residuals(y_meas, y_model, 'Bias', options.Bias, 'Scale', options.Scale, 'Whitening', options.Whitening, 'Window', options.Window);
[outl, hinfo] = NonGravForces.qc.hampel(r, 'Window', options.HampelWindow, 'Threshold', options.HampelThreshold);
% isolated outliers: flagged runs of at most 2 samples (a spike); longer runs may be the peaks of a pulse
isolated = false(K, M);
for m = 1:M
    d = diff([false; outl(:, m); false]); starts = find(d == 1); ends = find(d == -1) - 1;
    for q = reshape(find(ends - starts + 1 <= 2), 1, [])
        isolated(starts(q):ends(q), m) = true;
    end
end
rc = r; rc(outl) = hinfo.median(outl);              % all outliers replaced: the jump search sees no spikes and no pulse peaks
ri = r; ri(isolated) = hinfo.median(isolated);      % only spikes replaced: the matched filter sees the pulses intact
axes = "all"; if M == 3, axes = ["x","y","z"]; end
pulseMask = false(K, M);
info.isolatedOutliers = isolated;
if isempty(options.Response)
    evI = NonGravForces.qc.eventTable(); info.impulses = "skipped: no Response given";
else
    L = numel(options.Response);
    [evI, iinfo] = NonGravForces.qc.detectImpulses(ri, 'Epoch', epoch, 'Response', options.Response, 'Priors', options.Thrusters, 'Threshold', options.ImpulseThreshold);
    info.impulses = iinfo;
    for e = 1:height(evI)
        q = find(epoch == evI.tStart(e), 1);
        cols = 1:M; if M == 3, cols = find(axes == evI.axis(e)); end
        span = q:min(K, q + L - 1);
        pulseMask(span, cols) = true;
        outl(span, cols) = false;                     % the pulse owns its samples
        rc(span, cols) = r(span, cols);               % restore, then blank for the jump search
        rc(span, cols) = NaN;
    end
end
mask = outl | pulseMask;
% outlier events: one row per flagged sample
[iK, iM] = find(outl);
tO = epoch(iK); amO = r(outl) - hinfo.median(outl); sgO = hinfo.sigma(outl);
evO = NonGravForces.qc.eventTable(tO, tO, repmat("outlier", numel(iK), 1), axes(iM)', amO, sgO, repmat("unassigned", numel(iK), 1), min(1, abs(amO) ./ max(sgO, realmin) / options.HampelThreshold));
jargs = {'Epoch', epoch, 'Model', options.JumpModel, 'Priors', options.Priors, 'MinSegment', options.JumpMinSegment};
if ~isempty(options.JumpPenalty), jargs = [jargs, {'Penalty', options.JumpPenalty}]; end
[evJ, jinfo] = NonGravForces.qc.detectJumps(rc, jargs{:});
info.jumps = jinfo;
events = sortrows([evO; evI; evJ], 'tStart');
info.residual = r;
info.residualClean = rc;
info.residualInfo = rinfo;
info.outliers = outl;
info.counts = struct('outlier', height(evO), 'impulse', height(evI), 'jump', height(evJ));
info.provenance = NonGravForces.provenance("qc.screen", 'Options', struct('Whitening', options.Whitening, 'Window', options.Window, 'HampelWindow', options.HampelWindow, 'HampelThreshold', options.HampelThreshold, 'ImpulseThreshold', options.ImpulseThreshold, 'JumpModel', options.JumpModel, 'JumpMinSegment', options.JumpMinSegment, 'responseGiven', ~isempty(options.Response)));
end
