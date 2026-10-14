function [events, info] = detectImpulses(r, varargin)
%DETECTIMPULSES [DEPRECATED, forwards to timeSeriesPreparation.qc.detectImpulses] Matched-filter detection of pulses with a known shape (the CRN impulse response) in whitened residuals; thruster priors are assigned first.
%
%   [events, info] = NonGravForces.qc.detectImpulses(r, 'Epoch', t, 'Response', h)
%   [events, info] = NonGravForces.qc.detectImpulses(r, 'Epoch', t, 'Response', h, 'Priors', thr, 'Threshold', 5)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05): it ported this unit UNCHANGED, froze the port
%   and gates every later change against it with isequaln (X-16). This unit
%   now only forwards: it warns NonGravForces:detectImpulses:deprecated and calls
%   timeSeriesPreparation.qc.detectImpulses with the same arguments. It is deleted in
%   the release after the one that carries the forwarding. The target's errors
%   propagate unwrapped (timeSeriesPreparation:detectImpulses:*), as every foreign error
%   does in this toolbox. Call the target directly.
%
% INPUT ARGUMENTS
%   Name      Size     Type     Units Description
%   ----      ----     ----     ----- -----------
%   r         Kx1|Kx3  double   any   residual series, as timeSeriesPreparation.qc.detectImpulses
%   varargin  1xN      cell     -     its name-value arguments, passed through unchanged
%
% NAME-VALUE ARGUMENTS
%   Epoch, Response, Priors, Threshold, PriorWindow_s: those of timeSeriesPreparation.qc.detectImpulses,
%   with its defaults; the two argument blocks are identical (compared 2026-10-14).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type  Units Description
%   ----      ----  ----  ----- -----------
%   events    Nx8   table   -  as timeSeriesPreparation.qc.detectImpulses, unchanged
%   info      1x1   struct  -  as timeSeriesPreparation.qc.detectImpulses, unchanged
%
% ACCURACY
%   None of its own: the target's numbers, bit for bit. tests/tNgfQc*.m run the
%   pre-move suite through this forwarder.
%
% ERRORS
%   NonGravForces:detectImpulses:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:detectImpulses:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:detectImpulses:deprecated   on every call; call timeSeriesPreparation.qc.detectImpulses
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05, X-16.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   help timeSeriesPreparation.qc.detectImpulses
%
% See also timeSeriesPreparation.qc.detectImpulses.
%
% KEYWORDS deprecated forwarder quality control
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments
    r
end
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.detectImpulses'))
    NonGravForces.fail("detectImpulses", "depMissing", "qc.detectImpulses forwards to timeSeriesPreparation.qc.detectImpulses and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:detectImpulses:deprecated", "NonGravForces.qc.detectImpulses is deprecated and forwards to timeSeriesPreparation.qc.detectImpulses -> call that directly; this name is removed in the release after the one that carries the forwarding");
[events, info] = timeSeriesPreparation.qc.detectImpulses(r, varargin{:});
end
