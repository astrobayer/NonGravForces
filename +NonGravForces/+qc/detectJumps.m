function [events, info] = detectJumps(r, varargin)
%DETECTJUMPS [DEPRECATED, forwards to timeSeriesPreparation.qc.detectJumps] Change points in level or slope of whitened residuals by PELT with a BIC penalty; eclipse and mass priors are assigned.
%
%   [events, info] = NonGravForces.qc.detectJumps(r, 'Epoch', t)
%   [events, info] = NonGravForces.qc.detectJumps(r, 'Epoch', t, 'Model', "slope", 'Priors', ecl, 'Penalty', beta)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05): it ported this unit UNCHANGED, froze the port
%   and gates every later change against it with isequaln (X-16). This unit
%   now only forwards: it warns NonGravForces:detectJumps:deprecated and calls
%   timeSeriesPreparation.qc.detectJumps with the same arguments. It is deleted in
%   the release after the one that carries the forwarding. The target's errors
%   propagate unwrapped (timeSeriesPreparation:detectJumps:*), as every foreign error
%   does in this toolbox. Call the target directly.
%
% INPUT ARGUMENTS
%   Name      Size     Type     Units Description
%   ----      ----     ----     ----- -----------
%   r         Kx1|Kx3  double   any   residual series, as timeSeriesPreparation.qc.detectJumps
%   varargin  1xN      cell     -     its name-value arguments, passed through unchanged
%
% NAME-VALUE ARGUMENTS
%   Epoch, Model, Penalty, MinSegment, Priors, PriorWindow_s: those of timeSeriesPreparation.qc.detectJumps,
%   with its defaults; the two argument blocks are identical (compared 2026-10-14).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type  Units Description
%   ----      ----  ----  ----- -----------
%   events    Nx8   table   -  as timeSeriesPreparation.qc.detectJumps, unchanged
%   info      1x1   struct  -  as timeSeriesPreparation.qc.detectJumps, unchanged
%
% ACCURACY
%   None of its own: the target's numbers, bit for bit. tests/tNgfQc*.m run the
%   pre-move suite through this forwarder.
%
% ERRORS
%   NonGravForces:detectJumps:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:detectJumps:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:detectJumps:deprecated   on every call; call timeSeriesPreparation.qc.detectJumps
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05, X-16.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   help timeSeriesPreparation.qc.detectJumps
%
% See also timeSeriesPreparation.qc.detectJumps.
%
% KEYWORDS deprecated forwarder quality control
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments
    r
end
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.detectJumps'))
    NonGravForces.fail("detectJumps", "depMissing", "qc.detectJumps forwards to timeSeriesPreparation.qc.detectJumps and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:detectJumps:deprecated", "NonGravForces.qc.detectJumps is deprecated and forwards to timeSeriesPreparation.qc.detectJumps -> call that directly; this name is removed in the release after the one that carries the forwarding");
[events, info] = timeSeriesPreparation.qc.detectJumps(r, varargin{:});
end
