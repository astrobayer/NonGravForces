function [mask, info] = hampel(x, varargin)
%HAMPEL [DEPRECATED, forwards to timeSeriesPreparation.qc.hampel] Robust point-outlier detector: moving median and scaled MAD, per column.
%
%   mask = NonGravForces.qc.hampel(x)
%   [mask, info] = NonGravForces.qc.hampel(x, 'Window', 61, 'Threshold', 4)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05): it ported this unit UNCHANGED, froze the port
%   and gates every later change against it with isequaln (X-16). This unit
%   now only forwards: it warns NonGravForces:hampel:deprecated and calls
%   timeSeriesPreparation.qc.hampel with the same arguments. It is deleted in
%   the release after the one that carries the forwarding. The target's errors
%   propagate unwrapped (timeSeriesPreparation:hampel:*), as every foreign error
%   does in this toolbox. Call the target directly.
%
% INPUT ARGUMENTS
%   Name      Size     Type     Units Description
%   ----      ----     ----     ----- -----------
%   x         KxM      double   any   samples down the columns, as timeSeriesPreparation.qc.hampel
%   varargin  1xN      cell     -     its name-value arguments, passed through unchanged
%
% NAME-VALUE ARGUMENTS
%   Window, Threshold, MinSigma: those of timeSeriesPreparation.qc.hampel,
%   with its defaults; the two argument blocks are identical (compared 2026-10-14).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units Description
%   ----      ----  ----    ----- -----------
%   mask      KxM   logical -  as timeSeriesPreparation.qc.hampel, unchanged
%   info      1x1   struct  -  as timeSeriesPreparation.qc.hampel, unchanged
%
% ACCURACY
%   None of its own: the target's numbers, bit for bit. tests/tNgfQc*.m run the
%   pre-move suite through this forwarder.
%
% ERRORS
%   NonGravForces:hampel:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:hampel:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:hampel:deprecated   on every call; call timeSeriesPreparation.qc.hampel
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05, X-16.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   help timeSeriesPreparation.qc.hampel
%
% See also timeSeriesPreparation.qc.hampel.
%
% KEYWORDS deprecated forwarder quality control
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments
    x
end
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.hampel'))
    NonGravForces.fail("hampel", "depMissing", "qc.hampel forwards to timeSeriesPreparation.qc.hampel and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:hampel:deprecated", "NonGravForces.qc.hampel is deprecated and forwards to timeSeriesPreparation.qc.hampel -> call that directly; this name is removed in the release after the one that carries the forwarding");
[mask, info] = timeSeriesPreparation.qc.hampel(x, varargin{:});
end
