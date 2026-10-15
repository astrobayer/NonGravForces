function [r, info] = residuals(y_meas, y_model, varargin)
%RESIDUALS [DEPRECATED, forwards to timeSeriesPreparation.qc.residuals] Residuals of measured minus modelled observations under a preliminary calibration, optionally whitened.
%
%   r = NonGravForces.qc.residuals(y_meas, y_model)
%   [r, info] = NonGravForces.qc.residuals(y_meas, y_model, 'Bias', b, 'Scale', S, 'Whitening', "highpass", 'Window', 601)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05): it ported this unit UNCHANGED, froze the port
%   and gates every later change against it with isequaln (X-16). This unit
%   now only forwards: it warns NonGravForces:residuals:deprecated and calls
%   timeSeriesPreparation.qc.residuals with the same arguments. It is deleted in
%   the release after the one that carries the forwarding. The target's errors
%   propagate unwrapped (timeSeriesPreparation:residuals:*), as every foreign error
%   does in this toolbox. Call the target directly.
%
% INPUT ARGUMENTS
%   Name      Size     Type     Units Description
%   ----      ----     ----     ----- -----------
%   y_meas    KxM      double   any   measured series, as timeSeriesPreparation.qc.residuals
%   y_model   KxM      double   any   modelled series, as timeSeriesPreparation.qc.residuals
%   varargin  1xN      cell     -     its name-value arguments, passed through unchanged
%
% NAME-VALUE ARGUMENTS
%   Bias, Scale, Whitening, Window: those of timeSeriesPreparation.qc.residuals,
%   with its defaults; the two argument blocks are identical (compared 2026-10-14).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type  Units Description
%   ----      ----  ----  ----- -----------
%   r         KxM   double  -  as timeSeriesPreparation.qc.residuals, unchanged
%   info      1x1   struct  -  as timeSeriesPreparation.qc.residuals, unchanged
%
% ACCURACY
%   None of its own: the target's numbers, bit for bit. tests/tNgfQc*.m run the
%   pre-move suite through this forwarder.
%
% ERRORS
%   NonGravForces:residuals:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:residuals:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:residuals:deprecated   on every call; call timeSeriesPreparation.qc.residuals
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05, X-16.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   help timeSeriesPreparation.qc.residuals
%
% See also timeSeriesPreparation.qc.residuals.
%
% KEYWORDS deprecated forwarder quality control
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments
    y_meas
    y_model
end
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.residuals'))
    NonGravForces.fail("residuals", "depMissing", "qc.residuals forwards to timeSeriesPreparation.qc.residuals and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:residuals:deprecated", "NonGravForces.qc.residuals is deprecated and forwards to timeSeriesPreparation.qc.residuals -> call that directly; this name is removed in the release after the one that carries the forwarding");
[r, info] = timeSeriesPreparation.qc.residuals(y_meas, y_model, varargin{:});
end
