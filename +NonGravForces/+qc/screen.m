function [mask, events, info] = screen(y_meas, y_model, epoch, varargin)
%SCREEN [DEPRECATED, forwards to timeSeriesPreparation.qc.screen] Quality-control orchestrator: residuals -> outliers -> impulses -> jumps -> mask and event table.
%
%   [mask, events, info] = NonGravForces.qc.screen(y_meas, y_model, t)
%   [mask, events, info] = NonGravForces.qc.screen(y_meas, y_model, t, 'Bias', b, 'Scale', S, 'Response', crn, 'Thrusters', thr, 'Priors', ecl)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05): it ported this unit UNCHANGED, froze the port
%   and gates every later change against it with isequaln (X-16). This unit
%   now only forwards: it warns NonGravForces:screen:deprecated and calls
%   timeSeriesPreparation.qc.screen with the same arguments. It is deleted in
%   the release after the one that carries the forwarding. The target's errors
%   propagate unwrapped (timeSeriesPreparation:screen:*), as every foreign error
%   does in this toolbox. Call the target directly.
%
% INPUT ARGUMENTS
%   Name      Size     Type     Units Description
%   ----      ----     ----     ----- -----------
%   y_meas    KxM      double   any   measured series, as timeSeriesPreparation.qc.screen
%   y_model   KxM      double   any   modelled series, as timeSeriesPreparation.qc.screen
%   epoch     Kx1      datetime UTC   sample epochs, as timeSeriesPreparation.qc.screen
%   varargin  1xN      cell     -     its name-value arguments, passed through unchanged
%
% NAME-VALUE ARGUMENTS
%   Bias, Scale, Whitening, Window, HampelWindow, HampelThreshold, Response, ImpulseThreshold, Thrusters, Priors, JumpModel, JumpPenalty, JumpMinSegment: those of timeSeriesPreparation.qc.screen,
%   with its defaults; the two argument blocks are identical (compared 2026-10-14).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units Description
%   ----      ----  ----    ----- -----------
%   mask      KxM   logical -  as timeSeriesPreparation.qc.screen, unchanged
%   events    Nx8   table   -  as timeSeriesPreparation.qc.screen, unchanged
%   info      1x1   struct  -  as timeSeriesPreparation.qc.screen, unchanged
%
% ACCURACY
%   None of its own: the target's numbers, bit for bit. tests/tNgfQc*.m run the
%   pre-move suite through this forwarder.
%
% ERRORS
%   NonGravForces:screen:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:screen:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:screen:deprecated   on every call; call timeSeriesPreparation.qc.screen
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05, X-16.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   help timeSeriesPreparation.qc.screen
%
% See also timeSeriesPreparation.qc.screen.
%
% KEYWORDS deprecated forwarder quality control
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments
    y_meas
    y_model
    epoch
end
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.screen'))
    NonGravForces.fail("screen", "depMissing", "qc.screen forwards to timeSeriesPreparation.qc.screen and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:screen:deprecated", "NonGravForces.qc.screen is deprecated and forwards to timeSeriesPreparation.qc.screen -> call that directly; this name is removed in the release after the one that carries the forwarding");
[mask, events, info] = timeSeriesPreparation.qc.screen(y_meas, y_model, epoch, varargin{:});
end
