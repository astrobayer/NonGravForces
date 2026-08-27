function varargout = calibrate(varargin)
%CALIBRATE Accelerometer calibration bridge to lsqvce (bias spline + 3x3 scale + jumps). NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-19 (level L8). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-19.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:calibrate:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-19; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.app.calibrate(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-19
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
varargout = cell(1, nargout);   % set so a caller with outputs is refused by fail, not by an unset return
NonGravForces.fail("calibrate", "notImplemented", "calibrate is specified in WP-19 -> not delivered in 2.0.0, see HANDOVER.md");
end
