function varargout = detectJumps(varargin) %#ok<STOUT,INUSD>
%DETECTJUMPS Piecewise-constant/linear change-point detection (PELT + BIC). NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-31 (level L7). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-31.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:detectJumps:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-31; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.qc.detectJumps(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-31
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("detectJumps", "notImplemented", "detectJumps is specified in WP-31 -> not delivered in 2.0.0, see HANDOVER.md");
end
