function varargout = checkAttitudeSense(varargin) %#ok<STOUT,INUSD>
%CHECKATTITUDESENSE Physics test of quaternion sense (nadir and line-of-sight). NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-11 (level L3). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-11.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:checkAttitudeSense:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-11; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.geom.checkAttitudeSense(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-11
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("checkAttitudeSense", "notImplemented", "checkAttitudeSense is specified in WP-11 -> not delivered in 2.0.0, see HANDOVER.md");
end
