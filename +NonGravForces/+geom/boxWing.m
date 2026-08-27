function varargout = boxWing(varargin)
%BOXWING Generic box-wing macro model factory. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-09 (level L3). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-09.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:boxWing:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-09; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.geom.boxWing(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-09
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
varargout = cell(1, nargout);   % set so a caller with outputs is refused by fail, not by an unset return
NonGravForces.fail("boxWing", "notImplemented", "boxWing is specified in WP-09 -> not delivered in 2.0.0, see HANDOVER.md");
end
