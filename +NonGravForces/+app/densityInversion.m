function varargout = densityInversion(varargin) %#ok<STOUT,INUSD>
%DENSITYINVERSION Thermosphere density from calibrated accelerations (Doornbos). NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-32 (level L8). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-32.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:densityInversion:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-32; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.app.densityInversion(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-32
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("densityInversion", "notImplemented", "densityInversion is specified in WP-32 -> not delivered in 2.0.0, see HANDOVER.md");
end
