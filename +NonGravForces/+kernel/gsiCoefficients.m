function varargout = gsiCoefficients(varargin)
%GSICOEFFICIENTS Pressure and shear coefficients per panel: sentman cll maxwell schamberg. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-13 (level L4). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-13.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:gsiCoefficients:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-13; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.kernel.gsiCoefficients(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-13
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
varargout = cell(1, nargout);   % set so a caller with outputs is refused by fail, not by an unset return
NonGravForces.fail("gsiCoefficients", "notImplemented", "gsiCoefficients is specified in WP-13 -> not delivered yet, see HANDOVER.md");
end
