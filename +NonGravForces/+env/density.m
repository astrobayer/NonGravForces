function varargout = density(varargin) %#ok<STOUT,INUSD>
%DENSITY Thermosphere density and composition dispatcher (exponential msis00 msis21 dtm2020 jb2008). NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-05 (level L2). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-05.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:density:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-05; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.env.density(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-05
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("density", "notImplemented", "density is specified in WP-05 -> not delivered in 2.0.0, see HANDOVER.md");
end
