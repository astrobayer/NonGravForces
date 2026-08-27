function varargout = accommodation(varargin) %#ok<STOUT,INUSD>
%ACCOMMODATION Energy accommodation coefficient: fixed sesam langmuir. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-26 (level L4). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-26.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:accommodation:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-26; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.kernel.accommodation(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-26
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("accommodation", "notImplemented", "accommodation is specified in WP-26 -> not delivered in 2.0.0, see HANDOVER.md");
end
