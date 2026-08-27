function varargout = bodyTorques(varargin) %#ok<STOUT,INUSD>
%BODYTORQUES Gravity-gradient and magnetic body torques. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-24 (level L4). This stub exists so that
%   the package layout is complete from the first commit and so that a caller
%   written against the future API fails LOUDLY with a stable identifier
%   instead of resolving to something else on the path (rule R4).
%
% INPUT ARGUMENTS
%   varargin   any   -   Accepted and ignored; the signature is fixed in WP-24.
%
% OUTPUT ARGUMENTS
%   (none)     Never returned; the function always raises.
%
% ERRORS
%   NonGravForces:bodyTorques:notImplemented   Always.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-24; doc/NonGravForces_Theory_and_Implementation.pdf.
%
% EXAMPLE
%   try, NonGravForces.kernel.bodyTorques(); catch ME, disp(ME.identifier); end
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-24
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
NonGravForces.fail("bodyTorques", "notImplemented", "bodyTorques is specified in WP-24 -> not delivered in 2.0.0, see HANDOVER.md");
end
