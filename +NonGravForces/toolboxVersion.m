function v = toolboxVersion(format)
%TOOLBOXVERSION Semantic version of NonGravForces - the single authority.
%
%   v = NonGravForces.toolboxVersion            % struct
%   s = NonGravForces.toolboxVersion('string')  % '2.0.0'
%
%   The patch number IS the revision number of HANDOVER.md (process rule P-7):
%   both are changed in the same edit, so a code file and the document that
%   specifies it can never silently disagree about which revision they belong
%   to.
%
% INPUT ARGUMENTS
%   Name     Size  Type   Units  Description
%   ----     ----  ----   -----  -----------
%   format   1x1   string -      "struct" (default) or "string".
%
% OUTPUT ARGUMENTS
%   Name     Size  Type   Units  Description
%   ----     ----  ----   -----  -----------
%   v        1x1   struct -      .major .minor .patch (uint16), .semantic
%                                (string), .date (string, ISO), .concept
%                                (uint8, 3 = this redesign), .model (string:
%                                authoring model). With "string": char '2.0.0'.
%
% ERRORS
%   MATLAB:validators:mustBeMember   Unknown format keyword.
%
% EXAMPLE
%   NonGravForces.toolboxVersion('string')
%
% See also NonGravForces.dependencyReport.
%
% KEYWORDS version semver
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    format (1,1) string {mustBeMember(format, ["struct","string"])} = "struct"
end
v = struct('major', uint16(2), 'minor', uint16(0), 'patch', uint16(0), ...
           'semantic', "2.0.0", 'date', "2026-08-27", 'concept', uint8(3), ...
           'model', "Claude Fable 5 (Anthropic)");
if format == "string"
    v = char(v.semantic);
end
end
